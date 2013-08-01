/**
 *  @file   LArContent/src/LArTwoDSeed/SeedFindingAlgorithm.cc
 * 
 *  @brief  Implementation of the particle seed algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArTwoDSeed/SeedFindingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedFindingAlgorithm::Run()
{
    // Input lists
    const ClusterList *pVertexSeedClusterList = NULL;
    const ClusterList *pNonSeedClusterList = NULL;

    if (STATUS_CODE_SUCCESS == PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pVertexSeedClusterList))
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this,
            m_nonSeedClusterListName, pNonSeedClusterList));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pNonSeedClusterList));
    }

    // Particle seed formation
    ClusterList vertexSeedClusters((NULL != pVertexSeedClusterList) ? ClusterList(*pVertexSeedClusterList) : ClusterList());
    ClusterList nonSeedClusters((NULL != pNonSeedClusterList) ? ClusterList(*pNonSeedClusterList) : ClusterList());

    ParticleSeedVector particleSeedVector;
    this->GetParticleSeeds(vertexSeedClusters, nonSeedClusters, particleSeedVector);

    if (particleSeedVector.empty())
        return STATUS_CODE_SUCCESS;

    // List manipulation - ensure same behaviour if there are or are not separate vertex clusters
    if (NULL != pVertexSeedClusterList)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListName, vertexSeedClusters));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentClusterList(*this, m_nonSeedClusterListName));
    }

    // Make cluster merges and write out final lists
    ClusterList finalSeedClusters;
    this->MakeClusterMerges(particleSeedVector, finalSeedClusters);

    if ((NULL == pVertexSeedClusterList) && (NULL != pNonSeedClusterList))
    {
        ClusterList finalNonSeedClusters(*pNonSeedClusterList);

        for (ClusterList::const_iterator iter = finalSeedClusters.begin(), iterEnd = finalSeedClusters.end(); iter != iterEnd; ++iter)
            finalNonSeedClusters.erase(*iter);

        if (!finalNonSeedClusters.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListName, finalNonSeedClusters));
    }

    if (!finalSeedClusters.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName, finalSeedClusters));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentClusterList(*this, m_seedClusterListName));
    }

    for (ParticleSeedVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
        delete *iter;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedFindingAlgorithm::GetParticleSeeds(const ClusterList &vertexSeedClusters, const ClusterList &nonSeedClusters, ParticleSeedVector &particleSeedVector) const
{
    // Create vertex-based particle seeds, using clusters from earlier algorithm
    ParticleSeedVector vertexSeedVector;

    for (ClusterList::const_iterator iter = vertexSeedClusters.begin(), iterEnd = vertexSeedClusters.end(); iter != iterEnd; ++iter)
    {
        vertexSeedVector.push_back(new ParticleSeed(*iter));
    }

    // Create length-based particle seeds, using sliding window length cut
    ParticleSeedVector lengthSeedVector;
    int nIter(vertexSeedClusters.size());

    ClusterVector clusterVector(nonSeedClusters.begin(), nonSeedClusters.end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        int lengthCut((nIter < m_finalChangeIter) ? m_initialLengthCut : m_finalLengthCut);

        if ((nIter > m_initialChangeIter) && (nIter < m_finalChangeIter))
            lengthCut = m_initialLengthCut + (nIter - m_initialChangeIter) * (m_finalLengthCut - m_initialLengthCut) / (m_finalChangeIter - m_initialChangeIter);

        if ((*iter)->GetOrderedCaloHitList().size() >= lengthCut)
            lengthSeedVector.push_back(new ParticleSeed(*iter));

        ++nIter;
    }

    // Associate vertex and length-based seeds
    std::sort(vertexSeedVector.begin(), vertexSeedVector.end(), SeedFindingAlgorithm::SortByLayerSpan);
    std::sort(lengthSeedVector.begin(), lengthSeedVector.end(), SeedFindingAlgorithm::SortByLayerSpan);

    for (ParticleSeedVector::const_iterator iterP = vertexSeedVector.begin(), iterPEnd = vertexSeedVector.end(); iterP != iterPEnd; ++iterP)
    {
        ParticleSeed *pPrimarySeed = *iterP;
        particleSeedVector.push_back(pPrimarySeed);

        ParticleSeedVector associatedSeeds;
        this->FindAssociatedSeeds(pPrimarySeed, lengthSeedVector, associatedSeeds);

        for (ParticleSeedVector::iterator iterA = associatedSeeds.begin(), iterAEnd = associatedSeeds.end(); iterA != iterAEnd; ++iterA)
        {
            if ((*iterA) != pPrimarySeed)
            {
                pPrimarySeed->AddClusterList((*iterA)->GetClusterList());
                delete *iterA;
            }
        }
    }

    // Associations between remaining length-based seeds
    while (true)
    {
        ParticleSeedVector::iterator primarySeedIter = lengthSeedVector.end();

        for (ParticleSeedVector::iterator iterP = lengthSeedVector.begin(), iterPEnd = lengthSeedVector.end(); iterP != iterPEnd; ++iterP)
        {
            if ((*iterP) != NULL)
            {
                primarySeedIter = iterP;
                break;
            }
        }

        if (lengthSeedVector.end() == primarySeedIter)
            break;

        ParticleSeed *pPrimarySeed = *primarySeedIter;
        particleSeedVector.push_back(pPrimarySeed);
        *primarySeedIter = NULL;

        ParticleSeedVector associatedSeeds;
        this->FindAssociatedSeeds(pPrimarySeed, lengthSeedVector, associatedSeeds);

        for (ParticleSeedVector::iterator iterA = associatedSeeds.begin(), iterAEnd = associatedSeeds.end(); iterA != iterAEnd; ++iterA)
        {
            if ((*iterA) != pPrimarySeed)
            {
                pPrimarySeed->AddClusterList((*iterA)->GetClusterList());
                delete *iterA;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedFindingAlgorithm::FindAssociatedSeeds(ParticleSeed *const pParticleSeed, ParticleSeedVector &candidateSeeds,
    ParticleSeedVector &associatedSeeds) const
{
    ParticleSeedVector currentSeedAssociations, newSeedAssociations;
    currentSeedAssociations.push_back(pParticleSeed);

    while (!currentSeedAssociations.empty())
    {
        for (ParticleSeedVector::iterator iterI = candidateSeeds.begin(), iterIEnd = candidateSeeds.end(); iterI != iterIEnd; ++iterI)
        {
            ParticleSeed *pCandidateParticleSeed = *iterI;

            if (NULL == pCandidateParticleSeed)
                continue;

            for (ParticleSeedVector::iterator iterJ = currentSeedAssociations.begin(), iterJEnd = currentSeedAssociations.end(); iterJ != iterJEnd; ++iterJ)
            {
                ParticleSeed *pAssociatedParticleSeed = *iterJ;

                if (!this->AreParticleSeedsAssociated(pAssociatedParticleSeed, pCandidateParticleSeed))
                    continue;

                newSeedAssociations.push_back(pCandidateParticleSeed);
                *iterI = NULL;
                break;
            }
        }

        associatedSeeds.insert(associatedSeeds.end(), currentSeedAssociations.begin(), currentSeedAssociations.end());
        currentSeedAssociations = newSeedAssociations;
        newSeedAssociations.clear();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedFindingAlgorithm::AreParticleSeedsAssociated(const ParticleSeed *const pAssociatedSeed, const ParticleSeed *const pCandidateSeed) const
{
    const ClusterList &associatedClusterList(pAssociatedSeed->GetClusterList());
    const ClusterList &candidateClusterList(pCandidateSeed->GetClusterList());

    for (ClusterList::const_iterator iterI = associatedClusterList.begin(), iterIEnd = associatedClusterList.end(); iterI != iterIEnd; ++iterI)
    {
        for (ClusterList::const_iterator iterJ = candidateClusterList.begin(), iterJEnd = candidateClusterList.end(); iterJ != iterJEnd; ++iterJ)
        {
            if (this->AreSeedClustersAssociated(*iterI, *iterJ))
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedFindingAlgorithm::AreSeedClustersAssociated(const Cluster *const pClusterI, const Cluster *const pClusterJ) const
{
    // Direction measurements
    const bool currentVertexExists(LArVertexHelper::DoesCurrentVertexExist());

    const bool isForwardI = (currentVertexExists) ? LArVertexHelper::IsForwardInZ(pClusterI) : false;
    const bool isBackwardI = (currentVertexExists) ? LArVertexHelper::IsBackwardInZ(pClusterI) : false;
    const bool checkForwardI = (isForwardI || !isBackwardI);
    const bool checkBackwardI = (isBackwardI || !isForwardI);

    const bool isForwardJ = (currentVertexExists) ? LArVertexHelper::IsForwardInZ(pClusterJ) : false;
    const bool isBackwardJ = (currentVertexExists) ? LArVertexHelper::IsBackwardInZ(pClusterJ) : false;
    const bool checkForwardJ = (isForwardJ || !isBackwardJ);
    const bool checkBackwardJ = (isBackwardJ || !isForwardJ);

    // Consistent Longitudinal Directions
    const bool consistentForwardDirection = (checkForwardI && checkForwardJ);
    const bool consistentBackwardDirection = (checkBackwardI && checkBackwardJ);
    const bool consistentLongitudinalDirection = (consistentForwardDirection || consistentBackwardDirection);

    // Apply track-like selection criteria 
    const bool isTrackLikeI = (LArClusterHelper::LArTrackWidth(pClusterI) < 0.15);
    const bool isTrackLikeJ = (LArClusterHelper::LArTrackWidth(pClusterJ) < 0.15);
    const bool useTrackLikeCuts = (isTrackLikeI || isTrackLikeJ);

    // Calculate proximity variables
    const float rOuterI(LArClusterHelper::GetClosestDistance(pClusterI->GetCentroid(pClusterI->GetOuterPseudoLayer()), pClusterJ));
    const float rOuterJ(LArClusterHelper::GetClosestDistance(pClusterJ->GetCentroid(pClusterJ->GetOuterPseudoLayer()), pClusterI));
    const float rInnerI(LArClusterHelper::GetClosestDistance(pClusterI->GetCentroid(pClusterI->GetInnerPseudoLayer()), pClusterJ));
    const float rInnerJ(LArClusterHelper::GetClosestDistance(pClusterJ->GetCentroid(pClusterJ->GetInnerPseudoLayer()), pClusterI));

    // Association check 1, look for overlapping clusters
    if (consistentLongitudinalDirection)
    {
        if (useTrackLikeCuts)
        {
            // look for overlapping tracks
            if ((rInnerJ < 2.5 && rOuterJ < 2.5) && (!checkForwardI || rInnerI > 10.) && (!checkBackwardI || rOuterI > 10.))
                return true;

            if ((rInnerI < 2.5 && rOuterI < 2.5) && (!checkForwardJ || rInnerJ > 10.) && (!checkBackwardJ || rOuterJ > 10.))
                return true;

            if ((rInnerI < 2.5 && rOuterJ < 2.5) && (!checkForwardJ || rInnerJ > 10.) && (!checkBackwardI || rOuterI > 10.))
                return true;

            if ((rInnerJ < 2.5 && rOuterI < 2.5) && (!checkForwardI || rInnerI > 10.) && (!checkBackwardJ || rOuterJ > 10.))
                return true;
        }
        else
        {
            // look for overlapping showers
            if ((rInnerJ < 5. && rOuterJ < 5.) && (!checkForwardI || rInnerI > 5.) && (!checkBackwardI || rOuterI > 5.))
                return true;

            if ((rInnerI < 5. && rOuterI < 5.) && (!checkForwardJ || rInnerJ > 5.) && (!checkBackwardJ || rOuterJ > 5.))
                return true;

            if ((rInnerI < 5. && rOuterJ < 5.) && (!checkForwardJ || rInnerJ > 5.) && (!checkBackwardI || rOuterI > 5.))
                return true;

            if ((rInnerJ < 5. && rOuterI < 5.) && (!checkForwardI || rInnerI > 5.) && (!checkBackwardJ || rOuterJ > 5.))
                return true;
        }
    }

    // Association check 2, look for branching clusters (where J is a branch of I)
    if ((!checkForwardI || rInnerI > 20.) && (!checkBackwardI || rOuterI > 20.) && ((checkForwardJ && rInnerJ < 2.5) || (checkBackwardJ && rOuterJ < 2.5)))
        return true;

    // Association check 2, look for branching clusters (where I is a branch of J)
    if ((!checkForwardJ || rInnerJ > 20.) && (!checkBackwardJ || rOuterJ > 20.) && ((checkForwardI && rInnerI < 2.5) || (checkBackwardI && rOuterI < 2.5)))
        return true;

    // Association check 3, association through pointing

    // Bail out if clusters are well separated
    if (std::min(std::min(rInnerJ, rOuterJ), std::min(rInnerI, rOuterI) ) > 50.)
        return false;

    // Remove overlaps before checking pointing information
    if ((pClusterJ->GetOuterPseudoLayer() > pClusterI->GetInnerPseudoLayer()) && (pClusterI->GetOuterPseudoLayer() > pClusterJ->GetInnerPseudoLayer()))
        return false;

    // Calculate distances and angles
    static const unsigned int m_seedFitLayerWindow(30);

    ClusterHelper::ClusterFitResult innerFitI, innerFitJ, outerFitI, outerFitJ;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pClusterI, m_seedFitLayerWindow, innerFitI));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pClusterJ, m_seedFitLayerWindow, innerFitJ));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pClusterI, m_seedFitLayerWindow, outerFitI));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pClusterJ, m_seedFitLayerWindow, outerFitJ));

    const CartesianVector innerDirectionI(innerFitI.GetDirection());
    const CartesianVector outerDirectionI(outerFitI.GetDirection());
    const CartesianVector innerDirectionJ(innerFitJ.GetDirection());
    const CartesianVector outerDirectionJ(outerFitJ.GetDirection());

    const CartesianVector innerCentroidI(pClusterI->GetCentroid(pClusterI->GetInnerPseudoLayer()));
    const CartesianVector outerCentroidI(pClusterI->GetCentroid(pClusterI->GetOuterPseudoLayer()));
    const CartesianVector innerCentroidJ(pClusterJ->GetCentroid(pClusterJ->GetInnerPseudoLayer()));
    const CartesianVector outerCentroidJ(pClusterJ->GetCentroid(pClusterJ->GetOuterPseudoLayer()));

    // Pointing check 1, longitudinal distances
    if (pClusterI->GetInnerPseudoLayer() > pClusterJ->GetOuterPseudoLayer())
    {
        if (pClusterI->GetInnerPseudoLayer() - pClusterJ->GetOuterPseudoLayer() > 100)
            return false;

        if (std::max(outerDirectionI.GetDotProduct(innerCentroidI - outerCentroidJ), innerFitJ.GetDirection().GetDotProduct(innerCentroidI - outerCentroidJ)) > 30.)
            return false;
    }

    if (pClusterJ->GetInnerPseudoLayer() > pClusterI->GetOuterPseudoLayer())
    {
        if (pClusterJ->GetInnerPseudoLayer() - pClusterI->GetOuterPseudoLayer() > 100)
            return false;

        if (std::max(outerDirectionJ.GetDotProduct(innerCentroidJ - outerCentroidI), innerFitI.GetDirection().GetDotProduct(innerCentroidJ - outerCentroidI)) > 30.)
            return false;
    }

    // Pointing check 2, pointing angle and transverse distances
    // check 2(a): innerI->innerJ
    if (consistentForwardDirection)
    {
        if (innerDirectionI.GetCosOpeningAngle(innerDirectionJ) > 0.985)
        {
            if ((pClusterI->GetInnerPseudoLayer() < pClusterJ->GetInnerPseudoLayer()) && (innerDirectionI.GetDotProduct(innerCentroidJ - innerCentroidI) < 50.) && 
                ((innerDirectionI.GetCrossProduct(innerCentroidJ - innerCentroidI)).GetMagnitudeSquared() < 2.5 * 2.5))
            {
                return true;
            }

            if ((pClusterJ->GetInnerPseudoLayer() < pClusterI->GetInnerPseudoLayer()) && (innerDirectionJ.GetDotProduct(innerCentroidI - innerCentroidJ) < 50.) &&
                ((innerDirectionJ.GetCrossProduct(innerCentroidI - innerCentroidJ)).GetMagnitudeSquared() < 2.5 * 2.5))
            {
                return true;
            }
        }
    }

    // check 2(b): outerI->outerJ
    if (consistentBackwardDirection)
    {
        if (outerDirectionI.GetCosOpeningAngle(outerDirectionJ) > 0.985)
        {
            if ((pClusterI->GetOuterPseudoLayer() > pClusterJ->GetOuterPseudoLayer()) && (outerDirectionI.GetDotProduct(outerCentroidI - outerCentroidJ) < 50.) &&
                ((outerDirectionI.GetCrossProduct(outerCentroidI - outerCentroidJ)).GetMagnitudeSquared() < 2.5 * 2.5))
            {
                return true;
            }

            if ((pClusterJ->GetOuterPseudoLayer() > pClusterI->GetOuterPseudoLayer()) && (outerDirectionJ.GetDotProduct(outerCentroidJ - outerCentroidI) < 50.) &&
                ((outerDirectionJ.GetCrossProduct(outerCentroidJ - outerCentroidI)).GetMagnitudeSquared() < 2.5 * 2.5))
            {
                return true;
            }
        }
    }

    // check 2(c): outerJ->innerI
    if (consistentLongitudinalDirection)
    {
        if ((pClusterJ->GetOuterPseudoLayer() < pClusterI->GetInnerPseudoLayer()) && (innerDirectionI.GetCosOpeningAngle(outerDirectionJ) > 0.985))
        {
            if (((innerDirectionI.GetCrossProduct(innerCentroidI - outerCentroidJ)).GetMagnitudeSquared() < 2.5 * 2.5) ||
                ((outerDirectionJ.GetCrossProduct(innerCentroidI - outerCentroidJ)).GetMagnitudeSquared() < 2.5 * 2.5))
            {
                return true;
            }
        }
    }

    // check 2(d): outerI->innerJ
    if (consistentLongitudinalDirection)
    {
        if ((pClusterI->GetOuterPseudoLayer() < pClusterJ->GetInnerPseudoLayer()) && (outerDirectionI.GetCosOpeningAngle(innerDirectionJ) > 0.985))
        {
            if (((outerDirectionI.GetCrossProduct(innerCentroidJ - outerCentroidI)).GetMagnitudeSquared() < 2.5 * 2.5) ||
                ((innerDirectionJ.GetCrossProduct(innerCentroidJ - outerCentroidI)).GetMagnitudeSquared() < 2.5 * 2.5))
            {
                return true;
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedFindingAlgorithm::MakeClusterMerges(const ParticleSeedVector &particleSeedVector, pandora::ClusterList &finalClusterList) const
{
    for (ParticleSeedVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        const ParticleSeed *pParticleSeed = *iter;
        const ClusterList &particleSeedClusterList(pParticleSeed->GetClusterList());

        if (particleSeedClusterList.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Cluster *pParentCluster = *(particleSeedClusterList.begin());
        finalClusterList.insert(pParentCluster);

        for (ClusterList::const_iterator iterM = particleSeedClusterList.begin(), iterMEnd = particleSeedClusterList.end(); iterM != iterMEnd; ++iterM)
        {
            Cluster *pDaughterCluster = *iterM;

            if (pDaughterCluster == pParentCluster)
                continue;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedFindingAlgorithm::SortByLayerSpan(const ParticleSeed *const pLhs, const ParticleSeed *const pRhs)
{
    int minLayerLhs(std::numeric_limits<int>::max()), minLayerRhs(std::numeric_limits<int>::max()), maxLayerLhs(0), maxLayerRhs(0);
    float energyLhs(0.f), energyRhs(0.f);

    for (ClusterList::const_iterator iter = pLhs->GetClusterList().begin(), iterEnd = pLhs->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->GetOuterPseudoLayer() > maxLayerLhs) {maxLayerLhs = (*iter)->GetOuterPseudoLayer();}
        if ((*iter)->GetInnerPseudoLayer() < minLayerLhs) {minLayerLhs = (*iter)->GetInnerPseudoLayer();}
        energyLhs += (*iter)->GetHadronicEnergy();
    }

    for (ClusterList::const_iterator iter = pRhs->GetClusterList().begin(), iterEnd = pRhs->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->GetOuterPseudoLayer() > maxLayerRhs) {maxLayerRhs = (*iter)->GetOuterPseudoLayer();}
        if ((*iter)->GetInnerPseudoLayer() < minLayerRhs) {minLayerRhs = (*iter)->GetInnerPseudoLayer();}
        energyRhs += (*iter)->GetHadronicEnergy();
    }

    if ((maxLayerLhs - minLayerLhs) != (maxLayerRhs - minLayerRhs))
        return ((maxLayerLhs - minLayerLhs) > (maxLayerRhs - minLayerRhs));

    return (energyLhs > energyRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SeedFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    m_initialLengthCut = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "InitialLengthCut", m_initialLengthCut));

    m_finalLengthCut = 50;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "FinalLengthCut", m_finalLengthCut));

    m_initialChangeIter = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "InitialChangeIter", m_initialChangeIter));

    m_finalChangeIter = 6;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "FinalChangeIter", m_finalChangeIter));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
