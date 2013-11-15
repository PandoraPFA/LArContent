/**
 *  @file   LArContent/src/LArTwoDSeed/SeedMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the particle seed algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArTwoDSeed/SeedMergingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedMergingAlgorithm::Run()
{
    ParticleSeedVector particleSeedVector;

    try
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

        this->GetParticleSeeds(pClusterList, particleSeedVector);
        this->MakeClusterMerges(particleSeedVector);
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "SeedMergingAlgorithm: exception " << statusCodeException.ToString() << std::endl;
    }

    for (ParticleSeedVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        delete *iter;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedMergingAlgorithm::GetParticleSeeds(const ClusterList *const pClusterList, ParticleSeedVector &particleSeedVector) const
{
    ParticleSeedVector localSeedVector;

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        localSeedVector.push_back(new ParticleSeed(*iter));
    }

    std::sort(localSeedVector.begin(), localSeedVector.end(), SeedMergingAlgorithm::SortByLayerSpan);

    while (true)
    {
        ParticleSeedVector::iterator primarySeedIter = localSeedVector.end();

        for (ParticleSeedVector::iterator iterP = localSeedVector.begin(), iterPEnd = localSeedVector.end(); iterP != iterPEnd; ++iterP)
        {
            if ((*iterP) != NULL)
            {
                primarySeedIter = iterP;
                break;
            }
        }

        if (localSeedVector.end() == primarySeedIter)
            break;

        ParticleSeed *pPrimarySeed = *primarySeedIter;
        particleSeedVector.push_back(pPrimarySeed);
        *primarySeedIter = NULL;

        ParticleSeedVector associatedSeeds;
        this->FindAssociatedSeeds(pPrimarySeed, localSeedVector, associatedSeeds);

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

void SeedMergingAlgorithm::FindAssociatedSeeds(ParticleSeed *const pParticleSeed, ParticleSeedVector &candidateSeeds,
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

bool SeedMergingAlgorithm::AreParticleSeedsAssociated(const ParticleSeed *const pAssociatedSeed, const ParticleSeed *const pCandidateSeed) const
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

bool SeedMergingAlgorithm::AreSeedClustersAssociated(Cluster *const pClusterI, Cluster *const pClusterJ) const
{
    if (!LArVertexHelper::DoesCurrentVertexExist())
        return false;

    const CartesianVector &eventVertex(LArVertexHelper::GetCurrentVertex());
    const LArPointingCluster pointingClusterI(pClusterI);
    const LArPointingCluster pointingClusterJ(pClusterJ);

    const bool isVertexAssociatedI(LArPointingClusterHelper::IsNode(eventVertex, pointingClusterI.GetInnerVertex()) ||
        LArPointingClusterHelper::IsNode(eventVertex, pointingClusterI.GetOuterVertex()) ||
        LArPointingClusterHelper::IsEmission(eventVertex, pointingClusterI.GetInnerVertex()) ||
        LArPointingClusterHelper::IsEmission(eventVertex, pointingClusterI.GetOuterVertex()) );

    const bool isVertexAssociatedJ(LArPointingClusterHelper::IsNode(eventVertex, pointingClusterJ.GetInnerVertex()) ||
        LArPointingClusterHelper::IsNode(eventVertex, pointingClusterJ.GetOuterVertex()) ||
        LArPointingClusterHelper::IsEmission(eventVertex, pointingClusterJ.GetInnerVertex()) ||
        LArPointingClusterHelper::IsEmission(eventVertex, pointingClusterJ.GetOuterVertex()) );

    if (isVertexAssociatedI && isVertexAssociatedJ)
        return false;

    // Apply track-like selection criteria 
    const bool isTrackLikeI = (LArClusterHelper::LArTrackWidth(pClusterI) < 0.15f);
    const bool isTrackLikeJ = (LArClusterHelper::LArTrackWidth(pClusterJ) < 0.15f);

    if (isTrackLikeI || isTrackLikeJ)
        return false;

    // Direction measurements
    const bool currentVertexExists(LArVertexHelper::DoesCurrentVertexExist());

    const bool isForwardI(currentVertexExists ? LArVertexHelper::IsForwardInZ(pClusterI) : false);
    const bool isBackwardI(currentVertexExists ? LArVertexHelper::IsBackwardInZ(pClusterI) : false);
    const bool checkForwardI(isForwardI || !isBackwardI);
    const bool checkBackwardI(isBackwardI || !isForwardI);

    const bool isForwardJ(currentVertexExists ? LArVertexHelper::IsForwardInZ(pClusterJ) : false);
    const bool isBackwardJ(currentVertexExists ? LArVertexHelper::IsBackwardInZ(pClusterJ) : false);
    const bool checkForwardJ(isForwardJ || !isBackwardJ);
    const bool checkBackwardJ(isBackwardJ || !isForwardJ);

    const bool consistentLongitudinalDirection = ((checkForwardI && checkForwardJ) || (checkBackwardI && checkBackwardJ));

    if (!consistentLongitudinalDirection)
        return false;

    // Calculate proximity variables and look for overlapping showers
    const float rOuterI(LArClusterHelper::GetClosestDistance(pClusterI->GetCentroid(pClusterI->GetOuterPseudoLayer()), pClusterJ));
    const float rOuterJ(LArClusterHelper::GetClosestDistance(pClusterJ->GetCentroid(pClusterJ->GetOuterPseudoLayer()), pClusterI));
    const float rInnerI(LArClusterHelper::GetClosestDistance(pClusterI->GetCentroid(pClusterI->GetInnerPseudoLayer()), pClusterJ));
    const float rInnerJ(LArClusterHelper::GetClosestDistance(pClusterJ->GetCentroid(pClusterJ->GetInnerPseudoLayer()), pClusterI));

    if ((rInnerJ < 5. && rOuterJ < 5.) && (!checkForwardI || rInnerI > 5.) && (!checkBackwardI || rOuterI > 5.))
        return true;

    if ((rInnerI < 5. && rOuterI < 5.) && (!checkForwardJ || rInnerJ > 5.) && (!checkBackwardJ || rOuterJ > 5.))
        return true;

    if ((rInnerI < 5. && rOuterJ < 5.) && (!checkForwardJ || rInnerJ > 5.) && (!checkBackwardI || rOuterI > 5.))
        return true;

    if ((rInnerJ < 5. && rOuterI < 5.) && (!checkForwardI || rInnerI > 5.) && (!checkBackwardJ || rOuterJ > 5.))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedMergingAlgorithm::MakeClusterMerges(const ParticleSeedVector &particleSeedVector) const
{
    for (ParticleSeedVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        const ParticleSeed *pParticleSeed = *iter;
        const ClusterList &particleSeedClusterList(pParticleSeed->GetClusterList());

        if (particleSeedClusterList.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Cluster *pParentCluster = *(particleSeedClusterList.begin());

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

bool SeedMergingAlgorithm::SortByLayerSpan(const ParticleSeed *const pLhs, const ParticleSeed *const pRhs)
{
    unsigned int minLayerLhs(std::numeric_limits<unsigned int>::max()), minLayerRhs(std::numeric_limits<unsigned int>::max()), maxLayerLhs(0), maxLayerRhs(0);
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

    unsigned int layerSpanLhs((maxLayerLhs>minLayerLhs) ? maxLayerLhs - minLayerLhs : 0);
    unsigned int layerSpanRhs((maxLayerRhs>minLayerRhs) ? maxLayerRhs - minLayerRhs : 0);

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return (energyLhs > energyRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SeedMergingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
