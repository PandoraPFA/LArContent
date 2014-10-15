/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm() :
    m_minCaloHitsPerCluster(5),
    m_nearbyTrackDistance(1.f),
    m_nearbyClusterDistance(2.5f),
    m_remoteClusterDistance(10.f),
    m_useExistingPfosAsSeeds(true),
    m_useExistingPfosAsBranches(false),
    m_directionTanAngle(1.732f),
    m_directionApexShift(0.333f),
    m_useMCFigureOfMerit(false),
    m_useFirstImprovedSeed(false),
    m_shouldRemoveShowerPfos(true),
    m_showerLikeNBranches(5),
    m_showerLikeCaloHitRatio(2.f),
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexLongitudinalDistance(20.f),
    m_maxVertexTransverseDistance(1.5f),
    m_vertexAngularAllowance(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    PfoList pfoList;
    this->GetInputPfoList(pfoList);
    ClusterInfoMap nCaloHitsPerCluster, nBranchesPerCluster;

    for (StringVector::const_iterator listIter = m_inputClusterListNames.begin(), listIterEnd = m_inputClusterListNames.end(); listIter != listIterEnd; ++listIter)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *listIter, pClusterList));

        ClusterList usedClusters;
        Cluster *pSeedCluster(NULL);

        while (this->GetNextSeedCandidate(pClusterList, usedClusters, pSeedCluster))
        {
            SeedAssociationList seedAssociationList;
            this->GetSeedAssociationList(ClusterVector(1, pSeedCluster), pClusterList, seedAssociationList);

            if (seedAssociationList.size() != 1)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            SeedAssociationList finalSeedAssociationList;
            this->CheckSeedAssociationList(seedAssociationList.begin(), finalSeedAssociationList);

            for (SeedAssociationList::const_iterator iter = finalSeedAssociationList.begin(), iterEnd = finalSeedAssociationList.end(); iter != iterEnd; ++iter)
            {
                Cluster *pParentCluster(iter->first);
                const ClusterVector &branchClusters(iter->second);

                this->StoreNCaloHitsPerCluster(pParentCluster, nCaloHitsPerCluster);
                this->StoreNBranchesPerCluster(pParentCluster, branchClusters, nBranchesPerCluster);
                this->ProcessBranchClusters(pParentCluster, branchClusters, *listIter, pfoList);

                usedClusters.insert(pParentCluster);
                usedClusters.insert(branchClusters.begin(), branchClusters.end());
            }
        }
    }

    if (m_shouldRemoveShowerPfos)
        this->RemoveShowerPfos(nCaloHitsPerCluster, nBranchesPerCluster, pfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::GetNextSeedCandidate(const ClusterList *const pClusterList, const ClusterList &usedClusters,
    Cluster *&pSeedCluster) const
{
    pSeedCluster = NULL;

    ClusterVector clusterVector;
    clusterVector.insert(clusterVector.end(), pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), ClusterCharacterisationAlgorithm::SortClusters);

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (usedClusters.count(pCluster))
            continue;

        if (!m_useExistingPfosAsSeeds && !pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        pSeedCluster = pCluster;
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::GetSeedAssociationList(const ClusterVector &particleSeedVector, const ClusterList *const pClusterList,
    SeedAssociationList &seedAssociationList) const
{
    ClusterVector candidateClusters;
    ClusterList seedClusters(particleSeedVector.begin(), particleSeedVector.end());

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCandidateCluster = *iter;

        if (seedClusters.count(pCandidateCluster))
            continue;

        if (!m_useExistingPfosAsBranches && !pCandidateCluster->IsAvailable())
            continue;

        if (pCandidateCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        candidateClusters.push_back(pCandidateCluster);
    }

    std::sort(candidateClusters.begin(), candidateClusters.end(), ClusterCharacterisationAlgorithm::SortClusters);
    ClusterUsageMap forwardUsageMap, backwardUsageMap;

    for (ClusterVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        this->FindAssociatedClusters(*iter, candidateClusters, forwardUsageMap, backwardUsageMap);
    }

    this->IdentifyClusterMerges(particleSeedVector, backwardUsageMap, seedAssociationList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::CheckSeedAssociationList(SeedAssociationList::const_iterator seedIter, SeedAssociationList &finalSeedAssociationList) const
{
    ClusterList usedClusters, availableClusters;
    usedClusters.insert(seedIter->first);
    availableClusters.insert(seedIter->second.begin(), seedIter->second.end());

    SeedAssociationList seedAssociationList;
    seedAssociationList.insert(SeedAssociationList::value_type(seedIter->first, seedIter->second));
    const float originalFigureOfMerit(this->GetFigureOfMerit(seedAssociationList));
    float bestFigureOfMerit(originalFigureOfMerit);

    Cluster *pTrialSeedCluster = NULL;
    bool betterConfigurationFound(false);
    SeedAssociationList bestSeedAssociationList;

    while (this->GetNextSeedCandidate(&availableClusters, usedClusters, pTrialSeedCluster))
    {
        usedClusters.insert(pTrialSeedCluster);

        ClusterVector trialParticleSeedVector;
        trialParticleSeedVector.push_back(seedIter->first);
        trialParticleSeedVector.push_back(pTrialSeedCluster);

        SeedAssociationList trialSeedAssociationList;
        this->GetSeedAssociationList(trialParticleSeedVector, &availableClusters, trialSeedAssociationList);
        const float trialFigureOfMerit(this->GetFigureOfMerit(trialSeedAssociationList));

        if (trialFigureOfMerit > bestFigureOfMerit)
        {
            betterConfigurationFound = true;
            bestFigureOfMerit = trialFigureOfMerit;
            bestSeedAssociationList = trialSeedAssociationList;

            if (m_useFirstImprovedSeed)
                break;
        }
    }

    if (betterConfigurationFound)
    {
        for (SeedAssociationList::const_iterator iter = bestSeedAssociationList.begin(), iterEnd = bestSeedAssociationList.end(); iter != iterEnd; ++iter)
        {
            this->CheckSeedAssociationList(iter, finalSeedAssociationList);
        }
    }
    else
    {
        finalSeedAssociationList.insert(SeedAssociationList::value_type(seedIter->first, seedIter->second));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::ProcessBranchClusters(Cluster *const pParentCluster, const ClusterVector &branchClusters, const std::string &listName,
    PfoList &pfoList) const
{
    for (ClusterVector::const_iterator iter = branchClusters.begin(), iterEnd = branchClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pBranchCluster(*iter);

        if (!pBranchCluster->IsAvailable() && m_shouldRemoveShowerPfos)
        {
            Pfo *pTargetPfo(NULL);
            this->FindTargetPfo(pBranchCluster, pfoList, pTargetPfo);
            pfoList.erase(pTargetPfo);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pTargetPfo));
        }

        if (pBranchCluster->IsAvailable())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pBranchCluster, listName, listName));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterCharacterisationAlgorithm::AssociationType ClusterCharacterisationAlgorithm::AreClustersAssociated(const Cluster *const pClusterSeed, const Cluster *const pCluster) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *pVertex((pVertexList->size() == 1) ? *(pVertexList->begin()) : NULL);

    // Direction of seed cluster
    const LArVertexHelper::ClusterDirection seedDirection((NULL == pVertex) ? LArVertexHelper::DIRECTION_UNKNOWN :
        LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex, pClusterSeed, m_directionTanAngle, m_directionApexShift));
    const bool checkSeedForward(seedDirection != LArVertexHelper::DIRECTION_BACKWARD_IN_Z);
    const bool checkSeedBackward(seedDirection != LArVertexHelper::DIRECTION_FORWARD_IN_Z);

    // Direction of candidate cluster
    const LArVertexHelper::ClusterDirection candidateDirection((NULL == pVertex) ? LArVertexHelper::DIRECTION_UNKNOWN :
        LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex, pCluster, m_directionTanAngle, m_directionApexShift));
    const bool checkCandidateForward(candidateDirection != LArVertexHelper::DIRECTION_BACKWARD_IN_Z);
    const bool checkCandidateBackward(candidateDirection != LArVertexHelper::DIRECTION_FORWARD_IN_Z);

    // Calculate distances of association
    const float sOuter(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetOuterPseudoLayer()), pCluster));
    const float cOuter(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()), pClusterSeed));
    const float sInner(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetInnerPseudoLayer()), pCluster));
    const float cInner(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()), pClusterSeed));

    // Veto track-track end-to-end associations
    if (!pClusterSeed->IsAvailable() && !pCluster->IsAvailable() &&
        ((sOuter < m_nearbyTrackDistance) || (sInner < m_nearbyTrackDistance)))
    {
        return NONE;
    }

    // Association check 1(a), look for enclosed clusters
    if ((cOuter < m_nearbyClusterDistance && cInner < m_nearbyClusterDistance) &&
        (!checkSeedForward || (sInner > m_nearbyClusterDistance)) &&
        (!checkSeedBackward || (sOuter > m_nearbyClusterDistance)))
    {
        return STRONG;
    }

    // Association check 1(b), look for overlapping clusters
    if ((checkSeedForward == checkCandidateForward) && (checkSeedBackward == checkCandidateBackward))
    {
        if ((cInner < m_nearbyClusterDistance && sOuter < m_nearbyClusterDistance) &&
            (!checkSeedForward || (sInner > m_nearbyClusterDistance)) &&
            (!checkSeedBackward || (cOuter > m_nearbyClusterDistance)))
        {
            return STRONG;
        }

        if ((cOuter < m_nearbyClusterDistance && sInner < m_nearbyClusterDistance) &&
            (!checkSeedBackward || (sOuter > m_nearbyClusterDistance)) &&
            (!checkSeedForward || (cInner > m_nearbyClusterDistance)))
        {
            return STRONG;
        }
    }

    // Association check 2, look for branching clusters
    if ((!checkSeedForward || (sInner > m_remoteClusterDistance)) &&
        (!checkSeedBackward || (sOuter > m_remoteClusterDistance)) &&
        ((checkCandidateForward && (cInner < m_nearbyClusterDistance)) || (checkCandidateBackward && (cOuter < m_nearbyClusterDistance))))
    {
        return STANDARD;
    }

    // Association check 3, look any distance below threshold
    if ((sOuter < m_nearbyClusterDistance) || (cOuter < m_nearbyClusterDistance) || (sInner < m_nearbyClusterDistance) || (cInner < m_nearbyClusterDistance))
        return SINGLE_ORDER;

    return NONE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ClusterCharacterisationAlgorithm::GetFigureOfMerit(const SeedAssociationList &seedAssociationList) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *pVertex((1 == pVertexList->size()) ? *(pVertexList->begin()) : NULL);

    if (m_useMCFigureOfMerit)
    {
        return this->GetMCFigureOfMerit(seedAssociationList);
    }
    else if ((NULL != pVertex) && (VERTEX_3D == pVertex->GetVertexType()))
    {
        return this->GetRecoFigureOfMerit(pVertex, seedAssociationList);
    }
    else
    {
        // ATTN Consistently returning same value will accept all candidate cluster merges
        return -1.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ClusterCharacterisationAlgorithm::GetMCFigureOfMerit(const SeedAssociationList &seedAssociationList) const
{
    unsigned int nMatchedClusters(0), nClusters(0);

    for (SeedAssociationList::const_iterator iter1 = seedAssociationList.begin(), iter1End = seedAssociationList.end(); iter1 != iter1End; ++iter1)
    {
        Cluster *pParentCluster(iter1->first);
        ++nClusters;

        const MCParticle *pParentMCParticle(NULL);
        try {pParentMCParticle = MCParticleHelper::GetMainMCParticle(pParentCluster);} catch (StatusCodeException &) {}

        ++nMatchedClusters;
        const ClusterVector &associatedClusters(iter1->second);

        for (ClusterVector::const_iterator iter2 = associatedClusters.begin(), iter2End = associatedClusters.end(); iter2 != iter2End; ++iter2)
        {
            Cluster *pBranchCluster = *iter2;
            ++nClusters;

            const MCParticle *pDaughterMCParticle(NULL);
            try {pDaughterMCParticle = MCParticleHelper::GetMainMCParticle(pBranchCluster);} catch (StatusCodeException &) {}

            if (pParentMCParticle == pDaughterMCParticle)
            {
                ++nMatchedClusters;
            }
        }
    }

    if (0 == nClusters)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float figureOfMerit(static_cast<float>(nMatchedClusters) / static_cast<float>(nClusters));
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ClusterCharacterisationAlgorithm::GetRecoFigureOfMerit(const Vertex *const pVertex, const SeedAssociationList &seedAssociationList) const
{
    unsigned int nVertexAssociatedSeeds(0), nVertexAssociatedNonSeeds(0);

    for (SeedAssociationList::const_iterator iter = seedAssociationList.begin(), iterEnd = seedAssociationList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pSeedCluster(iter->first);
        const ClusterVector &associatedClusters(iter->second);

        const HitType hitType(LArClusterHelper::GetClusterHitType(pSeedCluster));
        const CartesianVector vertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

        LArPointingClusterList pointingClusterSeedList;
        try {pointingClusterSeedList.push_back(LArPointingCluster(pSeedCluster));} catch (StatusCodeException &) {}

        LArPointingClusterList pointingClusterNonSeedList;
        for (ClusterVector::const_iterator cIter = associatedClusters.begin(), cIterEnd = associatedClusters.end(); cIter != cIterEnd; ++cIter)
        {
            try {pointingClusterNonSeedList.push_back(LArPointingCluster(*cIter));} catch (StatusCodeException &) {}
        }

        nVertexAssociatedSeeds += this->GetNVertexConnections(vertex2D, pointingClusterSeedList);
        nVertexAssociatedNonSeeds += this->GetNVertexConnections(vertex2D, pointingClusterNonSeedList);
    }

    const float figureOfMerit(static_cast<float>(nVertexAssociatedSeeds) - static_cast<float>(nVertexAssociatedNonSeeds));
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------                                     

unsigned int ClusterCharacterisationAlgorithm::GetNVertexConnections(const CartesianVector &vertexPosition2D, const LArPointingClusterList &pointingClusterList) const
{
    unsigned int nConnections(0);

    for (LArPointingClusterList::const_iterator cIter = pointingClusterList.begin(), cIterEnd = pointingClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        if (LArPointingClusterHelper::IsNode(vertexPosition2D, cIter->GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
            LArPointingClusterHelper::IsNode(vertexPosition2D, cIter->GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
            LArPointingClusterHelper::IsEmission(vertexPosition2D, cIter->GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
            LArPointingClusterHelper::IsEmission(vertexPosition2D, cIter->GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance))
        {
            ++nConnections;
        }
    }

    return nConnections;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::SortClusters(const Cluster *const pLhs, const Cluster *const pRhs)
{
    CartesianVector innerCoordinateLhs(0.f, 0.f, 0.f), outerCoordinateLhs(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pLhs, innerCoordinateLhs, outerCoordinateLhs);
    const float dLhs2((outerCoordinateLhs - innerCoordinateLhs).GetMagnitudeSquared());

    CartesianVector innerCoordinateRhs(0.f, 0.f, 0.f), outerCoordinateRhs(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pRhs, innerCoordinateRhs, outerCoordinateRhs);
    const float dRhs2((outerCoordinateRhs - innerCoordinateRhs).GetMagnitudeSquared());

    return (dLhs2 > dRhs2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::GetInputPfoList(PfoList &pfoList) const
{
    for (StringVector::const_iterator iter = m_inputPfoListNames.begin(), iterEnd = m_inputPfoListNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList = NULL;

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *iter, pPfoList))
        {
            std::cout << "ClusterCharacterisationAlgorithm : pfo list " << *iter << " unavailable." << std::endl;
            continue;
        }

        pfoList.insert(pPfoList->begin(), pPfoList->end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::StoreNCaloHitsPerCluster(const Cluster *const pCluster, ClusterInfoMap &clusterInfoMap) const
{
    // ATTN Stores only first value provided per cluster
    if (clusterInfoMap.count(pCluster))
        return;

    if (!clusterInfoMap.insert(ClusterInfoMap::value_type(pCluster, pCluster->GetNCaloHits())).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::StoreNBranchesPerCluster(const Cluster *const pCluster, const ClusterVector &branchList,
    ClusterInfoMap &clusterInfoMap) const
{
    // ATTN Stores total number of branches linked to eventual parent cluster
    ClusterInfoMap::const_iterator iIter = clusterInfoMap.find(pCluster);
    unsigned int nBranchesSum((clusterInfoMap.end() == iIter) ? 0 : iIter->second);

    for (ClusterVector::const_iterator iter = branchList.begin(), iterEnd = branchList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pBranchCluster = *iter;
        ClusterInfoMap::iterator bIter = clusterInfoMap.find(pBranchCluster);

        const unsigned int nBranches((clusterInfoMap.end() == bIter) ? 0 : bIter->second);
        nBranchesSum += (1 + nBranches);

        if (clusterInfoMap.end() != bIter)
            clusterInfoMap.erase(bIter);
    }

    clusterInfoMap[pCluster] = nBranchesSum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::RemoveShowerPfos(const ClusterInfoMap &nCaloHitsPerCluster, const ClusterInfoMap &nBranchesPerCluster, PfoList &pfoList) const
{
    for (StringVector::const_iterator listIter = m_inputClusterListNames.begin(), listIterEnd = m_inputClusterListNames.end(); listIter != listIterEnd; ++listIter)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *listIter, pClusterList));

        for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        {
            try
            {
                Cluster *pCluster = *iter;

                if (pCluster->IsAvailable())
                    continue;

                ClusterInfoMap::const_iterator nCaloHitsIter = nCaloHitsPerCluster.find(pCluster);
                ClusterInfoMap::const_iterator nBranchesIter = nBranchesPerCluster.find(pCluster);

                if ((nCaloHitsPerCluster.end() == nCaloHitsIter) || (nBranchesPerCluster.end() == nBranchesIter))
                    continue;

                if (0 == nCaloHitsIter->second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const float nCaloHitsRatio(static_cast<float>(pCluster->GetNCaloHits()) / static_cast<float>(nCaloHitsIter->second));
                const unsigned int nBranches(nBranchesIter->second);

                if ((nBranches < m_showerLikeNBranches) && (nCaloHitsRatio < m_showerLikeCaloHitRatio))
                    continue;

                Pfo *pTargetPfo = NULL;
                this->FindTargetPfo(pCluster, pfoList, pTargetPfo);

                pfoList.erase(pTargetPfo);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pTargetPfo));
            }
            catch (StatusCodeException &)
            {
                std::cout << "ClusterCharacterisationAlgorithm: Unable to remove shower-like pfo." << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::FindTargetPfo(Cluster *const pCluster, const PfoList &pfoList, Pfo *&pTargetPfo) const
{
    pTargetPfo = NULL;

    for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        Pfo *pPfo = *iter;

        if (pPfo->GetClusterList().count(pCluster))
        {
            pTargetPfo = pPfo;
            return;
        }
    }

    if (NULL == pTargetPfo)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputPfoListNames", m_inputPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NearbyTrackDistance", m_nearbyTrackDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NearbyClusterDistance", m_nearbyClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RemoteClusterDistance", m_remoteClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseExistingPfosAsSeeds", m_useExistingPfosAsSeeds));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseExistingPfosAsBranches", m_useExistingPfosAsBranches));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DirectionTanAngle", m_directionTanAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DirectionApexShift", m_directionApexShift));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseMCFigureOfMerit", m_useMCFigureOfMerit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseFirstImprovedSeed", m_useFirstImprovedSeed));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldRemoveShowerPfos", m_shouldRemoveShowerPfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerLikeNBranches", m_showerLikeNBranches));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerLikeCaloHitRatio", m_showerLikeCaloHitRatio));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexLongitudinalDistance", m_maxVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexAngularAllowance", m_vertexAngularAllowance));

    return SeedGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
