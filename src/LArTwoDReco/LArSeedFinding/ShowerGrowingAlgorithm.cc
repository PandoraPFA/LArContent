/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/ShowerGrowingAlgorithm.cc
 * 
 *  @brief  Implementation of the shower growing algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArTwoDReco/LArSeedFinding/ShowerGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerGrowingAlgorithm::ShowerGrowingAlgorithm() :
    m_minCaloHitsPerCluster(5),
    m_nearbyTrackDistance(1.f),
    m_nearbyClusterDistance(2.5f),
    m_remoteClusterDistance(10.f),
    m_directionTanAngle(1.732f),
    m_directionApexShift(0.333f),
    m_recursiveMode(false),
    m_useFirstImprovedSeed(false),
    m_useMCFigureOfMerit(false),
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

StatusCode ShowerGrowingAlgorithm::Run()
{
    PfoList pfoList;
    this->GetInputPfoList(pfoList);
    ClusterInfoMap nCaloHitsPerCluster, nBranchesPerCluster;

    for (StringVector::const_iterator listIter = m_inputClusterListNames.begin(), listIterEnd = m_inputClusterListNames.end(); listIter != listIterEnd; ++listIter)
    {
        try
        {
            const ClusterList *pClusterList = NULL;
            const StatusCode listStatusCode(PandoraContentApi::GetList(*this, *listIter, pClusterList));

            if (STATUS_CODE_NOT_INITIALIZED == listStatusCode)
            {
                std::cout << "ShowerGrowingAlgorithm: cluster list not found " << *listIter << std::endl;
                continue;
            }

            if (STATUS_CODE_SUCCESS != listStatusCode)
                return listStatusCode;

            m_clusterDirectionMap.clear();

            if (!m_recursiveMode)
            {
                this->SimpleModeShowerGrowing(pClusterList, *listIter, pfoList, nCaloHitsPerCluster, nBranchesPerCluster);
            }
            else
            {
                this->RecursiveModeShowerGrowing(pClusterList, *listIter, pfoList, nCaloHitsPerCluster, nBranchesPerCluster);
            }

            m_clusterDirectionMap.clear();
        }
        catch (StatusCodeException &statusCodeException)
        {
            m_clusterDirectionMap.clear();
            throw statusCodeException;
        }
    }

    if (m_shouldRemoveShowerPfos && !nCaloHitsPerCluster.empty() && !nBranchesPerCluster.empty())
        this->RemoveShowerPfos(nCaloHitsPerCluster, nBranchesPerCluster, pfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::SimpleModeShowerGrowing(const ClusterList *const pClusterList, const std::string &clusterListName,
    PfoList &pfoList, ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex(((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : NULL);

    ClusterList usedClusters;

    // Pick up all showers starting at vertex
    if (NULL != pVertex)
    {
        ClusterVector seedClusters;
        this->GetAllVertexSeedCandidates(pClusterList, pVertex, seedClusters);

        SeedAssociationList vertexSeedAssociationList;
        this->GetSeedAssociationList(seedClusters, pClusterList, vertexSeedAssociationList);
        this->ProcessSeedAssociationDetails(vertexSeedAssociationList, clusterListName, pfoList, usedClusters, nCaloHitsPerCluster, nBranchesPerCluster);
    }

    // Non-vertex showers
    const Cluster *pSeedCluster(NULL);

    while (this->GetNextSeedCandidate(pClusterList, usedClusters, pSeedCluster))
    {
        SeedAssociationList seedAssociationList;
        this->GetSeedAssociationList(ClusterVector(1, pSeedCluster), pClusterList, seedAssociationList);
        this->ProcessSeedAssociationDetails(seedAssociationList, clusterListName, pfoList, usedClusters, nCaloHitsPerCluster, nBranchesPerCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::RecursiveModeShowerGrowing(const ClusterList *const pClusterList, const std::string &clusterListName,
    PfoList &pfoList, ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const
{
    ClusterList usedClusters;
    const Cluster *pSeedCluster(NULL);

    while (this->GetNextSeedCandidate(pClusterList, usedClusters, pSeedCluster))
    {
        SeedAssociationList seedAssociationList;
        this->GetSeedAssociationList(ClusterVector(1, pSeedCluster), pClusterList, seedAssociationList);

        SeedAssociationList finalSeedAssociationList;
        this->CheckSeedAssociationList(seedAssociationList.begin(), finalSeedAssociationList);
        this->ProcessSeedAssociationDetails(finalSeedAssociationList, clusterListName, pfoList, usedClusters, nCaloHitsPerCluster, nBranchesPerCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerGrowingAlgorithm::GetNextSeedCandidate(const ClusterList *const pClusterList, const ClusterList &usedClusters,
    const Cluster *&pSeedCluster) const
{
    pSeedCluster = NULL;

    ClusterVector clusterVector;
    clusterVector.insert(clusterVector.end(), pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), ShowerGrowingAlgorithm::SortClusters);

    for (const Cluster *const pCluster : clusterVector)
    {
        if (usedClusters.count(pCluster))
            continue;

        if (MU_MINUS == std::abs(pCluster->GetParticleIdFlag()))
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        pSeedCluster = pCluster;
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::GetAllVertexSeedCandidates(const ClusterList *const pClusterList, const Vertex *const pVertex,
    ClusterVector &seedClusters) const
{
    ClusterVector clusterVector;
    clusterVector.insert(clusterVector.end(), pClusterList->begin(), pClusterList->end());

    if (clusterVector.empty())
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(clusterVector.at(0)));
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

    for (const Cluster *const pCluster : clusterVector)
    {
        if (MU_MINUS == std::abs(pCluster->GetParticleIdFlag()))
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        try
        {
            if (this->IsVertexAssociated(LArPointingCluster(pCluster), vertexPosition2D))
                seedClusters.push_back(pCluster);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::GetSeedAssociationList(const ClusterVector &particleSeedVector, const ClusterList *const pClusterList,
    SeedAssociationList &seedAssociationList) const
{
    if (particleSeedVector.empty())
        return;

    ClusterVector candidateClusters;
    ClusterList seedClusters(particleSeedVector.begin(), particleSeedVector.end());

    const ClusterList clusterList(*pClusterList);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCandidateCluster = *iter;

        if (seedClusters.count(pCandidateCluster))
            continue;

        if (MU_MINUS == std::abs(pCandidateCluster->GetParticleIdFlag()))
            continue;

        if (pCandidateCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        candidateClusters.push_back(pCandidateCluster);
    }

    std::sort(candidateClusters.begin(), candidateClusters.end(), ShowerGrowingAlgorithm::SortClusters);
    ClusterUsageMap forwardUsageMap, backwardUsageMap;

    for (ClusterVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        this->FindAssociatedClusters(*iter, candidateClusters, forwardUsageMap, backwardUsageMap);
    }

    this->IdentifyClusterMerges(particleSeedVector, backwardUsageMap, seedAssociationList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::CheckSeedAssociationList(SeedAssociationList::const_iterator seedIter, SeedAssociationList &finalSeedAssociationList) const
{
    ClusterList usedClusters, availableClusters;
    usedClusters.insert(seedIter->first);
    availableClusters.insert(seedIter->second.begin(), seedIter->second.end());

    SeedAssociationList seedAssociationList;
    seedAssociationList.insert(SeedAssociationList::value_type(seedIter->first, seedIter->second));
    const float originalFigureOfMerit(this->GetFigureOfMerit(seedAssociationList));
    float bestFigureOfMerit(originalFigureOfMerit);

    const Cluster *pTrialSeedCluster = NULL;
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

void ShowerGrowingAlgorithm::ProcessBranchClusters(const Cluster *const pParentCluster, const ClusterVector &branchClusters, const std::string &listName,
    PfoList &pfoList) const
{
    m_clusterDirectionMap.erase(pParentCluster);

    for (ClusterVector::const_iterator iter = branchClusters.begin(), iterEnd = branchClusters.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pBranchCluster(*iter);

        if (!pBranchCluster->IsAvailable() && m_shouldRemoveShowerPfos)
        {
            const Pfo *pTargetPfo(NULL);
            this->FindTargetPfo(pBranchCluster, pfoList, pTargetPfo);
            pfoList.erase(pTargetPfo);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pTargetPfo));
        }

        if (pBranchCluster->IsAvailable())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pBranchCluster, listName, listName));
        }

        m_clusterDirectionMap.erase(pBranchCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerGrowingAlgorithm::AssociationType ShowerGrowingAlgorithm::AreClustersAssociated(const Cluster *const pClusterSeed, const Cluster *const pCluster) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex(((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : NULL);

    // Direction of seed cluster (cache for efficiency)
    ClusterDirectionMap::const_iterator seedIter = m_clusterDirectionMap.find(pClusterSeed);

    if (m_clusterDirectionMap.end() == seedIter)
    {
        const LArVertexHelper::ClusterDirection direction((NULL == pVertex) ? LArVertexHelper::DIRECTION_UNKNOWN :
            LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex, pClusterSeed, m_directionTanAngle, m_directionApexShift));
        seedIter = m_clusterDirectionMap.insert(ClusterDirectionMap::value_type(pClusterSeed, direction)).first;
    }

    const LArVertexHelper::ClusterDirection seedDirection(seedIter->second);
    const bool checkSeedForward(seedDirection != LArVertexHelper::DIRECTION_BACKWARD_IN_Z);
    const bool checkSeedBackward(seedDirection != LArVertexHelper::DIRECTION_FORWARD_IN_Z);

    // Direction of candidate cluster (cache for efficiency)
    ClusterDirectionMap::const_iterator candIter = m_clusterDirectionMap.find(pCluster);

    if (m_clusterDirectionMap.end() == candIter)
    {
        const LArVertexHelper::ClusterDirection direction((NULL == pVertex) ? LArVertexHelper::DIRECTION_UNKNOWN :
            LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex, pCluster, m_directionTanAngle, m_directionApexShift));
        candIter = m_clusterDirectionMap.insert(ClusterDirectionMap::value_type(pCluster, direction)).first;
    }

    const LArVertexHelper::ClusterDirection candidateDirection(candIter->second);
    const bool checkCandidateForward(candidateDirection != LArVertexHelper::DIRECTION_BACKWARD_IN_Z);
    const bool checkCandidateBackward(candidateDirection != LArVertexHelper::DIRECTION_FORWARD_IN_Z);

    // Calculate distances of association
    const float sOuter(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetOuterPseudoLayer()), pCluster));
    const float cOuter(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()), pClusterSeed));
    const float sInner(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetInnerPseudoLayer()), pCluster));
    const float cInner(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()), pClusterSeed));

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

float ShowerGrowingAlgorithm::GetFigureOfMerit(const SeedAssociationList &seedAssociationList) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex(((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : NULL);

    if (m_useMCFigureOfMerit)
    {
        return this->GetMCFigureOfMerit(seedAssociationList);
    }
    else if (NULL != pVertex)
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

float ShowerGrowingAlgorithm::GetMCFigureOfMerit(const SeedAssociationList &seedAssociationList) const
{
    unsigned int nMatchedClusters(0), nClusters(0);

    for (SeedAssociationList::const_iterator iter1 = seedAssociationList.begin(), iter1End = seedAssociationList.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pParentCluster(iter1->first);
        ++nClusters;

        const MCParticle *pParentMCParticle(NULL);
        try {pParentMCParticle = MCParticleHelper::GetMainMCParticle(pParentCluster);} catch (StatusCodeException &) {}

        ++nMatchedClusters;
        const ClusterVector &associatedClusters(iter1->second);

        for (ClusterVector::const_iterator iter2 = associatedClusters.begin(), iter2End = associatedClusters.end(); iter2 != iter2End; ++iter2)
        {
            const Cluster *const pBranchCluster = *iter2;
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

float ShowerGrowingAlgorithm::GetRecoFigureOfMerit(const Vertex *const pVertex, const SeedAssociationList &seedAssociationList) const
{
    unsigned int nVertexAssociatedSeeds(0), nVertexAssociatedNonSeeds(0);

    for (SeedAssociationList::const_iterator iter = seedAssociationList.begin(), iterEnd = seedAssociationList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pSeedCluster(iter->first);
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

unsigned int ShowerGrowingAlgorithm::GetNVertexConnections(const CartesianVector &vertexPosition2D, const LArPointingClusterList &pointingClusterList) const
{
    unsigned int nConnections(0);

    for (LArPointingClusterList::const_iterator cIter = pointingClusterList.begin(), cIterEnd = pointingClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        if (this->IsVertexAssociated(*cIter, vertexPosition2D))
            ++nConnections;
    }

    return nConnections;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerGrowingAlgorithm::IsVertexAssociated(const LArPointingCluster &pointingCluster, const CartesianVector &vertexPosition2D) const
{
    return (LArPointingClusterHelper::IsNode(vertexPosition2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(vertexPosition2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsEmission(vertexPosition2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(vertexPosition2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerGrowingAlgorithm::SortClusters(const Cluster *const pLhs, const Cluster *const pRhs)
{
    CartesianVector innerCoordinateLhs(0.f, 0.f, 0.f), outerCoordinateLhs(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pLhs, innerCoordinateLhs, outerCoordinateLhs);
    const float dLhs2((outerCoordinateLhs - innerCoordinateLhs).GetMagnitudeSquared());

    CartesianVector innerCoordinateRhs(0.f, 0.f, 0.f), outerCoordinateRhs(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pRhs, innerCoordinateRhs, outerCoordinateRhs);
    const float dRhs2((outerCoordinateRhs - innerCoordinateRhs).GetMagnitudeSquared());

    return (dLhs2 > dRhs2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::GetInputPfoList(PfoList &pfoList) const
{
    for (StringVector::const_iterator iter = m_inputPfoListNames.begin(), iterEnd = m_inputPfoListNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList = NULL;

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *iter, pPfoList))
        { 
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ShowerGrowingAlgorithm : pfo list " << *iter << " unavailable." << std::endl;
            continue;
        }

        pfoList.insert(pPfoList->begin(), pPfoList->end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::ProcessSeedAssociationDetails(const SeedAssociationList &seedAssociationList, const std::string &clusterListName,
    PfoList &pfoList, ClusterList &usedClusters, ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const
{
    for (SeedAssociationList::const_iterator iter = seedAssociationList.begin(), iterEnd = seedAssociationList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pParentCluster(iter->first);
        const ClusterVector &branchClusters(iter->second);

        this->StoreNCaloHitsPerCluster(pParentCluster, nCaloHitsPerCluster);
        this->StoreNBranchesPerCluster(pParentCluster, branchClusters, nBranchesPerCluster);
        this->ProcessBranchClusters(pParentCluster, branchClusters, clusterListName, pfoList);

        usedClusters.insert(pParentCluster);
        usedClusters.insert(branchClusters.begin(), branchClusters.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::StoreNCaloHitsPerCluster(const Cluster *const pCluster, ClusterInfoMap &clusterInfoMap) const
{
    // ATTN Stores only first value provided per cluster
    if (clusterInfoMap.count(pCluster))
        return;

    if (!clusterInfoMap.insert(ClusterInfoMap::value_type(pCluster, pCluster->GetNCaloHits())).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::StoreNBranchesPerCluster(const Cluster *const pCluster, const ClusterVector &branchList,
    ClusterInfoMap &clusterInfoMap) const
{
    // ATTN Stores total number of branches linked to eventual parent cluster
    ClusterInfoMap::const_iterator iIter = clusterInfoMap.find(pCluster);
    unsigned int nBranchesSum((clusterInfoMap.end() == iIter) ? 0 : iIter->second);

    for (ClusterVector::const_iterator iter = branchList.begin(), iterEnd = branchList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pBranchCluster = *iter;
        ClusterInfoMap::iterator bIter = clusterInfoMap.find(pBranchCluster);

        const unsigned int nBranches((clusterInfoMap.end() == bIter) ? 0 : bIter->second);
        nBranchesSum += (1 + nBranches);

        if (clusterInfoMap.end() != bIter)
            clusterInfoMap.erase(bIter);
    }

    clusterInfoMap[pCluster] = nBranchesSum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::RemoveShowerPfos(const ClusterInfoMap &nCaloHitsPerCluster, const ClusterInfoMap &nBranchesPerCluster, PfoList &pfoList) const
{
    for (StringVector::const_iterator listIter = m_inputClusterListNames.begin(), listIterEnd = m_inputClusterListNames.end(); listIter != listIterEnd; ++listIter)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *listIter, pClusterList));

        for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        {
            try
            {
                const Cluster *const pCluster = *iter;

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

                const Pfo *pTargetPfo = NULL;
                this->FindTargetPfo(pCluster, pfoList, pTargetPfo);

                pfoList.erase(pTargetPfo);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pTargetPfo));
            }
            catch (StatusCodeException &)
            {
                std::cout << "ShowerGrowingAlgorithm: Unable to remove shower-like pfo." << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::FindTargetPfo(const Cluster *const pCluster, const PfoList &pfoList, const Pfo *&pTargetPfo) const
{
    pTargetPfo = NULL;

    for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        const Pfo *const pPfo = *iter;

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

StatusCode ShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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
        "DirectionTanAngle", m_directionTanAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DirectionApexShift", m_directionApexShift));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RecursiveMode", m_recursiveMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseFirstImprovedSeed", m_useFirstImprovedSeed));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseMCFigureOfMerit", m_useMCFigureOfMerit));

    if (m_useMCFigureOfMerit && !m_recursiveMode)
    {
        std::cout << "ShowerGrowingAlgorithm: UseMCFigureOfMerit available only in recursive mode - set RecursiveMode to true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

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
