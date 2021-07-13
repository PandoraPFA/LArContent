/**
 *  @file   larpandoracontent/LarTrackShowerId/MetricBasedShowerGrowingAlgorithm.cc
 * 
 *  @brief  Implementation of the metric based shower growing algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTrackShowerId/MetricBasedShowerGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MetricBasedShowerGrowingAlgorithm::MetricBasedShowerGrowingAlgorithm() :
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

StatusCode MetricBasedShowerGrowingAlgorithm::Run()
{
    PfoList pfoList;
    this->GetInputPfoList(pfoList);
    ClusterInfoMap nCaloHitsPerCluster, nBranchesPerCluster;

    for (StringVector::const_iterator listIter = m_inputClusterListNames.begin(), listIterEnd = m_inputClusterListNames.end(); listIter != listIterEnd; ++listIter)
    {
        try
        {
            const ClusterList *pClusterList = NULL;
            PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, *listIter, pClusterList));

            if (!pClusterList || pClusterList->empty())
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "MetricBasedShowerGrowingAlgorithm: unable to find cluster list " << *listIter << std::endl;

                continue;
            }

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

void MetricBasedShowerGrowingAlgorithm::SimpleModeShowerGrowing(const ClusterList *const pClusterList, const std::string &clusterListName,
    PfoList &pfoList, ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex(((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : NULL);

    ClusterSet usedClusters;

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

void MetricBasedShowerGrowingAlgorithm::RecursiveModeShowerGrowing(const ClusterList *const pClusterList, const std::string &clusterListName,
    PfoList &pfoList, ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const
{
    ClusterSet usedClusters;
    const Cluster *pSeedCluster(NULL);
//std::cout << "RecursiveModeShowerGrowing, starting " << std::endl;
    while (this->GetNextSeedCandidate(pClusterList, usedClusters, pSeedCluster))
    {
//std::cout << "RecursiveModeShowerGrowing, new seed, initial association list " << std::endl;
//PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1., -1., 1.f);

        SeedAssociationList seedAssociationList;
        this->GetSeedAssociationList(ClusterVector(1, pSeedCluster), pClusterList, seedAssociationList);

//int counter(0);
//for (const auto &mapEntry : seedAssociationList)
//{
//ClusterList seeds, branches;
//seeds.push_back(mapEntry.first);
//branches.insert(branches.end(), mapEntry.second.begin(), mapEntry.second.end());
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &seeds, std::string("seed_" + TypeToString(counter)).c_str(), RED);
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &branches, std::string("branches_" + TypeToString(counter)).c_str(), GREEN);
//++counter;
//}
//PandoraMonitoringApi::ViewEvent(this->GetPandora());

        SeedAssociationList finalSeedAssociationList;
        this->CheckSeedAssociationList(seedAssociationList.begin(), finalSeedAssociationList);
//std::cout << "RecursiveModeShowerGrowing, final seed association list " << std::endl;
//counter = 0;
//for (const auto &mapEntry : finalSeedAssociationList)
//{
//ClusterList seeds, branches;
//seeds.push_back(mapEntry.first);
//branches.insert(branches.end(), mapEntry.second.begin(), mapEntry.second.end());
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &seeds, std::string("seed_" + TypeToString(counter)).c_str(), RED);
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &branches, std::string("branches_" + TypeToString(counter)).c_str(), GREEN);
//++counter;
//}
//PandoraMonitoringApi::ViewEvent(this->GetPandora());

        this->ProcessSeedAssociationDetails(finalSeedAssociationList, clusterListName, pfoList, usedClusters, nCaloHitsPerCluster, nBranchesPerCluster);
    }
//std::cout << "RecursiveModeShowerGrowing, done " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MetricBasedShowerGrowingAlgorithm::GetNextSeedCandidate(const ClusterList *const pClusterList, const ClusterSet &usedClusters,
    const Cluster *&pSeedCluster) const
{
    pSeedCluster = NULL;

    ClusterVector clusterVector;
    clusterVector.insert(clusterVector.end(), pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), MetricBasedShowerGrowingAlgorithm::SortClusters);

    for (const Cluster *const pCluster : clusterVector)
    {
        if (usedClusters.count(pCluster))
            continue;

        if (MU_MINUS == std::abs(pCluster->GetParticleId()))
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        pSeedCluster = pCluster;
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MetricBasedShowerGrowingAlgorithm::GetAllVertexSeedCandidates(const ClusterList *const pClusterList, const Vertex *const pVertex,
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
        if (MU_MINUS == std::abs(pCluster->GetParticleId()))
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

    std::sort(seedClusters.begin(), seedClusters.end(), MetricBasedShowerGrowingAlgorithm::SortClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MetricBasedShowerGrowingAlgorithm::GetSeedAssociationList(const ClusterVector &particleSeedVector, const ClusterList *const pClusterList,
    SeedAssociationList &seedAssociationList) const
{
    if (particleSeedVector.empty())
        return;

    ClusterVector candidateClusters;
    const ClusterList clusterList(*pClusterList);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCandidateCluster = *iter;

        if (particleSeedVector.end() != std::find(particleSeedVector.begin(), particleSeedVector.end(), pCandidateCluster))
            continue;

        if (MU_MINUS == std::abs(pCandidateCluster->GetParticleId()))
            continue;

        if (pCandidateCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        candidateClusters.push_back(pCandidateCluster);
    }

    std::sort(candidateClusters.begin(), candidateClusters.end(), MetricBasedShowerGrowingAlgorithm::SortClusters);
    ClusterUsageMap forwardUsageMap, backwardUsageMap;

    for (ClusterVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        this->FindAssociatedClusters(*iter, candidateClusters, forwardUsageMap, backwardUsageMap);
    }

    this->IdentifyClusterMerges(particleSeedVector, backwardUsageMap, seedAssociationList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MetricBasedShowerGrowingAlgorithm::CheckSeedAssociationList(SeedAssociationList::const_iterator seedIter, SeedAssociationList &finalSeedAssociationList) const
{
    ClusterSet usedClusters;
    usedClusters.insert(seedIter->first);
    
    ClusterList availableClusters(seedIter->second.begin(), seedIter->second.end());

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
        ClusterList clusterList;
        for (const auto &mapEntry : bestSeedAssociationList) clusterList.push_back(mapEntry.first);
        clusterList.sort(LArClusterHelper::SortByNHits);

        for (const Cluster *const pCluster : clusterList)
            this->CheckSeedAssociationList(bestSeedAssociationList.find(pCluster), finalSeedAssociationList);
    }
    else
    {
        finalSeedAssociationList.insert(SeedAssociationList::value_type(seedIter->first, seedIter->second));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MetricBasedShowerGrowingAlgorithm::ProcessBranchClusters(const Cluster *const pParentCluster, const ClusterVector &branchClusters, const std::string &listName,
    PfoList &pfoList) const
{
    m_clusterDirectionMap.erase(pParentCluster);

    for (ClusterVector::const_iterator iter = branchClusters.begin(), iterEnd = branchClusters.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pBranchCluster(*iter);

        if (!pBranchCluster->IsAvailable() && m_shouldRemoveShowerPfos)
        {
            PfoList::iterator targetIter(pfoList.end());
            this->FindTargetPfo(pBranchCluster, pfoList, targetIter);

            const Pfo *const pTargetPfo(*targetIter);
            pfoList.erase(targetIter);

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

MetricBasedShowerGrowingAlgorithm::AssociationType MetricBasedShowerGrowingAlgorithm::AreClustersAssociated(const Cluster *const pClusterSeed, const Cluster *const pCluster) const
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

float MetricBasedShowerGrowingAlgorithm::GetFigureOfMerit(const SeedAssociationList &seedAssociationList) const
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

float MetricBasedShowerGrowingAlgorithm::GetMCFigureOfMerit(const SeedAssociationList &seedAssociationList) const
{
    unsigned int nMatchedClusters(0), nClusters(0);

    ClusterList clusterList;
    for (const auto &mapEntry : seedAssociationList) clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : clusterList)
    {
        ++nClusters;

        const MCParticle *pParentMCParticle(NULL);
        try 
        {
            pParentMCParticle = MCParticleHelper::GetMainMCParticle(pParentCluster);
            const unsigned int pdgCode(std::abs(pParentMCParticle->GetParticleId()));

            if (pdgCode == E_MINUS || pdgCode == PHOTON)
                pParentMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pParentMCParticle);
            //std::cout << "Parent pdg code " << pParentMCParticle->GetParticleId() << ", Energy " << pParentMCParticle->GetEnergy() << std::endl;
        }
        catch (const StatusCodeException &)
        {
        }

        ++nMatchedClusters;
        const ClusterVector &associatedClusters(seedAssociationList.at(pParentCluster));

        for (const Cluster *const pBranchCluster : associatedClusters)
        {
            ++nClusters;

            const MCParticle *pDaughterMCParticle(NULL);
            try
            {
                pDaughterMCParticle = MCParticleHelper::GetMainMCParticle(pBranchCluster);
                const unsigned int pdgCode(std::abs(pDaughterMCParticle->GetParticleId()));

                if (pdgCode == E_MINUS || pdgCode == PHOTON)
                    pDaughterMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pDaughterMCParticle);
                //std::cout << "Daughter pdg code " << pDaughterMCParticle->GetParticleId() << ", Energy " << pDaughterMCParticle->GetEnergy() << std::endl;    
            }
            catch (const StatusCodeException &)
            {
            }

            if (pParentMCParticle == pDaughterMCParticle)
                ++nMatchedClusters;
        }
    }

    if (0 == nClusters)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float figureOfMerit(static_cast<float>(nMatchedClusters) / static_cast<float>(nClusters));
//std::cout << "figureOfMerit " << figureOfMerit << std::endl;
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MetricBasedShowerGrowingAlgorithm::GetRecoFigureOfMerit(const Vertex *const pVertex, const SeedAssociationList &seedAssociationList) const
{
    unsigned int nVertexAssociatedSeeds(0), nVertexAssociatedNonSeeds(0);

    ClusterList clusterList;
    for (const auto &mapEntry : seedAssociationList) clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pSeedCluster : clusterList)
    {
        const ClusterVector &associatedClusters(seedAssociationList.at(pSeedCluster));

        const HitType hitType(LArClusterHelper::GetClusterHitType(pSeedCluster));
        const CartesianVector vertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

        LArPointingClusterList pointingClusterSeedList;
        try {pointingClusterSeedList.push_back(LArPointingCluster(pSeedCluster));} catch (StatusCodeException &) {}

        LArPointingClusterList pointingClusterNonSeedList;
        for (const Cluster *const pAssociatedCluster : associatedClusters)
        {
            try {pointingClusterNonSeedList.push_back(LArPointingCluster(pAssociatedCluster));} catch (StatusCodeException &) {}
        }

        nVertexAssociatedSeeds += this->GetNVertexConnections(vertex2D, pointingClusterSeedList);
        nVertexAssociatedNonSeeds += this->GetNVertexConnections(vertex2D, pointingClusterNonSeedList);
    }

    const float figureOfMerit(static_cast<float>(nVertexAssociatedSeeds) - static_cast<float>(nVertexAssociatedNonSeeds));
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------                                     

unsigned int MetricBasedShowerGrowingAlgorithm::GetNVertexConnections(const CartesianVector &vertexPosition2D, const LArPointingClusterList &pointingClusterList) const
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

bool MetricBasedShowerGrowingAlgorithm::IsVertexAssociated(const LArPointingCluster &pointingCluster, const CartesianVector &vertexPosition2D) const
{
    return (LArPointingClusterHelper::IsNode(vertexPosition2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(vertexPosition2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsEmission(vertexPosition2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(vertexPosition2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MetricBasedShowerGrowingAlgorithm::SortClusters(const Cluster *const pLhs, const Cluster *const pRhs)
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

void MetricBasedShowerGrowingAlgorithm::GetInputPfoList(PfoList &pfoList) const
{
    for (StringVector::const_iterator iter = m_inputPfoListNames.begin(), iterEnd = m_inputPfoListNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList = NULL;

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *iter, pPfoList))
        { 
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "MetricBasedShowerGrowingAlgorithm : pfo list " << *iter << " unavailable." << std::endl;
            continue;
        }

        pfoList.insert(pfoList.end(), pPfoList->begin(), pPfoList->end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MetricBasedShowerGrowingAlgorithm::ProcessSeedAssociationDetails(const SeedAssociationList &seedAssociationList, const std::string &clusterListName,
    PfoList &pfoList, ClusterSet &usedClusters, ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const
{
    ClusterList clusterList;
    for (const auto &mapEntry : seedAssociationList) clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : clusterList)
    {
        const ClusterVector &branchClusters(seedAssociationList.at(pParentCluster));

        this->StoreNCaloHitsPerCluster(pParentCluster, nCaloHitsPerCluster);
        this->StoreNBranchesPerCluster(pParentCluster, branchClusters, nBranchesPerCluster);
        this->ProcessBranchClusters(pParentCluster, branchClusters, clusterListName, pfoList);

        usedClusters.insert(pParentCluster);
        usedClusters.insert(branchClusters.begin(), branchClusters.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MetricBasedShowerGrowingAlgorithm::StoreNCaloHitsPerCluster(const Cluster *const pCluster, ClusterInfoMap &clusterInfoMap) const
{
    // ATTN Stores only first value provided per cluster
    if (clusterInfoMap.count(pCluster))
        return;

    if (!clusterInfoMap.insert(ClusterInfoMap::value_type(pCluster, pCluster->GetNCaloHits())).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MetricBasedShowerGrowingAlgorithm::StoreNBranchesPerCluster(const Cluster *const pCluster, const ClusterVector &branchList,
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

void MetricBasedShowerGrowingAlgorithm::RemoveShowerPfos(const ClusterInfoMap &nCaloHitsPerCluster, const ClusterInfoMap &nBranchesPerCluster, PfoList &pfoList) const
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

                PfoList::iterator targetIter(pfoList.end());
                this->FindTargetPfo(pCluster, pfoList, targetIter);

                const Pfo *const pTargetPfo(*targetIter);
                pfoList.erase(targetIter);

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pTargetPfo));
            }
            catch (StatusCodeException &)
            {
                std::cout << "MetricBasedShowerGrowingAlgorithm: Unable to remove shower-like pfo." << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MetricBasedShowerGrowingAlgorithm::FindTargetPfo(const Cluster *const pCluster, PfoList &pfoList, PfoList::iterator &targetIter) const
{
    for (PfoList::iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        const ClusterList &clusterList((*iter)->GetClusterList());

        if (clusterList.end() != std::find(clusterList.begin(), clusterList.end(), pCluster))
        {
            targetIter = iter;
            return;
        }
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MetricBasedShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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
        std::cout << "MetricBasedShowerGrowingAlgorithm: UseMCFigureOfMerit available only in recursive mode - set RecursiveMode to true" << std::endl;
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

    return BranchGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
