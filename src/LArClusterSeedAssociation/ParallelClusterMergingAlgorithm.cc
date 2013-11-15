/**
 *  @file   LArContent/src/LArClusterSeedAssociation/ParallelClusterMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the bounded cluster merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSeedAssociation/ParallelClusterMergingAlgorithm.h"

#include "LArHelpers/LArParticleIdHelper.h"
#include "LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar
{

StatusCode ParallelClusterMergingAlgorithm::Run()
{
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    const ClusterList *pNonSeedClusterList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_INITIALIZED != statusCode))
        return statusCode;

    if (STATUS_CODE_NOT_INITIALIZED == statusCode)
        return STATUS_CODE_SUCCESS;

    bool foundAssociations(true);
    unsigned int nIterations(0);

    while (foundAssociations && (nIterations++ < m_maxNumIterations))
    {
        foundAssociations = false;
        ClusterMergeMap clusterMergeMap;

        // Calculate seed properties and cache results
        ClusterPropertiesList seedPropertiesList;

        for (ClusterList::const_iterator iter = pSeedClusterList->begin(), iterEnd = pSeedClusterList->end(); iter != iterEnd; ++iter)
        {
            seedPropertiesList.push_back(ClusterProperties(*iter));
        }

        // Associate each non-seed cluster with its preferred seed cluster
        for (ClusterList::const_iterator iterI = pNonSeedClusterList->begin(), iterEndI = pNonSeedClusterList->end(); iterI != iterEndI; ++iterI)
        {
            Cluster *pTargetCluster = *iterI;

            if (!this->IsCleanCluster(pTargetCluster))
                continue;

            Cluster *pBestSeedCluster(NULL);
            unsigned int nClusters(0);
            unsigned int nCloseClusters(0);
            float closestDistanceSquared(m_clusterWindowRadiusSquared);

            for (ClusterPropertiesList::const_iterator iterJ = seedPropertiesList.begin(), iterEndJ = seedPropertiesList.end(); iterJ != iterEndJ; ++iterJ)
            {
                const ClusterProperties &seedProperties(*iterJ);
                Cluster *pSeedCluster = seedProperties.GetCluster();

                if (!this->IsParallel(pSeedCluster, pTargetCluster))
                    continue;

                const float thisDistanceSquared(this->GetClosestDistanceSquared(seedProperties, pTargetCluster));
                const bool isGoodMerge(this->IsGoodMerge(seedProperties, pTargetCluster));

                bool foundCluster(false);

                if (thisDistanceSquared < m_clusterWindowRadiusSquared)
                {
                    if (!seedProperties.IsMuon())
                    {
                        ++nClusters;
                        foundCluster = true;
                    }

                    if (thisDistanceSquared < m_clusterWindowCloseRadiusSquared)
                    {
                        ++nCloseClusters;
                        foundCluster = true;
                    }
                }

                if (foundCluster && (thisDistanceSquared < closestDistanceSquared))
                {
                    closestDistanceSquared = thisDistanceSquared;

                    if (isGoodMerge)
                    {
                        pBestSeedCluster = pSeedCluster;
                    }
                    else
                    {
                        pBestSeedCluster = NULL;
                    }
                }
            }

            if ((NULL != pBestSeedCluster) && (nCloseClusters <= 1))
            {
                clusterMergeMap[pBestSeedCluster].insert(pTargetCluster);
                foundAssociations = true;
            }
        }

// for(ClusterMergeMap::const_iterator iter = clusterMergeMap.begin(), iterEnd = clusterMergeMap.end(); iter != iterEnd; ++iter) 
// {
// ClusterList pTempList, dTempList;
// pTempList.insert(iter->first);
// for( ClusterList::const_iterator dauIter = iter->second.begin(), dauIterEnd = iter->second.end(); dauIter != dauIterEnd; ++dauIter)
// {
// dTempList.insert(*dauIter);
// }
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&pTempList, "PCMParent", RED);
// PandoraMonitoringApi::VisualizeClusters(&dTempList, "PCMDaughter", BLUE);
// PandoraMonitoringApi::ViewEvent();
// }

        // Make Merges
        for (ClusterMergeMap::const_iterator iter = clusterMergeMap.begin(), iterEnd = clusterMergeMap.end(); iter != iterEnd; ++iter)
        {
            for (ClusterList::const_iterator dauIter = iter->second.begin(), dauIterEnd = iter->second.end(); dauIter != dauIterEnd; ++dauIter)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, iter->first, *dauIter,
                    m_seedClusterListName, m_nonSeedClusterListName));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParallelClusterMergingAlgorithm::IsCleanCluster(const Cluster* pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minClusterSize)
        return false;

    if ((1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer()) < m_minClusterSize)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParallelClusterMergingAlgorithm::IsParallel(const Cluster *const pSeedCluster, const Cluster *const pTargetCluster) const
{
    // Require an overlap between the seed and target cluster
    if ((pTargetCluster->GetOuterPseudoLayer() < pSeedCluster->GetInnerPseudoLayer()) ||
        (pSeedCluster->GetOuterPseudoLayer() < pTargetCluster->GetInnerPseudoLayer()))
    {
         return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParallelClusterMergingAlgorithm::IsGoodMerge(const ClusterProperties &seedProperties, const Cluster *const pTargetCluster) const
{
    const double closestDistanceToVertexSquared(this->GetClosestDistanceToVertexSquared(seedProperties, pTargetCluster));
    const double closestDistanceSquared(this->GetClosestDistanceSquared(seedProperties, pTargetCluster));

    if (closestDistanceToVertexSquared - closestDistanceSquared < m_vertexVetoRadiusSquared)
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ParallelClusterMergingAlgorithm::GetClosestDistanceToVertexSquared(const ClusterProperties &seedProperties, const Cluster *const pTargetCluster) const
{
    const Cluster *const pSeedCluster(seedProperties.GetCluster());
    const CartesianVector seedInnerCentroid(pSeedCluster->GetCentroid(pSeedCluster->GetInnerPseudoLayer()));
    const CartesianVector seedOuterCentroid(pSeedCluster->GetCentroid(pSeedCluster->GetOuterPseudoLayer()));
    const CartesianVector targetInnerCentroid(pTargetCluster->GetCentroid(pTargetCluster->GetInnerPseudoLayer()));
    const CartesianVector targetOuterCentroid(pTargetCluster->GetCentroid(pTargetCluster->GetOuterPseudoLayer()));

    float closestDistanceSquared( std::numeric_limits<float>::max() );

    // if seed is track-like, or could be going forwards, check the inner layer associations
    if (seedProperties.IsForwardInZ() || !seedProperties.IsBackwardInZ() || seedProperties.IsMuon())
    {
        const float innerDistanceSquared = std::min((seedInnerCentroid-targetInnerCentroid).GetMagnitudeSquared(),
            (seedInnerCentroid-targetOuterCentroid).GetMagnitudeSquared());

        if (innerDistanceSquared < closestDistanceSquared)
            closestDistanceSquared = innerDistanceSquared;
    }

    // if seed is track-like, or could be going backwards, check the outer layer associations
    if (!seedProperties.IsForwardInZ() || seedProperties.IsBackwardInZ() || seedProperties.IsMuon())
    {
        const float outerDistanceSquared = std::min((seedOuterCentroid-targetInnerCentroid).GetMagnitudeSquared(),
             (seedOuterCentroid-targetOuterCentroid).GetMagnitudeSquared());

        if (outerDistanceSquared < closestDistanceSquared)
            closestDistanceSquared = outerDistanceSquared;
    }

    return closestDistanceSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ParallelClusterMergingAlgorithm::GetClosestDistanceSquared(const ClusterProperties &seedProperties, const Cluster *const pTargetCluster) const
{
    const Cluster *const pSeedCluster(seedProperties.GetCluster()); 
    const CartesianVector seedInnerCentroid(pSeedCluster->GetCentroid(pSeedCluster->GetInnerPseudoLayer()));
    const CartesianVector seedOuterCentroid(pSeedCluster->GetCentroid(pSeedCluster->GetOuterPseudoLayer()));
    const CartesianVector targetInnerCentroid(pTargetCluster->GetCentroid(pTargetCluster->GetInnerPseudoLayer()));
    const CartesianVector targetOuterCentroid(pTargetCluster->GetCentroid(pTargetCluster->GetOuterPseudoLayer()));

    float closestDistanceSquared(std::numeric_limits<float>::max());

    bool checkTargetInnerLayer(false);
    bool checkTargetOuterLayer(false);

    // seed cluster could be going forwards
    if (seedProperties.IsForwardInZ() || !seedProperties.IsBackwardInZ())
    {
        if ((targetInnerCentroid - seedInnerCentroid).GetMagnitudeSquared() < (targetOuterCentroid - seedInnerCentroid).GetMagnitudeSquared())
        {
            checkTargetInnerLayer = true;
        }
        else
        {
            checkTargetOuterLayer = true;
        }
    }

    // seed cluster could be going backwards
    if (!seedProperties.IsForwardInZ() || seedProperties.IsBackwardInZ())
    {
        if ((targetOuterCentroid - seedOuterCentroid).GetMagnitudeSquared() < (targetInnerCentroid - seedOuterCentroid).GetMagnitudeSquared())
        {
            checkTargetOuterLayer = true;
        }
        else
        {
            checkTargetInnerLayer = true;
        }
    }

    // check the inner layer of the target cluster
    if (checkTargetInnerLayer)
    {
        const float innerDistanceSquared(this->GetClosestDistanceSquared(pSeedCluster, pTargetCluster->GetCentroid(pTargetCluster->GetInnerPseudoLayer())));

        if (innerDistanceSquared < closestDistanceSquared)
            closestDistanceSquared = innerDistanceSquared;
    }

    // check the outer layer of the target cluster
    if (checkTargetOuterLayer)
    {
        const float outerDistanceSquared(this->GetClosestDistanceSquared(pSeedCluster, pTargetCluster->GetCentroid(pTargetCluster->GetOuterPseudoLayer())));

        if (outerDistanceSquared<closestDistanceSquared)
            closestDistanceSquared = outerDistanceSquared;
    }

    return closestDistanceSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ParallelClusterMergingAlgorithm::GetClosestDistanceSquared(const Cluster *pCluster, const CartesianVector &position) const
{
    // Calculate a search window for the given cluster
    const unsigned int positionLayer(GeometryHelper::GetPseudoLayer(position));
    const unsigned int windowLayerMin((positionLayer >= m_clusterWindowLayers) ? positionLayer - m_clusterWindowLayers : 0);
    const unsigned int windowLayerMax(positionLayer + m_clusterWindowLayers);

    const unsigned int firstLayer(std::max(pCluster->GetInnerPseudoLayer(), windowLayerMin));
    const unsigned int lastLayer(std::min(pCluster->GetOuterPseudoLayer(), windowLayerMax));

    // Calculate closest distance from the hits in the search window
    float closestDistanceSquared(std::numeric_limits<float>::max());
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (unsigned int iLayer = firstLayer; iLayer <= lastLayer; ++iLayer)
    {
        OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.find(iLayer);

        if (orderedCaloHitList.end() == layerIter)
            continue;

        const CaloHitList *pCaloHitList = layerIter->second;

        for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const float distanceSquared = ((*hitIter)->GetPositionVector() - position).GetMagnitudeSquared();

            if (distanceSquared < closestDistanceSquared)
                closestDistanceSquared = distanceSquared;
        }
    }

    return closestDistanceSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ParallelClusterMergingAlgorithm::ClusterProperties::ClusterProperties(Cluster *const pCluster) :
    m_pCluster(pCluster),
    m_isForwardInZ((LArVertexHelper::DoesCurrentVertexExist()) ? LArVertexHelper::IsForwardInZ(pCluster) : false),
    m_isBackwardInZ((LArVertexHelper::DoesCurrentVertexExist()) ? LArVertexHelper::IsBackwardInZ(pCluster) : false),
    m_isMuon(LArParticleIdHelper::LArMuonId(pCluster))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParallelClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    m_vertexVetoRadius = 10.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexVetoRadius", m_vertexVetoRadius));
    m_vertexVetoRadiusSquared = m_vertexVetoRadius*m_vertexVetoRadius;

    m_clusterWindowRadius = 5.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterWindowRadius", m_clusterWindowRadius));
    m_clusterWindowRadiusSquared = m_clusterWindowRadius*m_clusterWindowRadius;

    m_clusterWindowCloseRadius = 2.5; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterWindowCloseRadius", m_clusterWindowCloseRadius));
    m_clusterWindowCloseRadiusSquared = m_clusterWindowCloseRadius*m_clusterWindowCloseRadius;

    m_clusterWindowLayers = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterWindowLayers", m_clusterWindowLayers));

    m_minClusterSize = 4;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterSize", m_minClusterSize));

    m_maxNumIterations = 3;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxIterations", m_maxNumIterations));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
