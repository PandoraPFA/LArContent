/**
 *  @file   KinkSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the kink splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "KinkSplittingAlgorithm.h"
#include "LArParticleId.h"

using namespace pandora;

namespace lar
{

StatusCode KinkSplittingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), Cluster::SortByInnerLayer);

    for (ClusterVector::iterator iter = clusterVector.begin(); iter != clusterVector.end(); ++iter)
    {
        Cluster* pCluster = *iter;

        if (!this->IsPossibleKink(pCluster))
            continue;

        LArParticleId::TwoDSlidingXZFitResult twoDSlidingXZFitResult;
        LArParticleId::LArTwoDSlidingXZFit(pCluster, twoDSlidingXZFitResult);

        unsigned int splitLayer(std::numeric_limits<unsigned int>::max());

        if (STATUS_CODE_SUCCESS != twoDSlidingXZFitResult.FindLargestScatter(splitLayer))
            continue;

        if ((splitLayer <= pCluster->GetInnerPseudoLayer()) || (splitLayer >= pCluster->GetOuterPseudoLayer()))
            continue;

//ClusterList tempList;
//tempList.insert(pCluster);
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//PandoraMonitoringApi::VisualizeClusters(&tempList, "SplitCluster", BLUE);
//PandoraMonitoringApi::ViewEvent();
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SplitCluster(pCluster, splitLayer));
        *iter = NULL;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool KinkSplittingAlgorithm::IsPossibleKink(const Cluster *const pCluster) const
{
    if ((1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer()) < m_minClusterLayers)
        return false;

    ClusterHelper::ClusterFitResult innerLayerFit, outerLayerFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pCluster, m_minClusterLayers, innerLayerFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pCluster, m_minClusterLayers, outerLayerFit));

    if (innerLayerFit.GetDirection().GetCosOpeningAngle(outerLayerFit.GetDirection()) > m_maxCosScatteringAngle)
        return false;

    if ((innerLayerFit.GetRms() > m_minScatteringRms) && (outerLayerFit.GetRms() > m_minScatteringRms))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkSplittingAlgorithm::SplitCluster(Cluster *const pCluster, const unsigned int splitLayer) const
{
    // Begin cluster fragmentation operations
    ClusterList clusterList;
    clusterList.insert(pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Create new clusters
    Cluster *pCluster1(NULL);
    Cluster *pCluster2(NULL);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        const unsigned int thisLayer(iter->first);

        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;
            Cluster *&pClusterToModify((thisLayer < splitLayer) ? pCluster1 : pCluster2);

            if (NULL == pClusterToModify)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, pCaloHit, pClusterToModify));
            }
            else
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddCaloHitToCluster(*this, pClusterToModify, pCaloHit));
            }
        }
    }

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    m_minScatteringRms = 0.15;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinScatteringRms", m_minScatteringRms));

    m_maxCosScatteringAngle = std::cos(M_PI * 15.f / 180.f);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxCosScatteringAngle", m_maxCosScatteringAngle));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
