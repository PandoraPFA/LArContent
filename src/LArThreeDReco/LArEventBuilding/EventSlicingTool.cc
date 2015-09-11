/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/EventSlicingTool.cc
 * 
 *  @brief  Implementation of the event slicing tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/EventSlicingTool.h"

using namespace pandora;

namespace lar_content
{

void EventSlicingTool::Slice(NeutrinoParentAlgorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
    const HitTypeToNameMap &clusterListNames, SliceList &sliceList)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    ClusterToPfoMap clusterToPfoMap;

    ClusterList trackClusters3D;
    this->GetThreeDClusters(pAlgorithm, m_trackPfoListName, trackClusters3D, clusterToPfoMap);

    ClusterList showerClusters3D;
    this->GetThreeDClusters(pAlgorithm, m_showerPfoListName, showerClusters3D, clusterToPfoMap);

    ClusterSliceList clusterSliceList;
    this->GetClusterSliceList(trackClusters3D, showerClusters3D, clusterSliceList);

    if (clusterSliceList.size() < 2)
    {
        this->CopyAllHitsToSingleSlice(pAlgorithm, caloHitListNames, sliceList);
    }
    else
    {
        ClusterToSliceIndexMap clusterToSliceIndexMap;        
        this->CreateSlices(clusterSliceList, sliceList, clusterToSliceIndexMap);

        ClusterList assignedClusters;
        this->CopyPfoHitsToSlices(clusterSliceList, clusterToSliceIndexMap, clusterToPfoMap, sliceList, assignedClusters);

        ClusterList remainingClusters;
        this->GetRemainingClusters(pAlgorithm, clusterListNames, assignedClusters, remainingClusters);

        this->AssignRemainingHitsToSlices(remainingClusters, clusterToSliceIndexMap, sliceList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::GetThreeDClusters(const Algorithm *const /*pAlgorithm*/, const std::string &/*pfoListName*/, ClusterList &/*clusters3D*/,
    ClusterToPfoMap &/*clusterToPfoMap*/) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::GetClusterSliceList(const ClusterList &/*trackClusters3D*/, const ClusterList &/*showerClusters3D*/,
    ClusterSliceList &/*clusterSliceList*/) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::CopyAllHitsToSingleSlice(const Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
    SliceList &sliceList) const
{
    if (!sliceList.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const CaloHitList *pCaloHitListU(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm,
        caloHitListNames.at(TPC_VIEW_U), pCaloHitListU));

    const CaloHitList *pCaloHitListV(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm,
        caloHitListNames.at(TPC_VIEW_V), pCaloHitListV));

    const CaloHitList *pCaloHitListW(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm,
        caloHitListNames.at(TPC_VIEW_W), pCaloHitListW));

    if (pCaloHitListU || pCaloHitListV || pCaloHitListW)
    {
        sliceList.push_back(NeutrinoParentAlgorithm::Slice());
        NeutrinoParentAlgorithm::Slice &slice(sliceList.at(0));

        if (pCaloHitListU) slice.m_caloHitListU = *pCaloHitListU;
        if (pCaloHitListV) slice.m_caloHitListV = *pCaloHitListV;
        if (pCaloHitListW) slice.m_caloHitListW = *pCaloHitListW;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::CreateSlices(const ClusterSliceList &/*clusterSliceList*/, SliceList &/*sliceList*/, ClusterToSliceIndexMap &/*clusterToSliceIndexMap*/) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::CopyPfoHitsToSlices(const ClusterSliceList &/*clusterSliceList*/, const ClusterToSliceIndexMap &/*clusterToSliceIndexMap*/,
    const ClusterToPfoMap &/*clusterToPfoMap*/, SliceList &/*sliceList*/, ClusterList &/*assignedClusters*/) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::GetRemainingClusters(const Algorithm *const /*pAlgorithm*/, const HitTypeToNameMap &/*clusterListNames*/,
    const ClusterList &/*assignedClusters*/, ClusterList &/*remainingClusters*/) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::AssignRemainingHitsToSlices(const ClusterList &/*remainingClusters*/, const ClusterToSliceIndexMap &/*clusterToSliceIndexMap*/,
    SliceList &/*sliceList*/) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventSlicingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackPfoListName", m_trackPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerPfoListName", m_showerPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
