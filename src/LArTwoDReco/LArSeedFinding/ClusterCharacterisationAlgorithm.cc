/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"

#include "LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "LArPlugins/LArTransformationPlugin.h"

#include "LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm() :
    m_slidingFitWindow(10),
    m_minHitsInCluster(20),
    m_maxLayerGapFraction(0.2f),
    m_maxWidthPerUnitLength(0.15f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    for (StringVector::const_iterator listIter = m_inputClusterListNames.begin(), listIterEnd = m_inputClusterListNames.end(); listIter != listIterEnd; ++listIter)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *listIter, pClusterList));

        for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        {
            Cluster *const pCluster(*iter);

            if (this->IsClearTrack(pCluster))
            {
                pCluster->SetIsFixedMuonFlag(true);
            }
            else
            {
                pCluster->SetIsFixedMuonFlag(false);
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minHitsInCluster)
        return false;

    try
    {
        const float slidingFitPitch(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        if (layerFitResultMapS.size() < 2)
            return false;

        CartesianVector globalMinLayerPositionOnAxis(0.f, 0.f, 0.f), globalMaxLayerPositionOnAxis(0.f, 0.f, 0.f);
        showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, globalMinLayerPositionOnAxis);
        showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, globalMaxLayerPositionOnAxis);
        const float straightLinePathLength((globalMaxLayerPositionOnAxis - globalMinLayerPositionOnAxis).GetMagnitude());

        if (straightLinePathLength < std::numeric_limits<float>::epsilon())
            return false;

        float widthSum(0.f);
        int biggestLayerGap(0), previousLayer(layerFitResultMapS.begin()->first);

        for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
        {
            const int layerGap(iterS->first - previousLayer);
            previousLayer = iterS->first;

            if (layerGap > biggestLayerGap)
                biggestLayerGap = layerGap;

            LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
            LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

            if ((layerFitResultMapP.end() == iterP) || (layerFitResultMapN.end() == iterN))
                continue;

            widthSum += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
        }

        const float layerGapFraction(static_cast<float>(biggestLayerGap) / static_cast<float>(layerFitResultMapS.size()));

        if (layerGapFraction > m_maxLayerGapFraction)
            return false;

        // ATTN: Could also try to recover long (sideways) tracks, with cuts on nHits and nHitsPerUnitLength
        const float widthPerUnitLength(widthSum / straightLinePathLength);

        if (widthPerUnitLength < m_maxWidthPerUnitLength)
            return true;
    }
    catch (StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsInCluster", m_minHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLayerGapFraction", m_maxLayerGapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxWidthPerUnitLength", m_maxWidthPerUnitLength));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
