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

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm()
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
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    const int nHits(pCluster->GetNCaloHits());

    if (nHits < 20) // TODO config
        return false;

    float hitsPerUnitLength(-std::numeric_limits<float>::max());
    float widthPerUnitLength(-std::numeric_limits<float>::max());

    try
    {
        const float slidingFitPitch(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, 10, slidingFitPitch); // TODO config

        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        if (layerFitResultMapS.empty())
            return false;

        CartesianVector globalMinLayerPositionOnAxis(0.f, 0.f, 0.f), globalMaxLayerPositionOnAxis(0.f, 0.f, 0.f);
        showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, globalMinLayerPositionOnAxis);
        showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, globalMaxLayerPositionOnAxis);
        const float straightLinePathLength((globalMaxLayerPositionOnAxis - globalMinLayerPositionOnAxis).GetMagnitude());

        float widthSum(0.f);

        for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
        {
            LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
            LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

            if ((layerFitResultMapP.end() == iterP) || (layerFitResultMapN.end() == iterN))
                continue;

            widthSum += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
        }

        if (straightLinePathLength > std::numeric_limits<float>::epsilon())
        {
            hitsPerUnitLength = static_cast<float>(nHits) / straightLinePathLength;
            widthPerUnitLength = widthSum / straightLinePathLength;
        }
    }
    catch (StatusCodeException &)
    {
        return false;
    }

    if ((widthPerUnitLength > 0.f) && ((widthPerUnitLength < 0.15f) ||
        ((nHits > 40) && (hitsPerUnitLength < 1.5f)) ||
        ((nHits > 60) && (widthPerUnitLength < 0.25f)) )) // TODO config + refine logic
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
