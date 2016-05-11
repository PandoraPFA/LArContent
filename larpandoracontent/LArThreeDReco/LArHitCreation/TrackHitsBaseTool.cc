/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.cc
 *
 *  @brief  Implementation of the track hits base tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArPlugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h"

using namespace pandora;

namespace lar_content
{

TrackHitsBaseTool::TrackHitsBaseTool() :
    m_minViews(2),
    m_slidingFitWindow(20),
    m_chiSquaredCut(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackHitsBaseTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitList &inputTwoDHits,
    CaloHitList &newThreeDHits)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    try
    {
        if (!LArPfoHelper::IsTrack(pPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        MatchedSlidingFitMap matchedSlidingFitMap;
        this->BuildSlidingFitMap(pPfo, matchedSlidingFitMap);

        if (matchedSlidingFitMap.size() < 2)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        this->CreateThreeDHits(pAlgorithm, inputTwoDHits, matchedSlidingFitMap, newThreeDHits);
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackHitsBaseTool::BuildSlidingFitMap(const ParticleFlowObject *const pPfo, MatchedSlidingFitMap &matchedSlidingFitMap) const
{
    const ClusterList &pfoClusterList(pPfo->GetClusterList());
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    ClusterVector pfoClusterVector;
    pfoClusterVector.insert(pfoClusterVector.end(), pfoClusterList.begin(), pfoClusterList.end());
    std::sort(pfoClusterVector.begin(), pfoClusterVector.end(), LArClusterHelper::SortByNHits);

    for (ClusterVector::const_iterator iter = pfoClusterVector.begin(), iterEnd = pfoClusterVector.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster(*iter);
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (TPC_3D == hitType)
            continue;

        if (matchedSlidingFitMap.end() != matchedSlidingFitMap.find(hitType))
            continue;

        try
        {
            const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

            if (!matchedSlidingFitMap.insert(MatchedSlidingFitMap::value_type(hitType, slidingFitResult)).second)
                throw StatusCodeException(STATUS_CODE_FAILURE);
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackHitsBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinViews", m_minViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    return HitCreationBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
