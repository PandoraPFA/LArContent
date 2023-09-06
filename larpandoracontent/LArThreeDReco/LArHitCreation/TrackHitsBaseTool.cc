/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.cc
 *
 *  @brief  Implementation of the track hits base tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h"

using namespace pandora;

namespace lar_content
{

TrackHitsBaseTool::TrackHitsBaseTool() : m_minViews(2), m_slidingFitWindow(20)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackHitsBaseTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfo,
    const CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    try
    {
        if (!LArPfoHelper::IsTrack(pPfo))
            return;

        MatchedSlidingFitMap matchedSlidingFitMap;
        this->BuildSlidingFitMap(pPfo, matchedSlidingFitMap);

        if (matchedSlidingFitMap.size() < 2)
            return;

        this->GetTrackHits3D(inputTwoDHits, matchedSlidingFitMap, protoHitVector);
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackHitsBaseTool::BuildSlidingFitMap(const ParticleFlowObject *const pPfo, MatchedSlidingFitMap &matchedSlidingFitMap) const
{
    const ClusterList &pfoClusterList(pPfo->GetClusterList());

    ClusterVector pfoClusterVector;
    pfoClusterVector.insert(pfoClusterVector.end(), pfoClusterList.begin(), pfoClusterList.end());
    std::sort(pfoClusterVector.begin(), pfoClusterVector.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : pfoClusterVector)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (TPC_3D == hitType)
            continue;

        if (matchedSlidingFitMap.end() != matchedSlidingFitMap.find(hitType))
            continue;

        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinViews", m_minViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    return HitCreationBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
