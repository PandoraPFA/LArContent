/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.cc
 *
 *  @brief  Implementation of the track hits base tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPfoHelper.h"
#include "LArHelpers/LArClusterHelper.h"

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h"

using namespace pandora;

namespace lar
{

void TrackHitsBaseTool::Run(ThreeDHitCreationAlgorithm *pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitList &inputTwoDHits,
    CaloHitList &newThreeDHits)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

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
    const ClusterList &clusterList(pPfo->GetClusterList());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (TPC_3D == hitType)
            continue;

        if (matchedSlidingFitMap.end() != matchedSlidingFitMap.find(hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        TwoDSlidingFitResult slidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(pCluster, m_slidingFitWindow, slidingFitResult);

        if (!matchedSlidingFitMap.insert(MatchedSlidingFitMap::value_type(hitType, slidingFitResult)).second)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackHitsBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minViews = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinViews", m_minViews));

    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    m_chiSquaredCut = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    return HitCreationBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar
