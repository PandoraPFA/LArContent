/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.cc
 *
 *  @brief  Implementation of the transverse track hit creation tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCalculators/LArTransformationCalculator.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h"

using namespace pandora;

namespace lar
{

void TrackHitsBaseTool::Run(ThreeDHitCreationAlgorithm *pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitList &inputTwoDHits,
    CaloHitList &newThreeDHits, CaloHitList &omittedTwoDHits)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    try
    {
        MatchedSlidingFitMap matchedSlidingFitMap;
        this->BuildSlidingFitMap(pPfo, matchedSlidingFitMap);

        if (matchedSlidingFitMap.size() < 2)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        this->CreateThreeDHits(pAlgorithm, inputTwoDHits, matchedSlidingFitMap, newThreeDHits, omittedTwoDHits);
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

void TrackHitsBaseTool::GetPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType1, const HitType hitType2,
    const CartesianVector &fitPosition1, const CartesianVector &fitPosition2, CartesianVector &position3D, float &chiSquared) const
{
    const float sigmaHit(LArGeometryHelper::GetLArTransformationCalculator()->GetSigmaUVW());
    const float sigmaFit(sigmaHit); // TODO: Input uncertainties into this method
    const HitType hitType(pCaloHit2D->GetHitType());

    if (m_useChiSquaredApproach)
    {
        const float u((TPC_VIEW_U == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_U == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
        const float v((TPC_VIEW_V == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_V == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
        const float w((TPC_VIEW_W == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_W == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());

        const float sigmaU((TPC_VIEW_U == hitType) ? sigmaHit : sigmaFit);
        const float sigmaV((TPC_VIEW_V == hitType) ? sigmaHit : sigmaFit);
        const float sigmaW((TPC_VIEW_W == hitType) ? sigmaHit : sigmaFit);

        float bestY(std::numeric_limits<float>::max()), bestZ(std::numeric_limits<float>::max());
        LArGeometryHelper::GetLArTransformationCalculator()->GetMinChiSquaredYZ(u, v, w, sigmaU, sigmaV, sigmaW, bestY, bestZ, chiSquared);
        position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), bestY, bestZ);
    }
    else
    {
        const LArTransformationCalculator::PositionAndType hitPositionAndType(pCaloHit2D->GetPositionVector().GetZ(), hitType);
        const LArTransformationCalculator::PositionAndType fitPositionAndType1(fitPosition1.GetZ(), hitType1);
        const LArTransformationCalculator::PositionAndType fitPositionAndType2(fitPosition2.GetZ(), hitType2);

        float bestY(std::numeric_limits<float>::max()), bestZ(std::numeric_limits<float>::max());
        LArGeometryHelper::GetLArTransformationCalculator()->GetProjectedYZ(hitPositionAndType, fitPositionAndType1, fitPositionAndType2, sigmaHit, sigmaFit, bestY, bestZ, chiSquared);
        position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), bestY, bestZ);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackHitsBaseTool::GetPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType, const CartesianVector &fitPosition,
    CartesianVector &position3D, float &chiSquared) const
{
    if (pCaloHit2D->GetHitType() == hitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    LArGeometryHelper::MergeTwoPositions3D(pCaloHit2D->GetHitType(), hitType, pCaloHit2D->GetPositionVector(), fitPosition,
        position3D, chiSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackHitsBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    m_useChiSquaredApproach = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseChiSquaredApproach", m_useChiSquaredApproach));

    m_chiSquaredCut = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
