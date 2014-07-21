/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/TransverseTrackHitCreationTool.cc
 * 
 *  @brief  Implementation of the transverse track hit creation tool.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCalculators/LArTransformationCalculator.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArThreeDReco/LArHitCreation/TransverseTrackHitCreationTool.h"

using namespace pandora;

namespace lar
{

void TransverseTrackHitCreationTool::Run(ThreeDHitCreationAlgorithm *pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitList &inputTwoDHits,
    CaloHitList &newThreeDHits, CaloHitList &omittedTwoDHits)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    try
    {
        Cluster *pClusterU(NULL), *pClusterV(NULL), *pClusterW(NULL);
        this->GetClusters(pPfo, pClusterU, pClusterV, pClusterW);

        TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
        LArClusterHelper::LArTwoDSlidingFit(pClusterU, m_slidingFitWindow, slidingFitResultU);
        LArClusterHelper::LArTwoDSlidingFit(pClusterV, m_slidingFitWindow, slidingFitResultV);
        LArClusterHelper::LArTwoDSlidingFit(pClusterW, m_slidingFitWindow, slidingFitResultW);

        CaloHitList caloHitListU, caloHitListV, caloHitListW;
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_U, caloHitListU);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_V, caloHitListV);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_W, caloHitListW);

        this->CreateThreeDHits(pAlgorithm, caloHitListU, slidingFitResultV, slidingFitResultW, newThreeDHits, omittedTwoDHits);
        this->CreateThreeDHits(pAlgorithm, caloHitListV, slidingFitResultU, slidingFitResultW, newThreeDHits, omittedTwoDHits);
        this->CreateThreeDHits(pAlgorithm, caloHitListW, slidingFitResultU, slidingFitResultV, newThreeDHits, omittedTwoDHits);
    }
    catch (StatusCodeException &)
    {
        std::cout << "TransverseTrackHitCreationTool: Unable to create 3D hits for provided cluster " << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseTrackHitCreationTool::GetClusters(const ParticleFlowObject *const pPfo, Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW) const
{
    pClusterU = NULL; pClusterV = NULL; pClusterW = NULL;
    const ClusterList &clusterList(pPfo->GetClusterList());

    if (3 != clusterList.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        Cluster *&pTargetCluster((TPC_VIEW_U == hitType) ? pClusterU : (TPC_VIEW_V == hitType) ? pClusterV : pClusterW);
        pTargetCluster = pCluster;
    }

    if ((NULL == pClusterU) || (NULL == pClusterV) || (NULL == pClusterW))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseTrackHitCreationTool::CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const CaloHitList &inputTwoDHits,
    const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2, CaloHitList &newThreeDHits, CaloHitList &omittedTwoDHits) const
{
    for (CaloHitList::const_iterator iter = inputTwoDHits.begin(), iterEnd = inputTwoDHits.end(); iter != iterEnd; ++iter)
    {
        try
        {
            CaloHit *pCaloHit2D(*iter);

            CartesianVector position3D(0.f, 0.f, 0.f);
            float chiSquared(std::numeric_limits<float>::max());
            this->GetPosition3D(pCaloHit2D, fitResult1, fitResult2, position3D, chiSquared);

            if (chiSquared > m_chiSquaredCut)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

            CaloHit *pCaloHit3D(NULL);
            pAlgorithm->CreateThreeDHit(pCaloHit2D, position3D, pCaloHit3D);
            newThreeDHits.insert(pCaloHit3D);
        }
        catch (StatusCodeException &)
        {
            omittedTwoDHits.insert(*iter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseTrackHitCreationTool::GetPosition3D(const CaloHit *const pCaloHit2D, const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    CartesianVector &position3D, float &chiSquared) const
{
    CartesianPointList fitPositions1, fitPositions2;
    fitResult1.GetGlobalFitPositionListAtX(pCaloHit2D->GetPositionVector().GetX(), fitPositions1);
    fitResult2.GetGlobalFitPositionListAtX(pCaloHit2D->GetPositionVector().GetX(), fitPositions2);

    const HitType hitType1(LArClusterHelper::GetClusterHitType(fitResult1.GetCluster()));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(fitResult2.GetCluster()));

    bool foundPosition(false);
    CartesianVector bestPosition3D(0.f, 0.f, 0.f);
    float bestChiSquared(std::numeric_limits<float>::max());

    for (CartesianPointList::const_iterator iter1 = fitPositions1.begin(), iter1End = fitPositions1.end(); iter1 != iter1End; ++iter1)
    {
        for (CartesianPointList::const_iterator iter2 = fitPositions2.begin(), iter2End = fitPositions2.end(); iter2 != iter2End; ++iter2)
        {
            CartesianVector thisPosition3D(0.f, 0.f, 0.f);
            float thisChiSquared(std::numeric_limits<float>::max());
            this->GetPosition3D(pCaloHit2D, hitType1, hitType2, *iter1, *iter2, thisPosition3D, thisChiSquared);

            if (thisChiSquared < bestChiSquared)
            {
                foundPosition = true;
                bestPosition3D = thisPosition3D;
                bestChiSquared = thisChiSquared;
            }
        }
    }

    if (!foundPosition)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    position3D = bestPosition3D;
    chiSquared = bestChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseTrackHitCreationTool::GetPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType1, const HitType hitType2,
    const CartesianVector &fitPosition1, const CartesianVector &fitPosition2, CartesianVector &position3D, float &chiSquared) const
{
    const float sigmaHit(LArGeometryHelper::GetLArTransformationCalculator()->GetSigmaUVW());
    const float sigmaFit(m_sigmaFitMultiplier * sigmaHit);
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

StatusCode TransverseTrackHitCreationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    m_useChiSquaredApproach = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseChiSquaredApproach", m_useChiSquaredApproach));

    m_sigmaFitMultiplier = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SigmaFitMultiplier", m_sigmaFitMultiplier));

    m_chiSquaredCut = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
