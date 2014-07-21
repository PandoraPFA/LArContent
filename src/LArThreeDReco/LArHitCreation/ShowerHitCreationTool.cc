/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ShowerHitCreationTool.cc
 * 
 *  @brief  Implementation of the shower hit creation tool.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCalculators/LArTransformationCalculator.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArHitCreation/ShowerHitCreationTool.h"

using namespace pandora;

namespace lar
{

void ShowerHitCreationTool::Run(ThreeDHitCreationAlgorithm *pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitList &inputTwoDHits,
    CaloHitList &newThreeDHits, CaloHitList &omittedTwoDHits)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    try
    {
        Cluster *pClusterU(NULL), *pClusterV(NULL), *pClusterW(NULL);
        this->GetClusters(pPfo, pClusterU, pClusterV, pClusterW);

        CaloHitList caloHitListU, caloHitListV, caloHitListW;
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_U, caloHitListU);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_V, caloHitListV);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_W, caloHitListW);

        this->CreateThreeDHits(pAlgorithm, caloHitListU, caloHitListV, caloHitListW, newThreeDHits, omittedTwoDHits);
        this->CreateThreeDHits(pAlgorithm, caloHitListV, caloHitListU, caloHitListW, newThreeDHits, omittedTwoDHits);
        this->CreateThreeDHits(pAlgorithm, caloHitListW, caloHitListU, caloHitListV, newThreeDHits, omittedTwoDHits);
    }
    catch (StatusCodeException &)
    {
        std::cout << "ShowerHitCreationTool: Unable to create 3D hits for provided cluster " << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitCreationTool::GetClusters(const ParticleFlowObject *const pPfo, Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW) const
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

void ShowerHitCreationTool::CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const CaloHitList &inputTwoDHits,
    const CaloHitList &caloHitList1, const CaloHitList &caloHitList2, CaloHitList &newThreeDHits, CaloHitList &omittedTwoDHits) const
{
    for (CaloHitList::const_iterator iter = inputTwoDHits.begin(), iterEnd = inputTwoDHits.end(); iter != iterEnd; ++iter)
    {
        try
        {
            CaloHit *pCaloHit2D(*iter);
            const float x(pCaloHit2D->GetPositionVector().GetX());

            CaloHitList filteredList1, filteredList2;
            this->FilterCaloHits(x, m_xTolerance, caloHitList1, filteredList1);
            this->FilterCaloHits(x, m_xTolerance, caloHitList2, filteredList2);

            CartesianVector position3D(0.f, 0.f, 0.f);
            float chiSquared(std::numeric_limits<float>::max());
            this->GetPosition3D(pCaloHit2D, filteredList1, filteredList2, position3D, chiSquared);

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

void ShowerHitCreationTool::FilterCaloHits(const float x, const float xTolerance, const CaloHitList &inputCaloHitList, CaloHitList &outputCaloHitList) const
{
    outputCaloHitList.clear();

    for (CaloHitList::const_iterator iter = inputCaloHitList.begin(), iterEnd = inputCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;
        const float deltaX(pCaloHit->GetPositionVector().GetX() - x);

        if (std::fabs(deltaX) < xTolerance)
            outputCaloHitList.insert(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitCreationTool::GetPosition3D(const CaloHit *const pCaloHit2D, const CaloHitList &caloHitList1, const CaloHitList &caloHitList2,
    CartesianVector &position3D, float &chiSquared) const
{
    bool foundPosition(false);
    CartesianVector bestPosition3D(0.f, 0.f, 0.f);
    float bestChiSquared(std::numeric_limits<float>::max());

    for (CaloHitList::const_iterator iter1 = caloHitList1.begin(), iter1End = caloHitList1.end(); iter1 != iter1End; ++iter1)
    {
        for (CaloHitList::const_iterator iter2 = caloHitList2.begin(), iter2End = caloHitList2.end(); iter2 != iter2End; ++iter2)
        {
            CartesianVector thisPosition3D(0.f, 0.f, 0.f);
            float thisChiSquared(std::numeric_limits<float>::max());
            this->GetPosition3D(pCaloHit2D, *iter1, *iter2, thisPosition3D, thisChiSquared);

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

void ShowerHitCreationTool::GetPosition3D(const CaloHit *const pCaloHit2D, const CaloHit *const pCaloHit1, const CaloHit *const pCaloHit2,
    CartesianVector &position3D, float &chiSquared) const
{
    const HitType hitType(pCaloHit2D->GetHitType());
    const HitType hitType1(pCaloHit1->GetHitType());
    const HitType hitType2(pCaloHit2->GetHitType());
    const float sigmaHit(LArGeometryHelper::GetLArTransformationCalculator()->GetSigmaUVW());

    if (m_useChiSquaredApproach)
    {
        const float u((TPC_VIEW_U == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_U == hitType1) ? pCaloHit1->GetPositionVector().GetZ() : pCaloHit2->GetPositionVector().GetZ());
        const float v((TPC_VIEW_V == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_V == hitType1) ? pCaloHit1->GetPositionVector().GetZ() : pCaloHit2->GetPositionVector().GetZ());
        const float w((TPC_VIEW_W == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_W == hitType1) ? pCaloHit1->GetPositionVector().GetZ() : pCaloHit2->GetPositionVector().GetZ());

        float bestY(std::numeric_limits<float>::max()), bestZ(std::numeric_limits<float>::max());
        LArGeometryHelper::GetLArTransformationCalculator()->GetMinChiSquaredYZ(u, v, w, sigmaHit, sigmaHit, sigmaHit, bestY, bestZ, chiSquared);
        position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), bestY, bestZ);
    }
    else
    {
        const LArTransformationCalculator::PositionAndType hitPositionAndType(pCaloHit2D->GetPositionVector().GetZ(), hitType);
        const LArTransformationCalculator::PositionAndType fitPositionAndType1(pCaloHit1->GetPositionVector().GetZ(), hitType1);
        const LArTransformationCalculator::PositionAndType fitPositionAndType2(pCaloHit2->GetPositionVector().GetZ(), hitType2);

        float bestY(std::numeric_limits<float>::max()), bestZ(std::numeric_limits<float>::max());
        LArGeometryHelper::GetLArTransformationCalculator()->GetProjectedYZ(hitPositionAndType, fitPositionAndType1, fitPositionAndType2, sigmaHit, sigmaHit, bestY, bestZ, chiSquared);
        position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), bestY, bestZ);
    }

    if (m_useDeltaXCorrection)
    {
        const float deltaX1(pCaloHit1->GetPositionVector().GetX() - pCaloHit2D->GetPositionVector().GetX());
        const float deltaX2(pCaloHit2->GetPositionVector().GetX() - pCaloHit2D->GetPositionVector().GetX());
        chiSquared += ((deltaX1 * deltaX1) / (m_sigmaX * m_sigmaX)) + ((deltaX2 * deltaX2) / (m_sigmaX * m_sigmaX));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerHitCreationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_xTolerance = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "XTolerance", m_xTolerance));

    m_useDeltaXCorrection = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDeltaXCorrection", m_useDeltaXCorrection));

    m_sigmaX = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SigmaX", m_sigmaX));

    m_useChiSquaredApproach = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseChiSquaredApproach", m_useChiSquaredApproach));

    m_chiSquaredCut = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
