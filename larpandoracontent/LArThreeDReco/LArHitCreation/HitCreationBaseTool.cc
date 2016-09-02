/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.cc
 *
 *  @brief  Implementation of the hit creation base tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArPlugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

using namespace pandora;

namespace lar_content
{

HitCreationBaseTool::HitCreationBaseTool() :
    m_useChiSquaredApproach(true),
    m_useDeltaXCorrection(true),
    m_sigmaX(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitCreationBaseTool::~HitCreationBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationBaseTool::GetBestPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType1, const HitType hitType2,
    const CartesianPointList &fitPositionList1, const CartesianPointList &fitPositionList2, CartesianVector &position3D, float &chiSquared) const
{
    if (fitPositionList1.empty() && fitPositionList2.empty())
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    else if (fitPositionList1.empty())
    {
        if (fitPositionList2.size() != 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const CartesianVector &fitPosition2 = *(fitPositionList2.begin());
        this->GetPosition3D(pCaloHit2D, hitType2, fitPosition2, position3D, chiSquared);
    }
    else if (fitPositionList2.empty())
    {
        if (fitPositionList1.size() != 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const CartesianVector &fitPosition1 = *(fitPositionList1.begin());
        this->GetPosition3D(pCaloHit2D, hitType1, fitPosition1, position3D, chiSquared);
    }
    else
    {
        chiSquared = std::numeric_limits<float>::max();

        for (CartesianPointList::const_iterator iter1 = fitPositionList1.begin(), iterEnd1 = fitPositionList1.end(); iter1 != iterEnd1; ++iter1)
        {
            const CartesianVector &fitPosition1 = *iter1;
            for (CartesianPointList::const_iterator iter2 = fitPositionList2.begin(), iterEnd2 = fitPositionList2.end(); iter2 != iterEnd2; ++iter2)
            {
                const CartesianVector &fitPosition2 = *iter2;

                CartesianVector thisPosition3D(0.f, 0.f, 0.f);
                float thisChiSquared(std::numeric_limits<float>::max());
                this->GetPosition3D(pCaloHit2D, hitType1, hitType2, fitPosition1, fitPosition2, thisPosition3D, thisChiSquared);

                if (thisChiSquared < chiSquared)
                {
                    chiSquared = thisChiSquared;
                    position3D = thisPosition3D;
                }
            }
        }
    }

    // TODO Add additional chi-squared at edge of detector
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationBaseTool::GetPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType1, const HitType hitType2,
    const CartesianVector &fitPosition1, const CartesianVector &fitPosition2, CartesianVector &position3D, float &chiSquared) const
{
    // TODO Input better uncertainties into this method (sigmaHit, sigmaFit, sigmaX)
    const float sigmaHit(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetSigmaUVW());
    const float sigmaFit(sigmaHit); 
    const HitType hitType(pCaloHit2D->GetHitType());

    if (m_useChiSquaredApproach)
    {
        const double u((TPC_VIEW_U == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_U == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
        const double v((TPC_VIEW_V == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_V == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
        const double w((TPC_VIEW_W == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_W == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());

        const double sigmaU((TPC_VIEW_U == hitType) ? sigmaHit : sigmaFit);
        const double sigmaV((TPC_VIEW_V == hitType) ? sigmaHit : sigmaFit);
        const double sigmaW((TPC_VIEW_W == hitType) ? sigmaHit : sigmaFit);

        double bestY(std::numeric_limits<double>::max()), bestZ(std::numeric_limits<double>::max()), chiSquaredD(std::numeric_limits<double>::max());
        LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetMinChiSquaredYZ(u, v, w, sigmaU, sigmaV, sigmaW, bestY, bestZ, chiSquaredD);
        position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), static_cast<float>(bestY), static_cast<float>(bestZ));
        chiSquared = static_cast<float>(chiSquaredD);
    }
    else
    {
        const LArTransformationPlugin::PositionAndType hitPositionAndType(pCaloHit2D->GetPositionVector().GetZ(), hitType);
        const LArTransformationPlugin::PositionAndType fitPositionAndType1(fitPosition1.GetZ(), hitType1);
        const LArTransformationPlugin::PositionAndType fitPositionAndType2(fitPosition2.GetZ(), hitType2);

        double bestY(std::numeric_limits<double>::max()), bestZ(std::numeric_limits<double>::max()), chiSquaredD(std::numeric_limits<double>::max());
        LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetProjectedYZ(hitPositionAndType, fitPositionAndType1, fitPositionAndType2, sigmaHit, sigmaFit, bestY, bestZ, chiSquaredD);
        position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), static_cast<float>(bestY), static_cast<float>(bestZ));
        chiSquared = static_cast<float>(chiSquaredD);
    }

    if (m_useDeltaXCorrection)
    {
        const float deltaX1(pCaloHit2D->GetPositionVector().GetX() - fitPosition1.GetX());
        const float deltaX2(pCaloHit2D->GetPositionVector().GetX() - fitPosition2.GetX());
        chiSquared += ((deltaX1 * deltaX1) / (m_sigmaX * m_sigmaX)) + ((deltaX2 * deltaX2) / (m_sigmaX * m_sigmaX));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationBaseTool::GetPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType, const CartesianVector &fitPosition,
    CartesianVector &position3D, float &chiSquared) const
{
    if (pCaloHit2D->GetHitType() == hitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), pCaloHit2D->GetHitType(), hitType, pCaloHit2D->GetPositionVector(),
        fitPosition, position3D, chiSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitCreationBaseTool::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseChiSquaredApproach", m_useChiSquaredApproach));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDeltaXCorrection", m_useDeltaXCorrection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SigmaX", m_sigmaX));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
