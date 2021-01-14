/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.cc
 *
 *  @brief  Implementation of the hit creation base tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Geometry/LArTPC.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

using namespace pandora;

namespace lar_content
{

HitCreationBaseTool::HitCreationBaseTool() :
    m_sigmaX2(1.),
    m_sigmaYZ2(10.),
    m_chiSquaredCut(1.)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitCreationBaseTool::~HitCreationBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationBaseTool::GetBestPosition3D(const HitType hitType1, const HitType hitType2, const CartesianPointVector &fitPositionList1,
    const CartesianPointVector &fitPositionList2, ProtoHit &protoHit) const
{
    if (fitPositionList1.empty() && fitPositionList2.empty())
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    else if (fitPositionList1.empty())
    {
        if (fitPositionList2.size() != 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        this->GetBestPosition3D(hitType2, fitPositionList2.front(), protoHit);
    }
    else if (fitPositionList2.empty())
    {
        if (fitPositionList1.size() != 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        this->GetBestPosition3D(hitType1, fitPositionList1.front(), protoHit);
    }
    else
    {
        for (const CartesianVector &fitPosition1 : fitPositionList1)
        {
            for (const CartesianVector &fitPosition2 : fitPositionList2)
            {
                ProtoHit thisProtoHit(protoHit.GetParentCaloHit2D());
                this->GetBestPosition3D(hitType1, hitType2, fitPosition1, fitPosition2, thisProtoHit);

                if (!protoHit.IsPositionSet() || (thisProtoHit.GetChi2() < protoHit.GetChi2()))
                    protoHit = thisProtoHit;
            }
        }
    }

    // TODO Add additional chi-squared at edge of detector
}

//------------------------------------------------------------------------------------------------------------------------------------------

double HitCreationBaseTool::GetDistanceToDetectorEdge(const LArTPCMap &larTPCMap, const CartesianVector &position3D) const
{
    double distanceToEdge = 0;
    const double bestY(position3D.GetY());
    const double bestZ(position3D.GetZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC* currentTPC = mapEntry.second;
        // Check if the given hit is contained in the detector. If it is outside, then we get a positive value which
        // will set the distanceToEdge value.
        //
        // We want to know the largest distance outside of the detector, so that we can make these 3D hits less favourable.
        distanceToEdge = std::max(distanceToEdge, (currentTPC->GetCenterY() - 0.5f * currentTPC->GetWidthY()) - bestY);
        distanceToEdge = std::max(distanceToEdge, bestY - (currentTPC->GetCenterY() + 0.5f * currentTPC->GetWidthY()));
        distanceToEdge = std::max(distanceToEdge, (currentTPC->GetCenterZ() - 0.5f * currentTPC->GetWidthZ()) - bestZ);
        distanceToEdge = std::max(distanceToEdge, bestZ - (currentTPC->GetCenterZ() + 0.5f * currentTPC->GetWidthZ()));
    }

    return distanceToEdge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationBaseTool::GetBestPosition3D(const HitType hitType1, const HitType hitType2, const CartesianVector &fitPosition1,
    const CartesianVector &fitPosition2, ProtoHit &protoHit) const
{
    // TODO Input better uncertainties into this method (sigmaHit, sigmaFit, sigmaX)
    const CaloHit *const pCaloHit2D(protoHit.GetParentCaloHit2D());
    const HitType hitType(pCaloHit2D->GetHitType());

    const double sigmaFit(LArGeometryHelper::GetSigmaUVW(this->GetPandora()));
    const double sigmaHit(sigmaFit);

    CartesianVector position3D(0.f, 0.f, 0.f);
    double chi2(std::numeric_limits<double>::max());

    const double u((TPC_VIEW_U == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_U == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
    const double v((TPC_VIEW_V == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_V == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
    const double w((TPC_VIEW_W == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_W == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());

    const double sigmaU((TPC_VIEW_U == hitType) ? sigmaHit : sigmaFit);
    const double sigmaV((TPC_VIEW_V == hitType) ? sigmaHit : sigmaFit);
    const double sigmaW((TPC_VIEW_W == hitType) ? sigmaHit : sigmaFit);

    double bestY(std::numeric_limits<double>::max()), bestZ(std::numeric_limits<double>::max());
    this->GetPandora().GetPlugins()->GetLArTransformationPlugin()->GetMinChiSquaredYZ(u, v, w, sigmaU, sigmaV, sigmaW, bestY, bestZ, chi2);
    position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), static_cast<float>(bestY), static_cast<float>(bestZ));

    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    double distanceToEdge = this->GetDistanceToDetectorEdge(larTPCMap, position3D);

    const double deltaX1(pCaloHit2D->GetPositionVector().GetX() - fitPosition1.GetX());
    const double deltaX2(pCaloHit2D->GetPositionVector().GetX() - fitPosition2.GetX());
    const double chi2X(((deltaX1 * deltaX1) / m_sigmaX2) + ((deltaX2 * deltaX2) / m_sigmaX2));
    const double chi2YZ((distanceToEdge * distanceToEdge) / m_sigmaYZ2);

    protoHit.SetPosition3D(position3D, chi2 + chi2X + chi2YZ);
    protoHit.AddTrajectorySample(TrajectorySample(fitPosition1, hitType1, sigmaFit));
    protoHit.AddTrajectorySample(TrajectorySample(fitPosition2, hitType2, sigmaFit));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationBaseTool::GetBestPosition3D(const HitType hitType, const CartesianVector &fitPosition, ProtoHit &protoHit) const
{
    // TODO Input better uncertainties into this method (sigmaHit, sigmaFit, sigmaX)
    const CaloHit *const pCaloHit2D(protoHit.GetParentCaloHit2D());
    const double sigmaFit(LArGeometryHelper::GetSigmaUVW(this->GetPandora()));

    if (pCaloHit2D->GetHitType() == hitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector position3D(0.f, 0.f, 0.f);
    float chi2(std::numeric_limits<float>::max());
    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), pCaloHit2D->GetHitType(), hitType, pCaloHit2D->GetPositionVector(),
        fitPosition, position3D, chi2);

    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    double distanceToEdge = this->GetDistanceToDetectorEdge(larTPCMap, position3D);

    // ATTN Replace chi2 from LArGeometryHelper for consistency with three-view treatment (purely a measure of delta x)
    const double deltaX(pCaloHit2D->GetPositionVector().GetX() - fitPosition.GetX());
    const double chi2X((deltaX * deltaX) / m_sigmaX2);
    const double chi2YZ((distanceToEdge * distanceToEdge) / m_sigmaYZ2);

    protoHit.SetPosition3D(position3D, chi2 + chi2X + chi2YZ);
    protoHit.AddTrajectorySample(TrajectorySample(fitPosition, hitType, sigmaFit));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitCreationBaseTool::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    double sigmaX(std::sqrt(m_sigmaX2));
    double sigmaYZ(std::sqrt(m_sigmaYZ2));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SigmaX", sigmaX));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SigmaYZ", sigmaYZ));

    m_sigmaX2 = sigmaX * sigmaX;
    m_sigmaYZ2 = sigmaYZ * sigmaYZ;

    if (m_sigmaX2 < std::numeric_limits<double>::epsilon())
    {
        std::cout << "HitCreationBaseTool - Invalid parameter, SigmaX: " << sigmaX << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_sigmaYZ2 < std::numeric_limits<double>::epsilon())
    {
        std::cout << "HitCreationBaseTool - Invalid parameter, SigmaYZ: " << sigmaYZ << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
