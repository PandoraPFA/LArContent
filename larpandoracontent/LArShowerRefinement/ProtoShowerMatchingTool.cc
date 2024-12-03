/**
 *  @file   larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.cc
 *
 *  @brief  Implementation of the ProtoShower matching tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h"

using namespace pandora;

namespace lar_content
{

ProtoShowerMatchingTool::ProtoShowerMatchingTool() :
    m_spineSlidingFitWindow(20),
    m_maxXSeparation(5.f),
    m_maxSeparation(5.f),
    m_xExtrapolation(5.f),
    m_maxAngularDeviation(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoShowerMatchingTool::Run(const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV,
    const ProtoShowerVector &protoShowerVectorW, ProtoShowerMatchVector &protoShowerMatchVector)
{
    IntVector usedProtoShowersU, usedProtoShowersV, usedProtoShowersW;

    bool added(false);

    for (unsigned int uIndex = 0; uIndex < protoShowerVectorU.size(); ++uIndex)
    {
        added = false;

        if (std::find(usedProtoShowersU.begin(), usedProtoShowersU.end(), uIndex) != usedProtoShowersU.end())
            continue;

        const ProtoShower &protoShowerU(protoShowerVectorU.at(uIndex));

        for (unsigned int vIndex = 0; vIndex < protoShowerVectorV.size(); ++vIndex)
        {
            if (std::find(usedProtoShowersV.begin(), usedProtoShowersV.end(), vIndex) != usedProtoShowersV.end())
                continue;

            const ProtoShower &protoShowerV(protoShowerVectorV.at(vIndex));

            for (unsigned int wIndex = 0; wIndex < protoShowerVectorW.size(); ++wIndex)
            {
                if (std::find(usedProtoShowersW.begin(), usedProtoShowersW.end(), wIndex) != usedProtoShowersW.end())
                    continue;

                const ProtoShower &protoShowerW(protoShowerVectorW.at(wIndex));
                Consistency consistency(Consistency::POSITION);

                if (!this->ArePathwaysConsistent(protoShowerU, protoShowerV, protoShowerW, consistency))
                    continue;

                usedProtoShowersU.push_back(uIndex);
                usedProtoShowersV.push_back(vIndex);
                usedProtoShowersW.push_back(wIndex);

                protoShowerMatchVector.push_back(ProtoShowerMatch(protoShowerU, protoShowerV, protoShowerW, consistency));

                added = true;
                break;
            }

            if (added)
                break;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::ArePathwaysConsistent(
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, Consistency &consistency) const
{
    if (this->AreShowerStartsConsistent(protoShowerU, protoShowerV, protoShowerW))
    {
        consistency = Consistency::POSITION;
    }
    else if (this->AreDirectionsConsistent(protoShowerU, protoShowerV, protoShowerW))
    {
        consistency = Consistency::DIRECTION;
    }
    else
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::AreShowerStartsConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW) const
{
    const CartesianVector &showerStartU(protoShowerU.GetShowerCore().GetStartPosition());
    const CartesianVector &showerStartV(protoShowerV.GetShowerCore().GetStartPosition());
    const CartesianVector &showerStartW(protoShowerW.GetShowerCore().GetStartPosition());

    const float xSeparationUV(std::fabs(showerStartU.GetX() - showerStartV.GetX()));
    const float xSeparationUW(std::fabs(showerStartU.GetX() - showerStartW.GetX()));
    const float xSeparationVW(std::fabs(showerStartV.GetX() - showerStartW.GetX()));

    if ((xSeparationUV > m_maxXSeparation) || (xSeparationUW > m_maxXSeparation) || (xSeparationVW > m_maxXSeparation))
        return false;

    float chi2(0.f);
    CartesianVector projectionU(0.f, 0.f, 0.f), projectionV(0.f, 0.f, 0.f), projectionW(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, showerStartV, showerStartW, projectionU, chi2);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, showerStartW, showerStartU, projectionV, chi2);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, showerStartU, showerStartV, projectionW, chi2);

    const float separationU((projectionU - showerStartU).GetMagnitude());
    const float separationV((projectionV - showerStartV).GetMagnitude());
    const float separationW((projectionW - showerStartW).GetMagnitude());

    const float metric((separationU + separationV + separationW) / 3.f);

    return (metric < m_maxSeparation);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::AreDirectionsConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW) const
{
    const CartesianVector &directionU1(protoShowerU.GetConnectionPathway().GetStartDirection());
    const CartesianVector &directionV1(protoShowerV.GetConnectionPathway().GetStartDirection());
    const CartesianVector &directionW1(protoShowerW.GetConnectionPathway().GetStartDirection());

    if (this->AreDirectionsConsistent(protoShowerU.GetConnectionPathway().GetStartPosition(), protoShowerV.GetConnectionPathway().GetStartPosition(), 
        protoShowerW.GetConnectionPathway().GetStartPosition(), directionU1, directionV1, directionW1))
    {
        return true;
    }
    else
    {
        const bool isDownstream(
            protoShowerW.GetShowerCore().GetStartPosition().GetZ() > protoShowerW.GetConnectionPathway().GetStartPosition().GetZ());

        CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

        for (const CaloHit *const pCaloHit : protoShowerU.GetSpineHitList())
            spinePositionsU.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerV.GetSpineHitList())
            spinePositionsV.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerW.GetSpineHitList())
            spinePositionsW.push_back(pCaloHit->GetPositionVector());

        const TwoDSlidingFitResult spineFitU(&spinePositionsU, m_spineSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, m_spineSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, m_spineSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W));

        const CartesianVector directionU2(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionV2(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionW2(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

        return this->AreDirectionsConsistent(protoShowerU.GetConnectionPathway().GetStartPosition(), protoShowerV.GetConnectionPathway().GetStartPosition(), 
            protoShowerW.GetConnectionPathway().GetStartPosition(), directionU2, directionV2, directionW2);
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::AreDirectionsConsistent(const CartesianVector &nuVertexU, const CartesianVector &nuVertexV, const CartesianVector &nuVertexW,
    const CartesianVector &directionU, const CartesianVector &directionV, const CartesianVector &directionW) const
{
    const CartesianVector wireAxis(0.f, 0.f, 1.f);

    float wireDeviationU(directionU.GetOpeningAngle(wireAxis));
    wireDeviationU = std::min(wireDeviationU, static_cast<float>(M_PI - wireDeviationU));

    float wireDeviationV(directionV.GetOpeningAngle(wireAxis));
    wireDeviationV = std::min(wireDeviationV, static_cast<float>(M_PI - wireDeviationV));

    float wireDeviationW(directionW.GetOpeningAngle(wireAxis));
    wireDeviationW = std::min(wireDeviationW, static_cast<float>(M_PI - wireDeviationW));

    float radians((2.f * M_PI) / 180.f);
    bool isIsochronous((wireDeviationU < radians) && (wireDeviationV < radians) && (wireDeviationW < radians));

    if (isIsochronous)
        return true;

    if (directionU.GetX() * directionV.GetX() < 0.f)
        return false;

    if (directionU.GetX() * directionW.GetX() < 0.f)
        return false;

    if (directionV.GetX() * directionW.GetX() < 0.f)
        return false;

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float cosOpeningAngleU(std::fabs(directionU.GetCosOpeningAngle(xAxis)));
    const float cosOpeningAngleV(std::fabs(directionV.GetCosOpeningAngle(xAxis)));
    const float cosOpeningAngleW(std::fabs(directionW.GetCosOpeningAngle(xAxis)));
    const CartesianVector uPoint(nuVertexU + (directionU * (m_xExtrapolation / cosOpeningAngleU)));
    const CartesianVector vPoint(nuVertexV + (directionV * (m_xExtrapolation / cosOpeningAngleV)));
    const CartesianVector wPoint(nuVertexW + (directionW * (m_xExtrapolation / cosOpeningAngleW)));

    float chiSquared(0.f);

    CartesianVector projectionPointU(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vPoint, wPoint, projectionPointU, chiSquared);

    CartesianVector projectionPointV(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, wPoint, uPoint, projectionPointV, chiSquared);

    CartesianVector projectionPointW(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, uPoint, vPoint, projectionPointW, chiSquared);

    // Project U
    const CartesianVector projectionU((projectionPointU - nuVertexU).GetUnitVector());

    // Project V
    const CartesianVector projectionV((projectionPointV - nuVertexV).GetUnitVector());

    // Project W
    const CartesianVector projectionW((projectionPointW - nuVertexW).GetUnitVector());

    float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180.f / M_PI);
    float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180.f / M_PI);
    float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180.f / M_PI);

    const float metric((openingAngleU + openingAngleV + openingAngleW) / 3.f);

    if (metric > m_maxAngularDeviation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoShowerMatchingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SpineSlidingFitWindow", m_spineSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSeparation", m_maxSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "XExtrapolation", m_xExtrapolation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAngularDeviation", m_maxAngularDeviation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
