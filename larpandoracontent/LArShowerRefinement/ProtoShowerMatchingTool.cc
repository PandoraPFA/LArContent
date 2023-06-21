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

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h"

using namespace pandora;

namespace lar_content
{

ProtoShowerMatchingTool::ProtoShowerMatchingTool() :
    m_spineSlidingFitWindow(20),
    m_maxXSeparation(5.f),
    m_maxSeparation(5.f),
    m_maxAngularDeviation(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoShowerMatchingTool::Run(const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
    const ProtoShowerVector &protoShowerVectorW, ProtoShowerMatchVector &protoShowerMatchVector)
{
    IntVector usedProtoShowersU, usedProtoShowersV, usedProtoShowersW; 

    for (unsigned int uIndex = 0; uIndex < protoShowerVectorU.size(); ++uIndex)
    {
        const ProtoShower &protoShowerU(protoShowerVectorU.at(uIndex));

        for (unsigned int vIndex = 0; vIndex < protoShowerVectorV.size(); ++vIndex)
        {
            const ProtoShower &protoShowerV(protoShowerVectorV.at(vIndex));

            for (unsigned int wIndex = 0; wIndex < protoShowerVectorW.size(); ++wIndex)
            {
                const ProtoShower &protoShowerW(protoShowerVectorW.at(wIndex));

                if (std::find(usedProtoShowersU.begin(), usedProtoShowersU.end(), uIndex) != usedProtoShowersU.end())
                    continue;

                if (std::find(usedProtoShowersV.begin(), usedProtoShowersV.end(), vIndex) != usedProtoShowersV.end())
                    continue;             

                if (std::find(usedProtoShowersW.begin(), usedProtoShowersW.end(), wIndex) != usedProtoShowersW.end())
                    continue;

                Consistency consistency(Consistency::POSITION);

                if (!this->ArePathwaysConsistent(protoShowerU, protoShowerV, protoShowerW, consistency))
                    continue;

                usedProtoShowersU.push_back(uIndex);
                usedProtoShowersV.push_back(vIndex);
                usedProtoShowersW.push_back(wIndex);

                protoShowerMatchVector.push_back(ProtoShowerMatch(protoShowerU, protoShowerV, protoShowerW, consistency));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::ArePathwaysConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
    Consistency &consistency) const
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

bool ProtoShowerMatchingTool::AreShowerStartsConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, 
    const ProtoShower &protoShowerW) const
{
    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

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

bool ProtoShowerMatchingTool::AreDirectionsConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, 
    const ProtoShower &protoShowerW) const
{
    const CartesianVector &directionU1(protoShowerU.m_connectionPathway.m_startDirection);
    const CartesianVector &directionV1(protoShowerV.m_connectionPathway.m_startDirection);
    const CartesianVector &directionW1(protoShowerW.m_connectionPathway.m_startDirection);

    if (this->AreDirectionsConsistent(directionU1, directionV1, directionW1))
    {
        return true;
    }
    else
    {
        const bool isDownstream(protoShowerW.m_showerCore.m_startPosition.GetZ() > protoShowerW.m_connectionPathway.m_startPosition.GetZ());

        CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

        for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
            spinePositionsU.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
            spinePositionsV.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
            spinePositionsW.push_back(pCaloHit->GetPositionVector());

        const TwoDSlidingFitResult spineFitU(&spinePositionsU, m_spineSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, m_spineSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, m_spineSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        const CartesianVector directionU2(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionV2(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionW2(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

        return this->AreDirectionsConsistent(directionU2, directionV2, directionW2);
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::AreDirectionsConsistent(CartesianVector directionU, CartesianVector directionV, CartesianVector directionW) const
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
    {
        int positiveCount(0);
        positiveCount += directionU.GetX() > 0.f ? 0 : 1;
        positiveCount += directionV.GetX() > 0.f ? 0 : 1;
        positiveCount += directionW.GetX() > 0.f ? 0 : 1;

        if (positiveCount >= 2)
        {
            directionU = CartesianVector(std::fabs(directionU.GetX()), 0.f, directionU.GetZ());
            directionV = CartesianVector(std::fabs(directionV.GetX()), 0.f, directionV.GetZ());
            directionW = CartesianVector(std::fabs(directionW.GetX()), 0.f, directionW.GetZ());
        }
    }

    if (directionU.GetX() * directionV.GetX() < 0.f)
        return false;

    if (directionU.GetX() * directionW.GetX() < 0.f)
        return false;

    if (directionV.GetX() * directionW.GetX() < 0.f)
        return false;
    
    const CartesianVector projectionU(LArGeometryHelper::MergeTwoDirections(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, directionV, directionW));
    const CartesianVector projectionV(LArGeometryHelper::MergeTwoDirections(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, directionW, directionU));
    const CartesianVector projectionW(LArGeometryHelper::MergeTwoDirections(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, directionU, directionV));

    float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180.f / M_PI);
    float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180.f / M_PI);
    float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180.f / M_PI);

    if (isIsochronous)
    {
        openingAngleU = std::min(openingAngleU, 180.f - openingAngleU);
        openingAngleV = std::min(openingAngleV, 180.f - openingAngleV);
        openingAngleW = std::min(openingAngleW, 180.f - openingAngleW);
    }

    if ((openingAngleU > m_maxAngularDeviation) || (openingAngleV > m_maxAngularDeviation) || (openingAngleW > m_maxAngularDeviation))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoShowerMatchingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SpineSlidingFitWindow", m_spineSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxSeparation", m_maxSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxAngularDeviation", m_maxAngularDeviation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
