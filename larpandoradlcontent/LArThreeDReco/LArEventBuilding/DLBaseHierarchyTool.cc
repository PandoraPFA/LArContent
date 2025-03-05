/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLBaseHierarchyTool.cc
 *
 *  @brief  Implementation of the DL base hierarchy tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLBaseHierarchyTool.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLBaseHierarchyTool::DLBaseHierarchyTool() :
    m_vertexRegionRadiusSq(25.f),
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D"}),
    m_areBoundariesSet(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLBaseHierarchyTool::SetDetectorBoundaries()
{
    if (m_areBoundariesSet)
        return;

    m_detectorBoundaries = LArGeometryHelper::GetDetectorBoundaries(this->GetPandora());
    m_areBoundariesSet = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<float, float> DLBaseHierarchyTool::GetParticleInfoAboutPfoPosition(
    const Algorithm *const pAlgorithm, const ParticleFlowObject *const pPfo, const CartesianVector &pointOfInterest) const
{
    int hitCount(0);
    int particleCount(0);

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        if (PandoraContentApi::GetList(*pAlgorithm, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject *const pOtherPfo : *pPfoList)
        {
            if (pPfo == pOtherPfo)
                continue;

            bool isClose(false);

            CartesianPointVector otherPfoPositions3D;
            LArPfoHelper::GetCoordinateVector(pOtherPfo, TPC_3D, otherPfoPositions3D);

            for (const CartesianVector &otherPfoPosition : otherPfoPositions3D)
            {
                const double sepSq((otherPfoPosition - pointOfInterest).GetMagnitudeSquared());

                if (sepSq < m_vertexRegionRadiusSq)
                {
                    isClose = true;
                    ++hitCount;
                }
            }

            if (isClose)
                ++particleCount;
        }
    }

    return std::pair<float, float>({hitCount, particleCount});
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLBaseHierarchyTool::NormaliseNetworkParam(const float minLimit, const float maxLimit, float &networkParam) const
{
    const float interval(std::fabs(maxLimit - minLimit));

    if (interval < std::numeric_limits<float>::epsilon())
    {
        networkParam = minLimit;
        return;
    }

    if (networkParam < minLimit)
        networkParam = minLimit;

    if (networkParam > maxLimit)
        networkParam = maxLimit;

    networkParam /= interval;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLBaseHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    StatusCode statusCode(XmlHelper::ReadValue(xmlHandle, "VertexRegionRadius", m_vertexRegionRadiusSq));

    if (statusCode == STATUS_CODE_SUCCESS)
    {
        m_vertexRegionRadiusSq = m_vertexRegionRadiusSq * m_vertexRegionRadiusSq;
    }
    else if (statusCode != STATUS_CODE_NOT_FOUND)
    {
        return statusCode;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
