/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPBaseHierarchyTool.cc
 *
 *  @brief  Implementation of the MLP base hierarchy tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPBaseHierarchyTool.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

MLPBaseHierarchyTool::MLPBaseHierarchyTool() :
    m_detectorMinX(std::numeric_limits<float>::max()),
    m_detectorMaxX(-std::numeric_limits<float>::max()),
    m_detectorMinY(std::numeric_limits<float>::max()),
    m_detectorMaxY(-std::numeric_limits<float>::max()),
    m_detectorMinZ(std::numeric_limits<float>::max()),
    m_detectorMaxZ(-std::numeric_limits<float>::max()),
    areBoundariesSet(false),
    m_vertexRegionRadius(5.f),
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D"})
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPBaseHierarchyTool::SetDetectorBoundaries()
{
    if (areBoundariesSet)
        return;

    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());

    for (const auto &entry : larTPCMap)
    {
        m_detectorMinX = std::min(m_detectorMinX, entry.second->GetCenterX() - (entry.second->GetWidthX() * 0.5f));
        m_detectorMaxX = std::max(m_detectorMaxX, entry.second->GetCenterX() + (entry.second->GetWidthX() * 0.5f));
        m_detectorMinY = std::min(m_detectorMinY, entry.second->GetCenterY() - (entry.second->GetWidthY() * 0.5f));
        m_detectorMaxY = std::max(m_detectorMaxY, entry.second->GetCenterY() + (entry.second->GetWidthY() * 0.5f));
        m_detectorMinZ = std::min(m_detectorMinZ, entry.second->GetCenterZ() - (entry.second->GetWidthZ() * 0.5f));
        m_detectorMaxZ = std::max(m_detectorMaxZ, entry.second->GetCenterZ() + (entry.second->GetWidthZ() * 0.5f));
    }

    areBoundariesSet = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPBaseHierarchyTool::IsInFV(const CartesianVector &position) const
{
    if ((position.GetX() < m_detectorMinX) or (position.GetX() > m_detectorMaxX))
        return false;

    if ((position.GetY() < m_detectorMinY) or (position.GetY() > m_detectorMaxY))
        return false;

    if ((position.GetZ() < m_detectorMinZ) or (position.GetZ() > m_detectorMaxZ))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPBaseHierarchyTool::GetNSpacepoints(const HierarchyPfo &hierarchyPfo)
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(hierarchyPfo.GetPfo(), clusterList3D);

    int total3DHits(0);

    for (const Cluster *const pCluster3D : clusterList3D)
        total3DHits += pCluster3D->GetNCaloHits();

    return total3DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<float, float> MLPBaseHierarchyTool::GetParticleInfoAboutPfoPosition(const Algorithm *const pAlgorithm, 
    const ParticleFlowObject *const pPfo, const CartesianVector &pointOfInterest) const
{ 
    int hitCount = 0;
    int particleCount = 0;    

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
                const double sepSq = (otherPfoPosition - pointOfInterest).GetMagnitudeSquared();

                if (sepSq < (m_vertexRegionRadius * m_vertexRegionRadius))
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

void MLPBaseHierarchyTool::NormaliseNetworkParam(const float minLimit, const float maxLimit, float &networkParam) const
{
    const float interval(std::fabs(minLimit) + std::fabs(maxLimit));

    if (networkParam < minLimit)
        networkParam = minLimit;

    if (networkParam > maxLimit)
        networkParam = maxLimit;

    networkParam /= interval;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPBaseHierarchyTool::IsVectorSet(const CartesianVector &vector)
{
    return (vector.GetZ() > -990.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int MLPBaseHierarchyTool::AddToInput(const int startIndex, const FloatVector &paramVector, LArDLHelper::TorchInput &modelInput)
{
    int insertIndex(startIndex);

    for (float param : paramVector)
    {
        modelInput[0][insertIndex] = param;
        ++insertIndex;
    }

    return insertIndex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPBaseHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionRadius", m_vertexRegionRadius));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
