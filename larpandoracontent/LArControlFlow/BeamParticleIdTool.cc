/**
 *  @file   larpandoracontent/LArControlFlow/BeamParticleIdTool.cc
 *
 *  @brief  Implementation of the beam particle id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/BeamParticleIdTool.h"

#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

BeamParticleIdTool::BeamParticleIdTool() :
    m_selectAllBeamParticles(false),
    m_selectOnlyFirstSliceBeamParticles(false),
    m_visualizeID(false),
    m_projectionIntersectionCut(100.f),
    m_beamTPCIntersection(CartesianVector(0.f,0.f,0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BeamParticleIdTool::Initialize()
{
    // Get global TPC geometry information
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);

    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    m_tpcMinX = parentMinX;
    m_tpcMaxX = parentMaxX;
    m_tpcMinY = parentMinY;
    m_tpcMaxY = parentMaxY;
    m_tpcMinZ = parentMinZ;
    m_tpcMaxZ = parentMaxZ;

    CartesianVector normalTop(0,0,1), pointTop(0,0,m_tpcMaxZ);
    const Plane *pPlaneTop = new Plane(normalTop, pointTop);
    m_tpcPlanes.push_back(pPlaneTop);
  
    CartesianVector normalBottom(0,0,-1), pointBottom(0,0,m_tpcMinZ);
    const Plane *pPlaneBottom = new Plane(normalBottom, pointBottom);
    m_tpcPlanes.push_back(pPlaneBottom);
  
    CartesianVector normalRight(1,0,0), pointRight(m_tpcMaxX,0,0);
    const Plane *pPlaneRight = new Plane(normalRight, pointRight);
    m_tpcPlanes.push_back(pPlaneRight);
  
    CartesianVector normalLeft(-1,0,0), pointLeft(m_tpcMinX,0,0);
    const Plane *pPlaneLeft = new Plane(normalLeft, pointLeft);
    m_tpcPlanes.push_back(pPlaneLeft);
  
    CartesianVector normalBack(0,1,0), pointBack(0,m_tpcMaxY,0);
    const Plane *pPlaneBack = new Plane(normalBack, pointBack);
    m_tpcPlanes.push_back(pPlaneBack);
  
    CartesianVector normalFront(0,-1,0), pointFront(0,m_tpcMinY,0);
    const Plane *pPlaneFront = new Plane(normalFront, pointFront);
    m_tpcPlanes.push_back(pPlaneFront);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamParticleIdTool::SelectOutputPfos(const SliceHypotheses &beamSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (beamSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (m_selectAllBeamParticles || m_selectOnlyFirstSliceBeamParticles)
    {
        // Deafult Logic
        for (unsigned int sliceIndex = 0, nSlices = beamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
        {
            const PfoList &sliceOutput((m_selectAllBeamParticles || (m_selectOnlyFirstSliceBeamParticles && (0 == sliceIndex))) ?
                beamSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));

            selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
        }
    }
    else
    {
        for (unsigned int sliceIndex = 0, nSlices = beamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
        {
            const PfoList pfoListBeam(beamSliceHypotheses.at(sliceIndex));

            PfoList allConnectedPfoList;
            LArPfoHelper::GetAllConnectedPfos(pfoListBeam, allConnectedPfoList);

            CaloHitList caloHitList3D;
            LArPfoHelper::GetCaloHits(allConnectedPfoList, TPC_3D, caloHitList3D);

            if (caloHitList3D.empty())
                continue;

            CartesianVector centroid(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVecs;
            LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
            LArPcaHelper::RunPca(caloHitList3D, centroid, eigenValues, eigenVecs);
            CartesianVector majorAxis(eigenVecs.at(0));
            CartesianVector interceptOne(0.f,0.f,0.f), interceptTwo(0.f,0.f,0.f);

            try
            {
                this->GetTPCIntercepts(centroid, majorAxis, interceptOne, interceptTwo);
            }
            catch (...)
            {
                 std::cout << "Unable to find intersection of major axis of slice with TPC" << std::endl;
            }

            // Find the intercept closest to the beam TPC intersection
            float separationOne(std::sqrt(interceptOne.GetDistanceSquared(m_beamTPCIntersection))), separationTwo(std::sqrt(interceptTwo.GetDistanceSquared(m_beamTPCIntersection)));

            if (std::min(separationOne, separationTwo) < m_projectionIntersectionCut)
            {
                selectedPfos.insert(selectedPfos.end(), pfoListBeam.begin(), pfoListBeam.end());
            }

#ifdef MONITORING
            if (m_visualizeID)
            {
                CaloHitList beamCaloHitList2du, cosmicCaloHitList2du, beamCaloHitList2dv, cosmicCaloHitList2dv, beamCaloHitList2dw, cosmicCaloHitList2dw;
                std::string name(std::to_string(sliceIndex));
                PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &m_beamTPCIntersection, "BeamTPCIntersection", RED, 1);
                PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &interceptOne, &interceptTwo, "Major_Axis_Slice_" + name, BLACK, 1, 2);
                PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoListBeam, "Beam_Slice_" + name, VIOLET);
                PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &beamCaloHitList2du, "Beam_CaloHits_2DU_" + name, BLUE);
                PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &cosmicCaloHitList2du, "Cosmic_CaloHits_2DU_" + name, RED);
                PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &beamCaloHitList2dv, "Beam_CaloHits_2DV_" + name, BLUE);
                PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &cosmicCaloHitList2dv, "Cosmic_CaloHits_2DV_" + name, RED);
                PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &beamCaloHitList2dw, "Beam_CaloHits_2DW_" + name, BLUE);
                PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &cosmicCaloHitList2dw, "Cosmic_CaloHits_2DW_" + name, RED);
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }
#endif
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode BeamParticleIdTool::GetTPCIntercepts(CartesianVector &a0, CartesianVector &lineDirection, CartesianVector &interceptOne,  CartesianVector &interceptTwo) const
{
    CartesianVector a(lineDirection.GetUnitVector());
    std::vector<CartesianVector> intercepts;

    for (const Plane *const pPlane : m_tpcPlanes)
    {
        CartesianVector intercept(pPlane->GetLineIntersection(a0, a));
        if (this->IsContained(intercept))
        {
            intercepts.push_back(intercept);
        }
    }

    if (intercepts.size() == 2)
    {
        interceptOne = intercepts.at(0);
        interceptTwo = intercepts.at(1);
    }
    else 
    {
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BeamParticleIdTool::IsContained(CartesianVector &spacePoint) const 
{
    // TODO: Make cleaner floating point comparisons 
    float safetyFactor(0.0001);
    bool isContainedX(m_tpcMinX - safetyFactor < spacePoint.GetX() && spacePoint.GetX() < m_tpcMaxX + safetyFactor);
    bool isContainedY(m_tpcMinY - safetyFactor < spacePoint.GetY() && spacePoint.GetY() < m_tpcMaxY + safetyFactor);
    bool isContainedZ(m_tpcMinZ - safetyFactor < spacePoint.GetZ() && spacePoint.GetZ() < m_tpcMaxZ + safetyFactor);

    return isContainedX && isContainedY && isContainedZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

BeamParticleIdTool::Plane::Plane(CartesianVector &normal, CartesianVector &point) : 
    m_unitNormal(0.f, 0.f, 0.f),
    m_point(0.f, 0.f, 0.f)
{
    m_point = point;
    m_unitNormal = normal.GetUnitVector();
    m_a = m_unitNormal.GetX();
    m_b = m_unitNormal.GetY();
    m_c = m_unitNormal.GetZ();
    m_d = -1.f * (normal.GetDotProduct(point));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector BeamParticleIdTool::Plane::GetLineIntersection(CartesianVector &a0, CartesianVector &a) const
{
    if (std::fabs(a.GetDotProduct(m_unitNormal)) < std::numeric_limits<float>::min())
        return CartesianVector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    float t(-1.f * (a0.GetDotProduct(m_unitNormal) + m_d) / (a.GetDotProduct(m_unitNormal)));
    return CartesianVector(a.GetX() * t + a0.GetX(), a.GetY() * t + a0.GetY(), a.GetZ() * t + a0.GetZ());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BeamParticleIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectAllBeamParticles", m_selectAllBeamParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectOnlyFirstSliceBeamParticles", m_selectOnlyFirstSliceBeamParticles));

    // TODO Check logic - maybe remove these options completely?
    if ((m_selectAllBeamParticles || m_selectOnlyFirstSliceBeamParticles) && (m_selectAllBeamParticles == m_selectOnlyFirstSliceBeamParticles))
    {
        std::cout << "BeamParticleIdTool::ReadSettings - exactly one of SelectAllBeamParticles and SelectOnlyFirstSliceBeamParticles must be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ProjectionIntersectionCut", m_projectionIntersectionCut));

    FloatVector beamTPCIntersection;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "BeamTPCIntersection", beamTPCIntersection));
 
    if (3 == beamTPCIntersection.size())
    {
        m_beamTPCIntersection.SetValues(beamTPCIntersection.at(0), beamTPCIntersection.at(1), beamTPCIntersection.at(2));
    }
    else
    {
        // Default for protoDUNE.
        m_beamTPCIntersection.SetValues(-33.051, 461.06, 0);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeID", m_visualizeID));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
