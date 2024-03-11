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

using namespace pandora;

namespace lar_content
{

BeamParticleIdTool::BeamParticleIdTool() :
    m_selectAllBeamParticles(false),
    m_selectOnlyFirstSliceBeamParticles(false),
    m_tpcMinX(std::numeric_limits<float>::max()),
    m_tpcMaxX(-std::numeric_limits<float>::max()),
    m_tpcMinY(std::numeric_limits<float>::max()),
    m_tpcMaxY(-std::numeric_limits<float>::max()),
    m_tpcMinZ(std::numeric_limits<float>::max()),
    m_tpcMaxZ(-std::numeric_limits<float>::max()),
    m_beamTPCIntersection(0.f, 0.f, 0.f),
    m_beamDirection(0.f, 0.f, 0.f),
    m_projectionIntersectionCut(100.f),
    m_closestDistanceCut(100.f),
    m_angleToBeamCut(150.f * M_PI / 180.f),
    m_selectedFraction(10.f),
    m_nSelectedHits(100)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamParticleIdTool::SelectOutputPfos(const Algorithm *const pAlgorithm, const SliceHypotheses &beamSliceHypotheses,
    const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (beamSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // First, simple approach
    if (m_selectAllBeamParticles || m_selectOnlyFirstSliceBeamParticles)
    {
        for (unsigned int sliceIndex = 0, nSlices = beamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
        {
            const PfoList &sliceOutput((m_selectAllBeamParticles || (m_selectOnlyFirstSliceBeamParticles && (0 == sliceIndex)))
                    ? beamSliceHypotheses.at(sliceIndex)
                    : crSliceHypotheses.at(sliceIndex));

            const float score(m_selectAllBeamParticles || (m_selectOnlyFirstSliceBeamParticles && (0 == sliceIndex)) ? 1.f : -1.f);

            for (const ParticleFlowObject *const pPfo : sliceOutput)
            {
                object_creation::ParticleFlowObject::Metadata metadata;
                metadata.m_propertiesToAdd["TestBeamScore"] = score;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
            }

            selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
        }

        return;
    }

    // Now start to examine topology of beam slice hypotheses
    for (unsigned int sliceIndex = 0, nSlices = beamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        bool usebeamHypothesis(false);

        try
        {
            PfoList allConnectedPfoList;
            LArPfoHelper::GetAllConnectedPfos(beamSliceHypotheses.at(sliceIndex), allConnectedPfoList);

            CaloHitList caloHitList3D;
            LArPfoHelper::GetCaloHits(allConnectedPfoList, TPC_3D, caloHitList3D);

            CaloHitList selectedCaloHitList;
            float closestDistance(std::numeric_limits<float>::max());
            this->GetSelectedCaloHits(caloHitList3D, selectedCaloHitList, closestDistance);

            if (!selectedCaloHitList.empty())
            {
                CartesianVector centroidSel(0.f, 0.f, 0.f);
                LArPcaHelper::EigenVectors eigenVecsSel;
                LArPcaHelper::EigenValues eigenValuesSel(0.f, 0.f, 0.f);
                LArPcaHelper::RunPca(selectedCaloHitList, centroidSel, eigenValuesSel, eigenVecsSel);

                const CartesianVector &majorAxisSel(eigenVecsSel.front());
                const float supplementaryAngleToBeam(majorAxisSel.GetOpeningAngle(m_beamDirection));

                CartesianVector interceptOne(0.f, 0.f, 0.f), interceptTwo(0.f, 0.f, 0.f);
                this->GetTPCIntercepts(centroidSel, majorAxisSel, interceptOne, interceptTwo);

                const float separationOne((interceptOne - m_beamTPCIntersection).GetMagnitude());
                const float separationTwo((interceptTwo - m_beamTPCIntersection).GetMagnitude());

                if ((std::min(separationOne, separationTwo) < m_projectionIntersectionCut) && (closestDistance < m_closestDistanceCut) &&
                    (supplementaryAngleToBeam > m_angleToBeamCut))
                {
                    usebeamHypothesis = true;
                }
            }
        }
        catch (const StatusCodeException &)
        {
            usebeamHypothesis = false;
        }

        const PfoList &sliceOutput(usebeamHypothesis ? beamSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));
        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());

        const float score(usebeamHypothesis ? 1.f : -1.f);

        for (const ParticleFlowObject *const pPfo : sliceOutput)
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["TestBeamScore"] = score;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }
    }
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

    const CartesianVector normalTop(0.f, 0.f, 1.f), pointTop(0.f, 0.f, m_tpcMaxZ);
    const CartesianVector normalBottom(0.f, 0.f, -1.f), pointBottom(0.f, 0.f, m_tpcMinZ);
    const CartesianVector normalRight(1.f, 0.f, 0.f), pointRight(m_tpcMaxX, 0.f, 0.f);
    const CartesianVector normalLeft(-1.f, 0.f, 0.f), pointLeft(m_tpcMinX, 0.f, 0.f);
    const CartesianVector normalBack(0.f, 1.f, 0.f), pointBack(0.f, m_tpcMaxY, 0.f);
    const CartesianVector normalFront(0.f, -1.f, 0.f), pointFront(0.f, m_tpcMinY, 0.f);
    m_tpcPlanes.emplace_back(normalTop, pointTop);
    m_tpcPlanes.emplace_back(normalBottom, pointBottom);
    m_tpcPlanes.emplace_back(normalRight, pointRight);
    m_tpcPlanes.emplace_back(normalLeft, pointLeft);
    m_tpcPlanes.emplace_back(normalBack, pointBack);
    m_tpcPlanes.emplace_back(normalFront, pointFront);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamParticleIdTool::GetSelectedCaloHits(const CaloHitList &inputCaloHitList, CaloHitList &outputCaloHitList, float &closestHitToFaceDistance) const
{
    if (inputCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    typedef std::pair<const CaloHit *, float> HitDistancePair;
    typedef std::vector<HitDistancePair> HitDistanceVector;
    HitDistanceVector hitDistanceVector;

    for (const CaloHit *const pCaloHit : inputCaloHitList)
        hitDistanceVector.emplace_back(pCaloHit, (pCaloHit->GetPositionVector() - m_beamTPCIntersection).GetMagnitudeSquared());

    std::sort(hitDistanceVector.begin(), hitDistanceVector.end(),
        [](const HitDistancePair &lhs, const HitDistancePair &rhs) -> bool { return (lhs.second < rhs.second); });
    closestHitToFaceDistance = std::sqrt(hitDistanceVector.front().second);

    const unsigned int nInputHits(inputCaloHitList.size());
    const unsigned int nSelectedCaloHits(nInputHits < m_nSelectedHits
            ? nInputHits
            : static_cast<unsigned int>(std::round(static_cast<float>(nInputHits) * m_selectedFraction / 100.f + 0.5f)));

    for (const HitDistancePair &hitDistancePair : hitDistanceVector)
    {
        outputCaloHitList.push_back(hitDistancePair.first);

        if (outputCaloHitList.size() >= nSelectedCaloHits)
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamParticleIdTool::GetTPCIntercepts(
    const CartesianVector &a0, const CartesianVector &lineDirection, CartesianVector &interceptOne, CartesianVector &interceptTwo) const
{
    CartesianPointVector intercepts;
    const CartesianVector lineUnitVector(lineDirection.GetUnitVector());

    for (const Plane &plane : m_tpcPlanes)
    {
        const CartesianVector intercept(plane.GetLineIntersection(a0, lineUnitVector));

        if (this->IsContained(intercept))
            intercepts.push_back(intercept);
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BeamParticleIdTool::IsContained(const CartesianVector &spacePoint) const
{
    if ((m_tpcMinX - spacePoint.GetX() > std::numeric_limits<float>::epsilon()) ||
        (spacePoint.GetX() - m_tpcMaxX > std::numeric_limits<float>::epsilon()) ||
        (m_tpcMinY - spacePoint.GetY() > std::numeric_limits<float>::epsilon()) ||
        (spacePoint.GetY() - m_tpcMaxY > std::numeric_limits<float>::epsilon()) ||
        (m_tpcMinZ - spacePoint.GetZ() > std::numeric_limits<float>::epsilon()) ||
        (spacePoint.GetZ() - m_tpcMaxZ > std::numeric_limits<float>::epsilon()))
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

BeamParticleIdTool::Plane::Plane(const CartesianVector &normal, const CartesianVector &point) :
    m_unitNormal(normal.GetUnitVector()),
    m_point(point),
    m_d(-1.f * (normal.GetDotProduct(point)))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector BeamParticleIdTool::Plane::GetLineIntersection(const CartesianVector &a0, const CartesianVector &a) const
{
    if (std::fabs(a.GetDotProduct(m_unitNormal)) < std::numeric_limits<float>::min())
        return CartesianVector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    const float denominator(a.GetDotProduct(m_unitNormal));

    if (std::fabs(denominator) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    const float t(-1.f * (a0.GetDotProduct(m_unitNormal) + m_d) / denominator);
    return (a0 + (a * t));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BeamParticleIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectAllBeamParticles", m_selectAllBeamParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SelectOnlyFirstSliceBeamParticles", m_selectOnlyFirstSliceBeamParticles));

    if (m_selectAllBeamParticles && m_selectOnlyFirstSliceBeamParticles)
    {
        std::cout << "BeamParticleIdTool::ReadSettings - cannot use both SelectAllBeamParticles and SelectOnlyFirstSliceBeamParticles simultaneously"
                  << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    FloatVector beamTPCIntersection;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "BeamTPCIntersection", beamTPCIntersection));

    if (3 == beamTPCIntersection.size())
    {
        m_beamTPCIntersection.SetValues(beamTPCIntersection.at(0), beamTPCIntersection.at(1), beamTPCIntersection.at(2));
    }
    else if (!beamTPCIntersection.empty())
    {
        std::cout << "BeamParticleIdTool::ReadSettings - invalid BeamTPCIntersection specified " << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    else
    {
        // Default for protoDUNE.
        m_beamTPCIntersection.SetValues(-33.051, 461.06, 0);
    }

    FloatVector beamDirection;
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "BeamDirection", beamDirection));

    if (3 == beamDirection.size())
    {
        m_beamDirection.SetValues(beamDirection.at(0), beamDirection.at(1), beamDirection.at(2));
    }
    else if (!beamDirection.empty())
    {
        std::cout << "BeamParticleIdTool::ReadSettings - invalid BeamDirection specified " << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    else
    {
        // Default for protoDUNE.
        const float thetaXZ0(-11.844f * M_PI / 180.f);
        m_beamDirection.SetValues(std::sin(thetaXZ0), 0, std::cos(thetaXZ0));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ProjectionIntersectionCut", m_projectionIntersectionCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClosestDistanceCut", m_closestDistanceCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AngleToBeamCut", m_angleToBeamCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectedFraction", m_selectedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSelectedHits", m_nSelectedHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
