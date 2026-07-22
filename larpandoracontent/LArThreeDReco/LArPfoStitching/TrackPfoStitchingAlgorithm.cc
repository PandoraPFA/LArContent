/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoStitching/TrackPfoStitchingAlgorithm.cc
 *
 *  @brief  Implementation of the track pfo stitching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "Geometry/DetectorGap.h"
#include "Geometry/LArTPC.h"

#include "larpandoracontent/LArThreeDReco/LArPfoStitching/TrackPfoStitchingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackPfoStitchingAlgorithm::InferredLArTPC::InferredLArTPC(const PandoraApi::Geometry::LArTPC::Parameters &parameters) :
    LArTPC(parameters)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackPfoStitchingAlgorithm::TrackPfoStitchingAlgorithm() :
    m_visualise{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackPfoStitchingAlgorithm::Run()
{
    PANDORA_RETURN_IF(STATUS_CODE_SUCCESS, m_stitchingToolVector.empty());

    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    InferredLArTPCVector inferredTPCs;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->InferLArTPCs(inferredTPCs));
    if (inferredTPCs.size() == 1) // Might be in a per-TPC worker
        return STATUS_CODE_SUCCESS;

    PfoToLArTPCMap pfoToLArTPCMap;
    this->GetPfoToLArTPCMap(pPfoList, inferredTPCs, pfoToLArTPCMap);
    if (pfoToLArTPCMap.empty())
        return STATUS_CODE_SUCCESS;

    PfoToFloatMap stitchedPfosToX0Map;
    for (StitchingBaseTool *const pStitchingTool : m_stitchingToolVector)
        pStitchingTool->Run(this, pPfoList, pfoToLArTPCMap, stitchedPfosToX0Map);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackPfoStitchingAlgorithm::InferLArTPCs(InferredLArTPCVector &inferredTPCs) const
{
    const LArTPC &parentLArTPC(this->GetPandora().GetGeometry()->GetLArTPC());
    const float parentMinX(parentLArTPC.GetCenterX() - 0.5f * parentLArTPC.GetWidthX());
    const float parentMaxX(parentLArTPC.GetCenterX() + 0.5f * parentLArTPC.GetWidthX());

    std::vector<std::pair<float, float>> driftGapXIntervals;
    for (const DetectorGap *const pDetectorGap : this->GetPandora().GetGeometry()->GetDetectorGapList())
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap *>(pDetectorGap));

        if (!pLineGap || (pLineGap->GetLineGapType() != TPC_DRIFT_GAP))
            continue;

        PANDORA_RETURN_IF(STATUS_CODE_INVALID_PARAMETER,
            pLineGap->GetLineStartX() > pLineGap->GetLineEndX() ||
            pLineGap->GetLineStartX() < parentMinX || pLineGap->GetLineEndX() > parentMaxX);
        driftGapXIntervals.emplace_back(pLineGap->GetLineStartX(), pLineGap->GetLineEndX());
    }

    std::sort(driftGapXIntervals.begin(), driftGapXIntervals.end(),
        [](const std::pair<float, float> &lhs, const std::pair<float, float> &rhs) { return lhs.first < rhs.first; });

    std::vector<std::pair<float, float>> inferredTPCXIntervals;
    float inferredTPCMinX{parentMinX};
    for (const std::pair<float, float> &interval : driftGapXIntervals)
    {
        inferredTPCXIntervals.emplace_back(inferredTPCMinX, interval.first);
        inferredTPCMinX = interval.second;
    }
    inferredTPCXIntervals.emplace_back(inferredTPCMinX, parentMaxX);

    unsigned int volumeId{0};
    bool isDriftInPositiveX{true};
    for (const std::pair<float, float> &interval : inferredTPCXIntervals)
    {
        PandoraApi::Geometry::LArTPC::Parameters parameters;
        parameters.m_larTPCVolumeId = volumeId++;
        parameters.m_centerX = 0.5f * (interval.first + interval.second);
        parameters.m_centerY = parentLArTPC.GetCenterY();
        parameters.m_centerZ = parentLArTPC.GetCenterZ();
        parameters.m_widthX = interval.second - interval.first;
        parameters.m_widthY = parentLArTPC.GetWidthY();
        parameters.m_widthZ = parentLArTPC.GetWidthZ();
        parameters.m_wirePitchU = parentLArTPC.GetWirePitchU();
        parameters.m_wirePitchV = parentLArTPC.GetWirePitchV();
        parameters.m_wirePitchW = parentLArTPC.GetWirePitchW();
        parameters.m_wireAngleU = parentLArTPC.GetWireAngleU();
        parameters.m_wireAngleV = parentLArTPC.GetWireAngleV();
        parameters.m_wireAngleW = parentLArTPC.GetWireAngleW();
        parameters.m_sigmaUVW = parentLArTPC.GetSigmaUVW();
        parameters.m_isDriftInPositiveX = isDriftInPositiveX;
        isDriftInPositiveX = !isDriftInPositiveX;

        inferredTPCs.push_back(std::unique_ptr<InferredLArTPC>(new InferredLArTPC(parameters)));
    }

    PANDORA_RETURN_IF(STATUS_CODE_FAILURE, inferredTPCs.empty());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackPfoStitchingAlgorithm::GetPfoToLArTPCMap(const PfoList *const pPfoList, const InferredLArTPCVector &inferredTPCs,
    PfoToLArTPCMap &pfoToLArTPCMap) const
{
    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (!LArPfoHelper::IsTrack(pPfo) || LArPfoHelper::IsNeutrino(pPfo))
            continue;

        CartesianPointVector positions;
        LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, positions);

        const LArTPC *pLArTPC{nullptr};
        for (const CartesianVector &pos : positions)
        {
            const LArTPC *pPosLArTPC{nullptr};
            for (const std::unique_ptr<InferredLArTPC> &pInferredLArTPC : inferredTPCs)
            {
                const float minX{pInferredLArTPC->GetCenterX() - 0.5f * pInferredLArTPC->GetWidthX()};
                const float maxX{pInferredLArTPC->GetCenterX() + 0.5f * pInferredLArTPC->GetWidthX()};

                if (pos.GetX() >= minX && pos.GetX() <= maxX)
                {
                    pPosLArTPC = pInferredLArTPC.get();
                    break;
                }
            }
            if (!pPosLArTPC) // 3D hit ended up in a gap or outside detector, ignoring it
                continue;

            if (pLArTPC && pLArTPC != pPosLArTPC)
            {
                pLArTPC = nullptr;
                break;
            }

            pLArTPC = pPosLArTPC;
        }

        if (pLArTPC) // The Pfo is within a single LArTPC
            pfoToLArTPCMap.emplace(pPfo, pLArTPC);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string TrackPfoStitchingAlgorithm::GetListName(const Cluster *const pCluster) const
{
    for (const std::string &listName : m_daughterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        if ((PandoraContentApi::GetList(*this, listName, pClusterList) == STATUS_CODE_SUCCESS) &&
            (pClusterList->end() != std::find(pClusterList->begin(), pClusterList->end(), pCluster)))
        {
            return listName;
        }
    }

    PANDORA_THROW(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string TrackPfoStitchingAlgorithm::GetListName(const Vertex *const pVertex) const
{
    for (const std::string &listName : m_daughterListNames)
    {
        const VertexList *pVertexList{nullptr};
        if ((PandoraContentApi::GetList(*this, listName, pVertexList) == STATUS_CODE_SUCCESS) &&
            (pVertexList->end() != std::find(pVertexList->begin(), pVertexList->end(), pVertex)))
        {
            return listName;
        }
    }

    PANDORA_THROW(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackPfoStitchingAlgorithm::ShiftPfoHierarchy([[maybe_unused]]const ParticleFlowObject *const pParentPfo,
    [[maybe_unused]]const PfoToLArTPCMap &pfoToLArTPCMap, [[maybe_unused]]const float x0) const
{
    std::cout << "TrackPfoStitchingAlgorithm: Stitching tools must not be configured to apply an x0 correction" << std::endl;
    PANDORA_THROW(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackPfoStitchingAlgorithm::StitchPfos(
    const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete, PfoToLArTPCMap &pfoToLArTPCMap) const
{
    PANDORA_THROW_IF(STATUS_CODE_NOT_ALLOWED, pPfoToEnlarge == pPfoToDelete);

    if (m_visualise)
    {
        PfoList pfoToEnlargeList{pPfoToEnlarge};
        PfoList pfoToDeleteList;
        LArPfoHelper::GetAllDownstreamPfos(pPfoToDelete, pfoToDeleteList);
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoToEnlargeList, "PfoToEnlarge", GREEN, true, false));
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoToDeleteList, "PfoToDelete", RED, true, false));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    pfoToLArTPCMap.erase(pPfoToEnlarge);
    pfoToLArTPCMap.erase(pPfoToDelete);

    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, m_inputPfoListName));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoToDelete, pDaughterPfo));
    }

    for (const Vertex *const pDaughterVertex : daughterVertices)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::Delete(*this, pDaughterVertex, this->GetListName(pDaughterVertex)));
    }

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *const pParentCluster(PfoMopUpBaseAlgorithm::GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::MergeAndDeleteClusters(
                    *this, pParentCluster, pDaughterCluster, this->GetListName(pParentCluster), this->GetListName(pDaughterCluster)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackPfoStitchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "DaughterListNames", m_daughterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "StitchingTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
    {
        StitchingBaseTool *const pStitchingTool(dynamic_cast<StitchingBaseTool *>(pAlgorithmTool));

        PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, !pStitchingTool);

        m_stitchingToolVector.push_back(pStitchingTool);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
