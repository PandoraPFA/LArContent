/**
 *  @file   LArContent/src/LArStitching/StitchingObjectCreationTool.cc
 * 
 *  @brief  Implementation of the stitching object creation tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArStitching/StitchingObjectCreationTool.h"

using namespace pandora;
using namespace lar_content;

namespace lar_content
{

void StitchingObjectCreationTool::Run(const StitchingAlgorithm *const pAlgorithm, StitchingInfo &stitchingInfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    std::string clusterListName;
    const ClusterList *pClusterList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*pAlgorithm, pClusterList, clusterListName));

    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*pAlgorithm, pVertexList, vertexListName));

    std::string pfoListName;
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*pAlgorithm, pPfoList, pfoListName));

    const PandoraInstanceList &pandoraInstances(MultiPandoraApi::GetDaughterPandoraInstanceList(&(this->GetPandora())));

    for (const Pandora *const pPandora : pandoraInstances)
        this->Recreate3DContent(pAlgorithm, pPandora, MultiPandoraApi::GetVolumeInfo(pPandora), stitchingInfo);

    if (!pClusterList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*pAlgorithm, m_newClusterListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, m_newClusterListName));
    }

    if (!pVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*pAlgorithm, m_newVertexListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*pAlgorithm, m_newVertexListName));
    }

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*pAlgorithm, m_newPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*pAlgorithm, m_newPfoListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingObjectCreationTool::Recreate3DContent(const Algorithm *const pAlgorithm, const Pandora *const pPandora, const VolumeInfo &volumeInfo,
    StitchingInfo &stitchingInfo) const
{
    const PfoList *pPfoList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pPandora, pPfoList));

    for (const ParticleFlowObject *const pInputPfo : *pPfoList)
    {
        if (!pInputPfo->GetParentPfoList().empty())
            continue;

        this->Recreate3DContent(pAlgorithm, pInputPfo, nullptr, pPandora, volumeInfo, stitchingInfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingObjectCreationTool::Recreate3DContent(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pInputPfo,
    const ParticleFlowObject *const pNewParentPfo, const Pandora *const pPandora, const VolumeInfo &volumeInfo, StitchingInfo &stitchingInfo) const
{
    ClusterList inputClusterList, newClusterList;
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList);

    for (const Cluster *const pInputCluster : inputClusterList)
    {
        CaloHitList inputCaloHitList, newCaloHitList;
        pInputCluster->GetOrderedCaloHitList().GetCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.insert(this->CreateCaloHit(pAlgorithm, pInputCaloHit, pInputPfo, volumeInfo));

        if (!newCaloHitList.empty())
            newClusterList.insert(this->CreateCluster(pAlgorithm, pInputCluster, newCaloHitList));
    }

    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
    {
        if (VERTEX_3D == pInputVertex->GetVertexType())
            newVertexList.insert(this->CreateVertex(pAlgorithm, pInputVertex, pInputPfo, volumeInfo));
    }

    const ParticleFlowObject *const pNewPfo = this->CreatePfo(pAlgorithm, pInputPfo, newClusterList, newVertexList);
    this->AddStitchingInfo(pNewPfo, pPandora, stitchingInfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*pAlgorithm, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate3DContent(pAlgorithm, pInputDaughterPfo, pNewPfo, pPandora, volumeInfo, stitchingInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *StitchingObjectCreationTool::CreateCaloHit(const Algorithm *const pAlgorithm, const CaloHit *const pInputCaloHit,
    const ParticleFlowObject *const pInputPfo, const VolumeInfo &volumeInfo) const
{
    if (TPC_3D != pInputCaloHit->GetHitType())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float x0(volumeInfo.GetParticleX0(pInputPfo));
    const CartesianVector xOffset(volumeInfo.IsDriftInPositiveX() ? -x0 : x0, 0.f, 0.f);

    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pInputCaloHit->GetPositionVector() + volumeInfo.GetCenter() + xOffset;
    parameters.m_expectedDirection = pInputCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pInputCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pInputCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pInputCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pInputCaloHit->GetCellSize1();
    parameters.m_cellThickness = pInputCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pInputCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pInputCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pInputCaloHit->GetTime();
    parameters.m_inputEnergy = pInputCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pInputCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pInputCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pInputCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pInputCaloHit->IsDigital();
    parameters.m_hitType = pInputCaloHit->GetHitType();
    parameters.m_hitRegion = pInputCaloHit->GetHitRegion();
    parameters.m_layer = pInputCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pInputCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = pInputCaloHit->GetParentCaloHitAddress();

    const CaloHit *pNewCaloHit(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*pAlgorithm, parameters, pNewCaloHit));

    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = pInputCaloHit->IsIsolated();
    metadata.m_isPossibleMip = pInputCaloHit->IsPossibleMip();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AlterMetadata(*pAlgorithm, pNewCaloHit, metadata));

    return pNewCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *StitchingObjectCreationTool::CreateCluster(const Algorithm *const pAlgorithm, const Cluster *const pInputCluster,
    const CaloHitList &newCaloHitList) const
{
    if (TPC_3D != LArClusterHelper::GetClusterHitType(pInputCluster))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = newCaloHitList;

    const Cluster *pNewCluster(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, parameters, pNewCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = pInputCluster->GetParticleIdFlag();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AlterMetadata(*pAlgorithm, pNewCluster, metadata));

    return pNewCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *StitchingObjectCreationTool::CreateVertex(const Algorithm *const pAlgorithm, const Vertex *const pInputVertex,
    const ParticleFlowObject *const pInputPfo, const VolumeInfo &volumeInfo) const
{
    if (VERTEX_3D != pInputVertex->GetVertexType())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float x0(volumeInfo.GetParticleX0(pInputPfo));
    const CartesianVector xOffset(volumeInfo.IsDriftInPositiveX() ? -x0 : x0, 0.f, 0.f);

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pInputVertex->GetPosition() + volumeInfo.GetCenter() + xOffset;
    parameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    parameters.m_vertexType = pInputVertex->GetVertexType();

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*pAlgorithm, parameters, pNewVertex));

    return pNewVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *StitchingObjectCreationTool::CreatePfo(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pInputPfo,
    const ClusterList &newClusterList, const VertexList &newVertexList) const
{
    PandoraContentApi::ParticleFlowObject::Parameters parameters;
    parameters.m_particleId = pInputPfo->GetParticleId();
    parameters.m_charge = pInputPfo->GetCharge();
    parameters.m_mass = pInputPfo->GetMass();
    parameters.m_energy = pInputPfo->GetEnergy();
    parameters.m_momentum = pInputPfo->GetMomentum();
    parameters.m_clusterList = newClusterList;
    parameters.m_trackList.clear();
    parameters.m_vertexList = newVertexList;

    const ParticleFlowObject *pNewPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*pAlgorithm, parameters, pNewPfo));

    return pNewPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingObjectCreationTool::AddStitchingInfo(const ParticleFlowObject *const /*pNewPfo*/, const Pandora *const /*pPandora*/, StitchingInfo &/*stitchingInfo*/) const
{
    // TODO - work out what kind of information needs to be recorded
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingObjectCreationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NewClusterListName", m_newClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NewVertexListName", m_newVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NewPfoListName", m_newPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
