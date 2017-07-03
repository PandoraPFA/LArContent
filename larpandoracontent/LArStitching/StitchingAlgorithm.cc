/**
 *  @file   larpandoracontent/LArStitching/StitchingAlgorithm.cc
 *
 *  @brief  Implementation of the Stitching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

#include "larpandoracontent/LArStitching/StitchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StitchingAlgorithm::StitchingAlgorithm() :
    m_recreateTwoDContent(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingAlgorithm::Run()
{
    StitchingInfo stitchingInfo;

    for (StitchingTool *const pStitchingTool : m_algorithmToolVector)
    {
        pStitchingTool->Run(this, stitchingInfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *StitchingAlgorithm::CreateCaloHit(const CaloHit *const pInputCaloHit, const VolumeInfo &volumeInfo, const float x0) const
{
    return this->CreateCaloHit(pInputCaloHit, nullptr, volumeInfo, x0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *StitchingAlgorithm::CreateCaloHit(const CaloHit *const pInputCaloHit, const CaloHit *const pParentCaloHit,
    const VolumeInfo &volumeInfo, const float x0) const
{
    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = LArStitchingHelper::GetCorrectedPosition(volumeInfo, x0, pInputCaloHit->GetPositionVector());
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

    if (nullptr == pParentCaloHit)
    {
        parameters.m_pParentAddress = pInputCaloHit->GetParentAddress();
    }
    else
    {
        parameters.m_pParentAddress = static_cast<const void*>(pParentCaloHit);
    }

    const CaloHit *pNewCaloHit(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pNewCaloHit));

    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = pInputCaloHit->IsIsolated();
    metadata.m_isPossibleMip = pInputCaloHit->IsPossibleMip();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pNewCaloHit, metadata));

    return pNewCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *StitchingAlgorithm::CreateCluster(const Cluster *const pInputCluster, const CaloHitList &newCaloHitList,
    const CaloHitList &newIsolatedCaloHitList) const
{
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = newCaloHitList;
    parameters.m_isolatedCaloHitList = newIsolatedCaloHitList;

    const Cluster *pNewCluster(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = pInputCluster->GetParticleId();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pNewCluster, metadata));

    return pNewCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *StitchingAlgorithm::CreateVertex(const Vertex *const pInputVertex, const VolumeInfo &volumeInfo, const float x0) const
{
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = LArStitchingHelper::GetCorrectedPosition(volumeInfo, x0, pInputVertex->GetPosition());
    parameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    parameters.m_vertexType = pInputVertex->GetVertexType();

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));

    return pNewVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *StitchingAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo,
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
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, parameters, pNewPfo));

    return pNewPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingAlgorithm::ShiftPfoHierarchy(const ParticleFlowObject *const pParentPfo, const StitchingAlgorithm::StitchingInfo &stitchingInfo,
    const float x0) const
{
    // Verify that this is a top-level Pfo
    if (!pParentPfo->GetParentPfoList().empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Get the volume information for this Pfo
    const StitchingAlgorithm::PfoToVolumeIdMap &pfoToVolumeIdMap(stitchingInfo.m_pfoToVolumeIdMap);
    PfoToVolumeIdMap::const_iterator iter = pfoToVolumeIdMap.find(pParentPfo);

    if (pfoToVolumeIdMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const VolumeInfo &volumeInfo(MultiPandoraApi::GetVolumeInfo(&this->GetPandora(), iter->second));

    // Shift the Pfo hierarchy
    PfoList pfoList;
    LArPfoHelper::GetAllDownstreamPfos(pParentPfo, pfoList);

    for (const ParticleFlowObject *const pDaughterPfo : pfoList)
        this->ShiftPfo(pDaughterPfo, volumeInfo, x0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingAlgorithm::ShiftPfo(const ParticleFlowObject *const pInputPfo, const VolumeInfo &volumeInfo, const float x0) const
{
    // Shift x positions of all calo hits associated with a Pfo. We do this by recreating the calo hits in their shifted positions.
    // ATTN: need to keep careful track of parent addresses when recreating both 2D and 3D hits

    typedef std::map<const CaloHit*, const CaloHit*> CaloHitMap;
    CaloHitMap newParentAddresses;

    ClusterList inputClusterList2D, inputClusterList3D;

    LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);

    for (const Cluster *const pInputCluster : inputClusterList2D)
    {
        CaloHitList inputCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pNewCaloHit = this->CreateCaloHit(pInputCaloHit, volumeInfo, x0);
            newParentAddresses.insert(CaloHitMap::value_type(pInputCaloHit, pNewCaloHit));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pInputCluster, pNewCaloHit));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pInputCluster, pInputCaloHit));
        }

        const CaloHitList isolatedCaloHitList(pInputCluster->GetIsolatedCaloHitList().begin(), pInputCluster->GetIsolatedCaloHitList().end());

        for (const CaloHit *const pInputCaloHit : isolatedCaloHitList)
        {
            const CaloHit *const pNewCaloHit = this->CreateCaloHit(pInputCaloHit, volumeInfo, x0);
            newParentAddresses.insert(CaloHitMap::value_type(pInputCaloHit, pNewCaloHit));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedToCluster(*this, pInputCluster, pNewCaloHit));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveIsolatedFromCluster(*this, pInputCluster, pInputCaloHit));
        }
    }

    for (const Cluster *const pInputCluster : inputClusterList3D)
    {
        CaloHitList inputCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pParentCaloHit = static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress());

            CaloHitMap::const_iterator iter = newParentAddresses.find(pParentCaloHit);
            if (m_recreateTwoDContent && newParentAddresses.end() == iter)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            const CaloHit *const pNewParentCaloHit = (m_recreateTwoDContent ? iter->second : pParentCaloHit);
            const CaloHit *const pNewCaloHit = this->CreateCaloHit(pInputCaloHit, pNewParentCaloHit, volumeInfo, x0);

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pInputCluster, pNewCaloHit));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pInputCluster, pInputCaloHit));
        }

        const CaloHitList isolatedCaloHitList(pInputCluster->GetIsolatedCaloHitList().begin(), pInputCluster->GetIsolatedCaloHitList().end());

        for (const CaloHit *const pInputCaloHit : isolatedCaloHitList)
        {
            const CaloHit *const pParentCaloHit = static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress());

            CaloHitMap::const_iterator iter = newParentAddresses.find(pParentCaloHit);
            if (m_recreateTwoDContent && newParentAddresses.end() == iter)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            const CaloHit *const pNewParentCaloHit = (m_recreateTwoDContent ? iter->second : pParentCaloHit);
            const CaloHit *const pNewCaloHit = this->CreateCaloHit(pInputCaloHit, pNewParentCaloHit, volumeInfo, x0);

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedToCluster(*this, pInputCluster, pNewCaloHit));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveIsolatedFromCluster(*this, pInputCluster, pInputCaloHit));
        }
    }

    // Shift x positions of all vertices (by creating a new vertex and deleting the old vertex)
    const VertexList inputVertexList(pInputPfo->GetVertexList().begin(), pInputPfo->GetVertexList().end());

    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    for (const Vertex *const pInputVertex : inputVertexList)
    {
        const Vertex *const pNewVertex = this->CreateVertex(pInputVertex, volumeInfo, x0);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pInputPfo, pNewVertex));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo<Vertex>(*this, pInputPfo, pInputVertex));
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, this->GetNewVertexListName()));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, this->GetNewVertexListName()));

    for (const Vertex *const pInputVertex : inputVertexList)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pInputVertex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingAlgorithm::StitchPfos(const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete,
    StitchingAlgorithm::StitchingInfo &stitchingInfo) const
{
    // Get the mapping between Pfos and their drift volume IDs
    StitchingAlgorithm::PfoToVolumeIdMap &pfoToVolumeIdMap(stitchingInfo.m_pfoToVolumeIdMap);

    // Update stitching information here (ATTN: We assume that -1 is an invalid volume ID)
    pfoToVolumeIdMap[pPfoToEnlarge] = -1;
    pfoToVolumeIdMap.erase(pPfoToDelete);

    this->MergeAndDeletePfos(pPfoToEnlarge, pPfoToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "StitchingTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        StitchingTool *const pStitchingTool(dynamic_cast<StitchingTool*>(*iter));

        if (nullptr == pStitchingTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pStitchingTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "RecreateTwoDContent", m_recreateTwoDContent));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NewClusterListName", m_newClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NewVertexListName", m_newVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NewPfoListName", m_newPfoListName));

    if (m_newClusterListName.empty() || m_newVertexListName.empty() || m_newPfoListName.empty())
    {
        std::cout << "StitchingAlgorithm::ReadSettings - Invalid list name for new/recreated objects." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    m_daughterListNames.push_back(m_newClusterListName);
    m_daughterListNames.push_back(m_newVertexListName);
    m_daughterListNames.push_back(m_newPfoListName);

    return PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
