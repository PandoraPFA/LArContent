/**
 *  @file   larpandoracontent/LArStitching/StitchingObjectCreationTool.cc
 *
 *  @brief  Implementation of the stitching object creation tool class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArStitching/StitchingObjectCreationTool.h"

using namespace pandora;
using namespace lar_content;

namespace lar_content
{

StitchingObjectCreationTool::StitchingObjectCreationTool() :
    m_recreateTwoDContent(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingObjectCreationTool::Run(const StitchingAlgorithm *const pAlgorithm, StitchingInfo &stitchingInfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

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
    {
        this->Recreate3DContent(pAlgorithm, pPandora, MultiPandoraApi::GetVolumeInfo(pPandora), stitchingInfo);
    }

    if (!pClusterList->empty())
    {
        const std::string &newListName(pAlgorithm->GetNewClusterListName());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*pAlgorithm, newListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, newListName));
    }

    if (!pVertexList->empty())
    {
        const std::string &newListName(pAlgorithm->GetNewVertexListName());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*pAlgorithm, newListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*pAlgorithm, newListName));
    }

    if (!pPfoList->empty())
    {
        const std::string &newListName(pAlgorithm->GetNewPfoListName());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*pAlgorithm, newListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*pAlgorithm, newListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingObjectCreationTool::Recreate3DContent(const StitchingAlgorithm *const pAlgorithm, const Pandora *const pPandora, const VolumeInfo &volumeInfo,
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

void StitchingObjectCreationTool::Recreate3DContent(const StitchingAlgorithm *const pAlgorithm, const ParticleFlowObject *const pInputPfo,
    const ParticleFlowObject *const pNewParentPfo, const Pandora *const pPandora, const VolumeInfo &volumeInfo, StitchingInfo &stitchingInfo) const
{
    float x0(0.f);

    try
    {
        x0 = volumeInfo.GetParticleX0(pInputPfo);
    }
    catch (pandora::StatusCodeException& )
    {
    }

    ClusterList inputClusterList, newClusterList;

    if (!m_recreateTwoDContent)
    {
        LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList);
    }
    else
    {
        inputClusterList = pInputPfo->GetClusterList();
    }

    for (const Cluster *const pInputCluster : inputClusterList)
    {
        CaloHitList inputCaloHitList, newCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
          newCaloHitList.push_back(pAlgorithm->CreateCaloHit(pInputCaloHit, volumeInfo, x0));

        if (!newCaloHitList.empty())
            newClusterList.push_back(pAlgorithm->CreateCluster(pInputCluster, newCaloHitList));
    }

    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
    {
        if (VERTEX_3D == pInputVertex->GetVertexType())
          newVertexList.push_back(pAlgorithm->CreateVertex(pInputVertex, volumeInfo, x0));
    }

    const ParticleFlowObject *const pNewPfo = pAlgorithm->CreatePfo(pInputPfo, newClusterList, newVertexList);
    this->AddStitchingInfo(pNewPfo, pPandora, volumeInfo, stitchingInfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*pAlgorithm, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate3DContent(pAlgorithm, pInputDaughterPfo, pNewPfo, pPandora, volumeInfo, stitchingInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingObjectCreationTool::AddStitchingInfo(const ParticleFlowObject *const pNewPfo, const Pandora *const /*pPandora*/,
    const VolumeInfo &volumeInfo, StitchingInfo &stitchingInfo) const
{
    // TODO - work out what kind of information needs to be recorded
    if (!stitchingInfo.m_pfoToVolumeIdMap.insert(StitchingAlgorithm::PfoToVolumeIdMap::value_type(pNewPfo, volumeInfo.GetIdNumber())).second)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingObjectCreationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RecreateTwoDContent", m_recreateTwoDContent));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
