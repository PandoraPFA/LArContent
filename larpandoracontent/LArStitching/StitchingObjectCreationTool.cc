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

StitchingObjectCreationTool::StitchingObjectCreationTool()
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
        this->Recreate3DContent(pAlgorithm, pPandora, pPandora->GetGeometry()->GetLArTPC(), stitchingInfo);
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

void StitchingObjectCreationTool::Recreate3DContent(const StitchingAlgorithm *const pAlgorithm, const Pandora *const pPandora, const LArTPC &larTPC,
    StitchingInfo &stitchingInfo) const
{
    const PfoList *pPfoList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pPandora, pPfoList));

    for (const ParticleFlowObject *const pInputPfo : *pPfoList)
    {
        if (!pInputPfo->GetParentPfoList().empty())
            continue;

        this->Recreate3DContent(pAlgorithm, pInputPfo, nullptr, larTPC, stitchingInfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingObjectCreationTool::Recreate3DContent(const StitchingAlgorithm *const pAlgorithm, const ParticleFlowObject *const pInputPfo,
    const ParticleFlowObject *const pNewParentPfo, const LArTPC &larTPC, StitchingInfo &stitchingInfo) const
{
    // Assume zero input/offset X0
    const float x0(0.f);

    // Recreate clusters
    typedef std::map<const CaloHit*, const CaloHit*> CaloHitMap;
    CaloHitMap newParentAddresses;

    ClusterList inputClusterList2D, inputClusterList3D, newClusterList;

    if (pAlgorithm->RecreateTwoDContent())
        LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);

    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);

    for (const Cluster *const pInputCluster : inputClusterList2D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pNewCaloHit = pAlgorithm->CreateCaloHit(pInputCaloHit, larTPC, x0);
            newCaloHitList.push_back(pNewCaloHit);
            newParentAddresses.insert(CaloHitMap::value_type(pInputCaloHit, pNewCaloHit));
        }

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
        {
            const CaloHit *const pNewCaloHit = pAlgorithm->CreateCaloHit(pInputCaloHit, larTPC, x0);
            newIsolatedCaloHitList.push_back(pNewCaloHit);
            newParentAddresses.insert(CaloHitMap::value_type(pInputCaloHit, pNewCaloHit));
        }

        if (!newCaloHitList.empty())
            newClusterList.push_back(pAlgorithm->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    for (const Cluster *const pInputCluster : inputClusterList3D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pParentCaloHit = static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress());

            CaloHitMap::const_iterator iter = newParentAddresses.find(pParentCaloHit);
            if (pAlgorithm->RecreateTwoDContent() && newParentAddresses.end() == iter)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            const CaloHit *const pNewParentCaloHit = (pAlgorithm->RecreateTwoDContent() ? iter->second : pParentCaloHit);
            newCaloHitList.push_back(pAlgorithm->CreateCaloHit(pInputCaloHit, pNewParentCaloHit, larTPC, x0));
        }

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
        {
            const CaloHit *const pParentCaloHit = static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress());

            CaloHitMap::const_iterator iter = newParentAddresses.find(pParentCaloHit);
            if (pAlgorithm->RecreateTwoDContent() && newParentAddresses.end() == iter)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            const CaloHit *const pNewParentCaloHit = (pAlgorithm->RecreateTwoDContent() ? iter->second : pParentCaloHit);
            newIsolatedCaloHitList.push_back(pAlgorithm->CreateCaloHit(pInputCaloHit, pNewParentCaloHit, larTPC, x0));
        }

        if (!newCaloHitList.empty())
            newClusterList.push_back(pAlgorithm->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }


    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
    {
        if (VERTEX_3D == pInputVertex->GetVertexType())
          newVertexList.push_back(pAlgorithm->CreateVertex(pInputVertex, larTPC, x0));
    }

    const ParticleFlowObject *const pNewPfo = pAlgorithm->CreatePfo(pInputPfo, newClusterList, newVertexList);
    this->AddStitchingInfo(pNewPfo, larTPC, stitchingInfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*pAlgorithm, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate3DContent(pAlgorithm, pInputDaughterPfo, pNewPfo, larTPC, stitchingInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingObjectCreationTool::AddStitchingInfo(const ParticleFlowObject *const pNewPfo, const LArTPC &larTPC, StitchingInfo &stitchingInfo) const
{
    if (!stitchingInfo.m_pfoToLArTPCMap.insert(StitchingAlgorithm::PfoToLArTPCMap::value_type(pNewPfo, &larTPC)).second)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingObjectCreationTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
