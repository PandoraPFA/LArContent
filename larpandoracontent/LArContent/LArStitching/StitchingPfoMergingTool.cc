/**
 *  @file   LArContent/src/LArStitching/StitchingPfoMergingTool.cc
 * 
 *  @brief  Implementation of the stitching pfo merging tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArContent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArContent/LArStitching/StitchingPfoMergingTool.h"

using namespace pandora;

namespace lar_content
{

void StitchingPfoMergingTool::Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    const PandoraInstanceList &pandoraInstances(MultiPandoraApi::GetDaughterPandoraInstanceList(&(this->GetPandora())));

    if (2 != pandoraInstances.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // ATTN Handy display for initial debug purposes
    const PfoList *pPfoList1(nullptr), *pPfoList2(nullptr), *pRecreatedPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pandoraInstances.at(0), pPfoList1));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pandoraInstances.at(1), pPfoList2));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pRecreatedPfoList));

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList1, "InputPfoList1", BLUE));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList2, "InputPfoList2", GREEN));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedPfoList, "RecreatedPfoList", RED));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    // ATTN Example, place-holder implementation
    const PfoList *pNewPfoList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pNewPfoList));

    PfoList primaryPfos;

    for (const ParticleFlowObject *const pPfo : *pNewPfoList)
        primaryPfos.insert(LArPfoHelper::GetParentPfo(pPfo));

    if (2 != primaryPfos.size())
        return;

    // ATTN Can add other information to StitchingInfo (prototype in StitchingAlgorithm.h; values are stored in StitchingObjectCreationTool)
    StitchingAlgorithm::PfoToVolumeIdMap &pfoToVolumeIdMap(stitchingInfo.m_pfoToVolumeIdMap);

    for (const ParticleFlowObject *const pPfo1 : primaryPfos)
    {
        for (const ParticleFlowObject *const pPfo2 : primaryPfos)
        {
            // ATTN Could ultimately check to see whether pfos are in adjacent volumes, not just different volumes
            if ((pPfo1 == pPfo2) || (pfoToVolumeIdMap.at(pPfo1) == pfoToVolumeIdMap.at(pPfo2)))
                continue;

            // ATTN Add key pattern-recognition decision #1 here. Write a function to return this boolean.
            const bool shouldMergePfos(true);

            if (shouldMergePfos)
            {
                // ATTN Add key pattern-recognition decision #2 here. Write a function to return this boolean.
                const ParticleFlowObject *const pPfoToEnlarge(pPfo1);
                const ParticleFlowObject *const pPfoToDelete(pPfo2);

                // ATTN Update stitching information here, just placeholder/example update shown below
                pfoToVolumeIdMap[pPfoToEnlarge] = -1;
                pfoToVolumeIdMap.erase(pPfoToDelete);

                PfoList enlargePfoList, deletePfoList;
                enlargePfoList.insert(pPfoToEnlarge);
                deletePfoList.insert(pPfoToDelete);
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &enlargePfoList, "enlargePfoList", BLUE));
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &deletePfoList, "deletePfoList", GREEN));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

                this->MergeAndDeletePfos(pAlgorithm, pPfoToEnlarge, pPfoToDelete);

                PfoList mergedPfoList;
                mergedPfoList.insert(pPfoToEnlarge);
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &mergedPfoList, "mergedPfoList", RED));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

                // ATTN Some decisions could be made in very different ways, e.g. which vertices to keep, hierarchy, etc. This is all just a suggestion.

                // Now need to break-out from this example approach as we have invalidated current iterators. Care with subsequent implementation.
                return;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingPfoMergingTool::MergeAndDeletePfos(const StitchingAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfoToEnlarge,
    const ParticleFlowObject *const pPfoToDelete) const
{
    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*pAlgorithm, pPfoToDelete));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*pAlgorithm, pPfoToEnlarge, pDaughterPfo));
    }

    for (const  Vertex *const pDaughterVertex : daughterVertices)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*pAlgorithm, pDaughterVertex));
    }

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *pParentCluster(this->GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pParentCluster, pDaughterCluster));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*pAlgorithm, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *StitchingPfoMergingTool::GetParentCluster(const ClusterList &clusterList, const HitType hitType) const
{
    unsigned int mostHits(0);
    const Cluster *pBestParentCluster(nullptr);

    for (const Cluster *const pParentCluster : clusterList)
    {
        if (hitType != LArClusterHelper::GetClusterHitType(pParentCluster))
            continue;

        const unsigned int nParentHits(pParentCluster->GetNCaloHits());

        if (nParentHits > mostHits)
        {
            mostHits = nParentHits;
            pBestParentCluster = pParentCluster;
        }
    }

    return pBestParentCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingPfoMergingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
