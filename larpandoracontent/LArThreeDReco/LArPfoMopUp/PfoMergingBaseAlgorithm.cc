/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/PfoMergingBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the pfo merging algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/PfoMergingBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void PfoMergingBaseAlgorithm::MergeAndDeletePfos(const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete) const
{
    if (pPfoToEnlarge == pPfoToDelete)
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, this->GetListName(pPfoToDelete)));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoToEnlarge, pDaughterPfo));
    }

    for (const  Vertex *const pDaughterVertex : daughterVertices)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pDaughterVertex, this->GetListName(pDaughterVertex)));
    }

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *pParentCluster(this->GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster,
                this->GetListName(pParentCluster), this->GetListName(pDaughterCluster)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *PfoMergingBaseAlgorithm::GetParentCluster(const ClusterList &clusterList, const HitType hitType) const
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

template <typename T>
const std::string PfoMergingBaseAlgorithm::GetListName(const T *const pT) const
{
    std::string currentListName;
    const std::MANAGED_CONTAINER<const T*> *pCurrentList(nullptr);
    (void) PandoraContentApi::GetCurrentList(*this, pCurrentList, currentListName);

    if (pCurrentList && (pCurrentList->count(pT)))
        return currentListName;

    for (const std::string &listName : m_daughterListNames)
    {
        const std::MANAGED_CONTAINER<const T*> *pList(nullptr);
        (void) PandoraContentApi::GetList(*this, listName, pList);

        if (pList && (pList->count(pT)))
            return listName;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoMergingBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "DaughterListNames", m_daughterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
