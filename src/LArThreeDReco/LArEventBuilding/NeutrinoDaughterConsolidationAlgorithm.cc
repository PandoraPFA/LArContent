/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/NeutrinoDaughterConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the neutrino daughter consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArObjects/LArThreeDSlidingFitResult.h"

#include "LArThreeDReco/LArEventBuilding/NeutrinoDaughterConsolidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NeutrinoDaughterConsolidationAlgorithm::NeutrinoDaughterConsolidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoDaughterConsolidationAlgorithm::Run()
{
    const Vertex *pNeutrinoVertex(nullptr);
    const ParticleFlowObject *pNeutrinoPfo(nullptr);
    this->GetNeutrinoProperties(pNeutrinoPfo, pNeutrinoVertex);

    while (true)
    {
        PfoList parentCandidates, daughtersToMerge;
        this->ClassifyPrimaryParticles(pNeutrinoPfo, pNeutrinoVertex, parentCandidates, daughtersToMerge);

        PfoAssociationList pfoAssociationList;
        this->GetPfoAssociations(pNeutrinoVertex, parentCandidates, daughtersToMerge, pfoAssociationList);

        std::sort(pfoAssociationList.begin(), pfoAssociationList.end());
        const bool pfoMergeMade(this->ProcessPfoAssociations(pfoAssociationList));

        if (!pfoMergeMade)
            break;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterConsolidationAlgorithm::GetNeutrinoProperties(const ParticleFlowObject *&pNeutrinoPfo, const Vertex *&pNeutrinoVertex) const
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    PfoList neutrinoPfoList;
    VertexList neutrinoVertexList;

    for (const Pfo *const pPfo : *pPfoList)
    {
        if (!LArPfoHelper::IsNeutrino(pPfo))
            continue;

        neutrinoPfoList.insert(pPfo);
        neutrinoVertexList.insert(pPfo->GetVertexList().begin(), pPfo->GetVertexList().end());
    }

    if ((1 != neutrinoPfoList.size()) || (1 != neutrinoVertexList.size()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    pNeutrinoPfo = *(neutrinoPfoList.begin());
    pNeutrinoVertex = *(neutrinoVertexList.begin());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterConsolidationAlgorithm::ClassifyPrimaryParticles(const ParticleFlowObject *const pNeutrinoPfo, const Vertex *const pNeutrinoVertex,
    PfoList &parentCandidates, PfoList &daughtersToMerge) const
{
    const PfoList candidateList(pNeutrinoPfo->GetDaughterPfoList());

    for (const Pfo *const pCandidatePfo : candidateList)
    {
        const VertexList &vertexList(pCandidatePfo->GetVertexList());

        if (vertexList.empty())
            continue;

        if (1 != vertexList.size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const Vertex *const pCandidateVertex(*(vertexList.begin()));
        const float distance3D((pCandidateVertex->GetPosition() - pNeutrinoVertex->GetPosition()).GetMagnitude());

        if (distance3D > 5.f * 14.f) // TODO
        {
            daughtersToMerge.insert(pCandidatePfo);
        }
        else
        {
            parentCandidates.insert(pCandidatePfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterConsolidationAlgorithm::GetPfoAssociations(const Vertex *const pNeutrinoVertex, const PfoList &parentCandidates,
    const PfoList &daughtersToMerge, PfoAssociationList &pfoAssociationList) const
{
    const CartesianVector neutrinoVertexPosition(pNeutrinoVertex->GetPosition());
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (const Pfo *const pParentPfo : parentCandidates)
    {
        ClusterList threeDClusters;
        LArPfoHelper::GetThreeDClusterList(pParentPfo, threeDClusters);

        if (threeDClusters.empty())
            continue;

        if (1 != threeDClusters.size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        CartesianVector parentDirection(0.f, 0.f, 0.f);
        const float trackFactor(LArPfoHelper::IsTrack(pParentPfo) ? 0.001f : 1.f); //TODO

        try
        {
            const ThreeDSlidingFitResult slidingFitResult(*(threeDClusters.begin()), 20, slidingFitPitch); //TODO
            const bool useMinLayer((slidingFitResult.GetGlobalMinLayerPosition() - neutrinoVertexPosition).GetMagnitudeSquared() <
                (slidingFitResult.GetGlobalMaxLayerPosition() - neutrinoVertexPosition).GetMagnitudeSquared());
            parentDirection = useMinLayer ? slidingFitResult.GetGlobalMinLayerDirection() : slidingFitResult.GetGlobalMaxLayerDirection();
        }
        catch (StatusCodeException &)
        {
            continue;
        }

        for (const Pfo *const pDaughterPfo : daughtersToMerge)
        {
            if (1 != pDaughterPfo->GetVertexList().size())
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            const Vertex *const pDaughterVertex(*(pDaughterPfo->GetVertexList().begin()));
            const CartesianVector vertexSeparation(pDaughterVertex->GetPosition() - neutrinoVertexPosition);

            const float rL(vertexSeparation.GetDotProduct(parentDirection));
            const float rT(vertexSeparation.GetCrossProduct(parentDirection).GetMagnitude());
            const float score(std::exp(-rL / 14.f) * std::exp(-rT / 10.f) * trackFactor); //TODO
            pfoAssociationList.push_back(PfoAssociation(pParentPfo, pDaughterPfo, score));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoDaughterConsolidationAlgorithm::ProcessPfoAssociations(const PfoAssociationList &pfoAssociationList) const
{
    for (const PfoAssociation &pfoAssociation : pfoAssociationList)
    {
        //if (score)
        //    continue;
        this->MergeAndDeletePfos(pfoAssociation.GetParentPfo(), pfoAssociation.GetDaughterPfo());
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterConsolidationAlgorithm::MergeAndDeletePfos(const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete) const
{
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

const Cluster *NeutrinoDaughterConsolidationAlgorithm::GetParentCluster(const ClusterList &clusterList, const HitType hitType) const
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
const std::string &NeutrinoDaughterConsolidationAlgorithm::GetListName(const T *const pT) const
{
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

StatusCode NeutrinoDaughterConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NeutrinoPfoListName", m_neutrinoPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "DaughterListNames", m_daughterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
