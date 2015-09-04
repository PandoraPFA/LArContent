/**
 *  @file   LArContent/src/LArUtility/NeutrinoListManagementAlgorithm.cc
 * 
 *  @brief  Implementation of the neutrino list management algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPfoHelper.h"

#include "LArUtility/NeutrinoListManagementAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode NeutrinoListManagementAlgorithm::Run()
{
    this->PfoListManagement();
    this->VertexListManagement();
    this->ClusterListManagement();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoListManagementAlgorithm::PfoListManagement() const
{
    const PfoList *pPfoList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoListManagementAlgorithm: unable to find current pfo list " << std::endl;

        return;
    }

    PfoList neutrinoPfos, trackPfos, showerPfos;

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            neutrinoPfos.insert(pPfo);
        }
        else if (LArPfoHelper::IsTrack(pPfo))
        {
            trackPfos.insert(pPfo);
        }
        else if (LArPfoHelper::IsShower(pPfo))
        {
            showerPfos.insert(pPfo);
        }
    }

    if (neutrinoPfos.empty() || (trackPfos.empty() && showerPfos.empty()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_neutrinoPfoListName, neutrinoPfos));

    if (!trackPfos.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, trackPfos));

    if (!showerPfos.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, showerPfos));

    if (!pPfoList->empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoListManagementAlgorithm::VertexListManagement() const
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    if (!pVertexList || pVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoListManagementAlgorithm: unable to find current vertex list " << std::endl;

        return;
    }

    if (1 != pVertexList->size())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_neutrinoVertexListName));

    if (!pVertexList->empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoListManagementAlgorithm::ClusterListManagement() const
{
    const ClusterList *pClusterList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    
    if (!pClusterList || pClusterList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoListManagementAlgorithm: unable to find current cluster list " << std::endl;

        return;
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_clusterListName));

    if (!pClusterList->empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoListManagementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NeutrinoPfoListName", m_neutrinoPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackPfoListName", m_trackPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerPfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NeutrinoVertexListName", m_neutrinoVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListName", m_clusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
