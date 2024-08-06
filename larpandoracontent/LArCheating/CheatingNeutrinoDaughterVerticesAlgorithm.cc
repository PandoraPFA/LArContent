/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoDaughterVerticesAlgorithm.cc
 *
 *  @brief  Implementation of the cheating neutrino daughter vertices algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoDaughterVerticesAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingNeutrinoDaughterVerticesAlgorithm::CheatingNeutrinoDaughterVerticesAlgorithm() :
    m_collapseToPrimaryMCParticles(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoDaughterVerticesAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoListName, pPfoList));

    if (!pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingNeutrinoDaughterVerticesAlgorithm: pfo list unavailable." << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    LArMCParticleHelper::MCRelationMap mcPrimaryMap;
    this->GetMCPrimaryMap(mcPrimaryMap);

    PfoList neutrinoPfos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, neutrinoPfos);

    this->ProcessRecoNeutrinos(neutrinoPfos, mcPrimaryMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoDaughterVerticesAlgorithm::GetMCPrimaryMap(LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const
{
    if (m_collapseToPrimaryMCParticles)
    {
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoDaughterVerticesAlgorithm::ProcessRecoNeutrinos(const PfoList &neutrinoPfos, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const
{
    for (const ParticleFlowObject *const pNeutrinoPfo : neutrinoPfos)
    {
        PfoList daughterPfos;
        LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, daughterPfos);

        PfoList::iterator neutrinoIter(std::find(daughterPfos.begin(), daughterPfos.end(), pNeutrinoPfo));

        if (daughterPfos.end() != neutrinoIter)
            daughterPfos.erase(neutrinoIter);

        for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
        {
            try
            {
                this->ProcessDaughterPfo(pDaughterPfo, mcPrimaryMap);
            }
            catch (const StatusCodeException &)
            {
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoDaughterVerticesAlgorithm::ProcessDaughterPfo(
    const ParticleFlowObject *const pDaughterPfo, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const
{
    const MCParticle *pMCParticle(LArMCParticleHelper::GetMainMCParticle(pDaughterPfo));

    if (m_collapseToPrimaryMCParticles)
    {
        LArMCParticleHelper::MCRelationMap::const_iterator primaryIter = mcPrimaryMap.find(pMCParticle);

        if (mcPrimaryMap.end() == primaryIter)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        pMCParticle = primaryIter->second;
    }

    const VertexList *pVertexList(nullptr);
    std::string vertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pMCParticle->GetVertex();
    parameters.m_vertexLabel = VERTEX_INTERACTION;
    parameters.m_vertexType = VERTEX_3D;

    const Vertex *pVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));

    if (!pVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_vertexListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pDaughterPfo, pVertex));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoDaughterVerticesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CollapseToPrimaryMCParticles", m_collapseToPrimaryMCParticles));

    if (m_collapseToPrimaryMCParticles)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_vertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
