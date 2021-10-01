/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoDaughterVisibleVerticesAlgorithm.cc
 *
 *  @brief  Implementation of the cheating neutrino daughter vertices algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoDaughterVisibleVerticesAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingNeutrinoDaughterVisibleVerticesAlgorithm::CheatingNeutrinoDaughterVisibleVerticesAlgorithm() : m_collapseToPrimaryMCParticles(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoDaughterVisibleVerticesAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoListName, pPfoList));

    if (!pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingNeutrinoDaughterVisibleVerticesAlgorithm: pfo list unavailable." << std::endl;

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

void CheatingNeutrinoDaughterVisibleVerticesAlgorithm::GetMCPrimaryMap(LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const
{
    if (m_collapseToPrimaryMCParticles)
    {
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoDaughterVisibleVerticesAlgorithm::ProcessRecoNeutrinos(const PfoList &neutrinoPfos, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap) const
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

void CheatingNeutrinoDaughterVisibleVerticesAlgorithm::ProcessDaughterPfo(
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

    const CartesianVector mcVertexPosition = pMCParticle->GetVertex();
    CartesianVector closestPosition(mcVertexPosition);

    CaloHitList spacePoints;
    LArPfoHelper::GetCaloHits(pDaughterPfo, TPC_3D, spacePoints);

    std::cout << "Cheating neutrino visible vertices: " << spacePoints.size() << std::endl;

    float closestDistance(std::numeric_limits<float>::max());
                    
    for (const CaloHit *const spacePoint : spacePoints)
    {
        float distance((mcVertexPosition - spacePoint->GetPositionVector()).GetMagnitude());

        if (distance < closestDistance)
        {
            closestDistance = distance;
            closestPosition = spacePoint->GetPositionVector();
        }
    }

    const VertexList *pVertexList(nullptr);
    std::string vertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = closestPosition;
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

StatusCode CheatingNeutrinoDaughterVisibleVerticesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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
