/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating neutrino creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingNeutrinoCreationAlgorithm::CheatingNeutrinoCreationAlgorithm() :
    m_collapseToPrimaryMCParticles(false),
    m_vertexTolerance(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoCreationAlgorithm::Run()
{
    MCParticleVector mcNeutrinoVector;
    this->GetMCNeutrinoVector(mcNeutrinoVector);

    for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
    {
        const ParticleFlowObject *pNeutrinoPfo(nullptr);
        this->CreateAndSaveNeutrinoPfo(pMCNeutrino, pNeutrinoPfo);

        if (!pNeutrinoPfo)
            return STATUS_CODE_FAILURE;

        this->AddNeutrinoVertex(pMCNeutrino, pNeutrinoPfo);

        MCParticleToPfoMap mcParticleToPfoMap;
        this->GetMCParticleToDaughterPfoMap(mcParticleToPfoMap);
        this->CreatePfoHierarchy(pMCNeutrino, pNeutrinoPfo, mcParticleToPfoMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::GetMCNeutrinoVector(MCParticleVector &mcNeutrinoVector) const
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    MCParticleVector allMCNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, allMCNeutrinoVector);

    for (const MCParticle *const pMCNeutrino : allMCNeutrinoVector)
    {
        if (LArMCParticleHelper::IsNeutrino(pMCNeutrino) && (MC_3D == pMCNeutrino->GetMCParticleType()))
            mcNeutrinoVector.push_back(pMCNeutrino);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::CreateAndSaveNeutrinoPfo(const MCParticle *const pMCNeutrino, const ParticleFlowObject *&pNeutrinoPfo) const
{
    pNeutrinoPfo = nullptr;

    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = pMCNeutrino->GetParticleId();
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = pMCNeutrino->GetEnergy();
    pfoParameters.m_momentum = pMCNeutrino->GetMomentum();

    std::string neutrinoPfoListName;
    const PfoList *pNeutrinoPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNeutrinoPfoList, neutrinoPfoListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNeutrinoPfo));

    if (!pNeutrinoPfoList || pNeutrinoPfoList->empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_neutrinoPfoListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_neutrinoPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::AddNeutrinoVertex(const MCParticle *const pMCNeutrino, const ParticleFlowObject *const pNeutrinoPfo) const
{
    const VertexList *pVertexList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    const Vertex *pNeutrinoVertex(nullptr);
    float closestVertexDistance(std::numeric_limits<float>::max());

    for (const Vertex *const pVertex : *pVertexList)
    {
        const float distance((pVertex->GetPosition() - pMCNeutrino->GetEndpoint()).GetMagnitude());

        if (distance < closestVertexDistance)
        {
            pNeutrinoVertex = pVertex;
            closestVertexDistance = distance;
        }
    }

    if (!pNeutrinoVertex || (VERTEX_3D != pNeutrinoVertex->GetVertexType()) ||
        ((pNeutrinoVertex->GetPosition() - pMCNeutrino->GetEndpoint()).GetMagnitude() > m_vertexTolerance))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pNeutrinoPfo, pNeutrinoVertex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::GetMCParticleToDaughterPfoMap(MCParticleToPfoMap &mcParticleToPfoMap) const
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;

    if (m_collapseToPrimaryMCParticles)
    {
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    }

    for (const std::string &daughterPfoListName : m_daughterPfoListNames)
    {
        const PfoList *pDaughterPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, daughterPfoListName, pDaughterPfoList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CheatingNeutrinoCreationAlgorithm: pfo list " << daughterPfoListName << " unavailable." << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pDaughterPfo : *pDaughterPfoList)
        {
            const MCParticle *pMCParticle(LArMCParticleHelper::GetMainMCParticle(pDaughterPfo));

            if (m_collapseToPrimaryMCParticles)
            {
                LArMCParticleHelper::MCRelationMap::const_iterator primaryIter = mcPrimaryMap.find(pMCParticle);

                if (mcPrimaryMap.end() == primaryIter)
                    throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                pMCParticle = primaryIter->second;
            }

            if (!mcParticleToPfoMap.insert(MCParticleToPfoMap::value_type(pMCParticle, pDaughterPfo)).second)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::CreatePfoHierarchy(const MCParticle *const pParentMCParticle,
    const ParticleFlowObject *const pParentPfo, const MCParticleToPfoMap &mcParticleToPfoMap) const
{
    for (const MCParticle *const pDaughterMCParticle : pParentMCParticle->GetDaughterList())
    {
        MCParticleToPfoMap::const_iterator mapIter = mcParticleToPfoMap.find(pDaughterMCParticle);

        if (mcParticleToPfoMap.end() == mapIter)
        {
            this->CreatePfoHierarchy(pDaughterMCParticle, pParentPfo, mcParticleToPfoMap);
            continue;
        }

        const ParticleFlowObject *const pDaughterPfo(mapIter->second);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pParentPfo, pDaughterPfo));
        this->CreatePfoHierarchy(pDaughterMCParticle, pDaughterPfo, mcParticleToPfoMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CollapseToPrimaryMCParticles", m_collapseToPrimaryMCParticles));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DaughterPfoListNames", m_daughterPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexTolerance", m_vertexTolerance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
