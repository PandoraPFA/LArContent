/**
 *  @file   LArContent/src/LArCheating/CheatingNeutrinoCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating neutrino creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArCheating/CheatingNeutrinoCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingNeutrinoCreationAlgorithm::CheatingNeutrinoCreationAlgorithm() :
    m_vertexTolerance(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoCreationAlgorithm::Run()
{
    const MCParticle *pMCNeutrino(NULL);
    this->GetMCNeutrino(pMCNeutrino);

    const ParticleFlowObject *pNeutrinoPfo(NULL);
    this->CreateAndSaveNeutrinoPfo(pMCNeutrino, pNeutrinoPfo);

    this->AddNeutrinoVertex(pMCNeutrino, pNeutrinoPfo);
    this->AddDaughterPfos(pMCNeutrino, pNeutrinoPfo);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::GetMCNeutrino(const MCParticle *&pMCNeutrino) const
{
    pMCNeutrino = NULL;

    const MCParticleList *pMCParticleList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    MCParticleVector mcNeutrinoList;
    LArMCParticleHelper::GetNeutrinoMCParticleList(pMCParticleList, mcNeutrinoList);
    pMCNeutrino = ((1 == mcNeutrinoList.size()) ? *(mcNeutrinoList.begin()) : NULL);

    if (!pMCNeutrino || (MC_3D != pMCNeutrino->GetMCParticleType()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::CreateAndSaveNeutrinoPfo(const MCParticle *const pMCNeutrino, const ParticleFlowObject *&pNeutrinoPfo) const
{
    pNeutrinoPfo = NULL;

    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = pMCNeutrino->GetParticleId();
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = pMCNeutrino->GetEnergy();
    pfoParameters.m_momentum = pMCNeutrino->GetMomentum();

    std::string neutrinoPfoListName;
    const PfoList *pNeutrinoPfoList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNeutrinoPfoList, neutrinoPfoListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNeutrinoPfo));

    if ((NULL == pNeutrinoPfoList) || pNeutrinoPfoList->empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_neutrinoPfoListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_neutrinoPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::AddNeutrinoVertex(const MCParticle *const pMCNeutrino, const ParticleFlowObject *const pNeutrinoPfo) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pNeutrinoVertex((1 == pVertexList->size()) ? *(pVertexList->begin()) : NULL);

    if (!pNeutrinoVertex || (VERTEX_3D != pNeutrinoVertex->GetVertexType()) || ((pNeutrinoVertex->GetPosition() - pMCNeutrino->GetEndpoint()).GetMagnitude() > m_vertexTolerance))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pNeutrinoPfo, pNeutrinoVertex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoCreationAlgorithm::AddDaughterPfos(const MCParticle *const pMCNeutrino, const ParticleFlowObject *const pNeutrinoPfo) const
{
    for (const std::string &daughterPfoListName : m_daughterPfoListNames)
    {
        const PfoList *pDaughterPfoList(NULL);

        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, daughterPfoListName, pDaughterPfoList))
        {
            for (const ParticleFlowObject *const pDaughterPfo : *pDaughterPfoList)
            {
                bool isDaughter(true);
                ClusterList twoDClusterList;
                LArPfoHelper::GetTwoDClusterList(pDaughterPfo, twoDClusterList);

                for (const Cluster *const pTwoDCluster : twoDClusterList)
                {
                    try
                    {
                        // TODO consider hierarchy of daughter particles if clustering hasn't already collapsed to mc primaries
                        if (pMCNeutrino != LArMCParticleHelper::GetParentNeutrino(MCParticleHelper::GetMainMCParticle(pTwoDCluster)))
                            isDaughter = false;
                    }
                    catch (StatusCodeException &)
                    {
                        isDaughter = false;
                    }
                }

                if (isDaughter)
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pDaughterPfo));
            }
        }
        else if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        {
            std::cout << "CheatingNeutrinoCreationAlgorithm: pfo list " << daughterPfoListName << " unavailable." << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NeutrinoPfoListName", m_neutrinoPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "DaughterPfoListNames", m_daughterPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexTolerance", m_vertexTolerance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
