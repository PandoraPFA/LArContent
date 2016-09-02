/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the neutrino creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode NeutrinoCreationAlgorithm::Run()
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_vertexListName, 
        pVertexList));

    if (!pVertexList || pVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoCreationAlgorithm: unable to find vertex list " << m_vertexListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    std::string neutrinoPfoListName;
    const PfoList *pNeutrinoPfoList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNeutrinoPfoList, 
        neutrinoPfoListName));

    for (VertexList::const_iterator vIter = pVertexList->begin(), vIterEnd = pVertexList->end(); vIter != vIterEnd; ++vIter)
    {
        const Vertex *const pVertex = *vIter;

        if (VERTEX_3D != pVertex->GetVertexType())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // ATTN Placeholder properties to be set by a later neutrino properties algorithm
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = NU_MU;
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
        pfoParameters.m_vertexList.insert(pVertex);

        const ParticleFlowObject *pNeutrinoPfo = NULL;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNeutrinoPfo));
    }

    if ((NULL == pNeutrinoPfoList) || pNeutrinoPfoList->empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_neutrinoPfoListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_neutrinoPfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
