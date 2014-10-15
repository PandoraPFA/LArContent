/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/NeutrinoBuildingAlgorithm.cc
 * 
 *  @brief  Implementation of the neutrino building algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArThreeDReco/LArEventBuilding/NeutrinoBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode NeutrinoBuildingAlgorithm::Run()
{
    try
    {
        PfoList daughterPfoList;
        this->GetDaughterPfoList(daughterPfoList);

        ParticleFlowObject *pNeutrinoPfo(NULL);
        this->CreateNeutrinoPfo(pNeutrinoPfo);

        this->AddDaughtersAndSetId(pNeutrinoPfo, daughterPfoList);
    }
    catch (StatusCodeException &)
    {
        std::cout << "NeutrinoBuildingAlgorithm: unable to build neutrino." << std::endl;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoBuildingAlgorithm::GetDaughterPfoList(PfoList &pfoList) const
{
    for (StringVector::const_iterator iter = m_daughterPfoListNames.begin(), iterEnd = m_daughterPfoListNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList(NULL);

        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, *iter, pPfoList))
        {
            pfoList.insert(pPfoList->begin(), pPfoList->end());
        }
        else
        {
            std::cout << "NeutrinoBuildingAlgorithm: pfo list " << *iter << " unavailable." << std::endl;
        }
    }

    if (pfoList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoBuildingAlgorithm::CreateNeutrinoPfo(ParticleFlowObject *&pNeutrinoPfo) const
{
    // TODO Correct these placeholder parameters
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = NU_MU;
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);

    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    Vertex *const pVertex((1 == pVertexList->size()) ? *(pVertexList->begin()) : NULL);

    if ((NULL != pVertex) && (VERTEX_3D == pVertex->GetVertexType()))
        pfoParameters.m_vertexList.insert(pVertex);

    std::string neutrinoPfoListName;
    const PfoList *pNeutrinoPfoList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNeutrinoPfoList, neutrinoPfoListName));

    pNeutrinoPfo = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNeutrinoPfo));

    if ((NULL == pNeutrinoPfoList) || pNeutrinoPfoList->empty() || (NULL == pNeutrinoPfo))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoBuildingAlgorithm::AddDaughtersAndSetId(ParticleFlowObject *const pNeutrinoPfo, const PfoList &daughterPfoList) const
{
    unsigned int nPrimaryTwoDHits(0);
    ParticleFlowObject *pPrimaryDaughter(NULL);

    for (PfoList::const_iterator iter = daughterPfoList.begin(), iterEnd = daughterPfoList.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pDaughterPfo(*iter);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pDaughterPfo))

        const unsigned int nTwoDHits(this->GetNTwoDHitsInPfo(pDaughterPfo));

        if (!pPrimaryDaughter || (nTwoDHits > nPrimaryTwoDHits))
        {
            nPrimaryTwoDHits = nTwoDHits;
            pPrimaryDaughter = pDaughterPfo;
        }
    }

    if (NULL == pPrimaryDaughter)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (E_MINUS == std::abs(pPrimaryDaughter->GetParticleId()))
        pNeutrinoPfo->SetParticleId(NU_E);

    if (MU_MINUS == std::abs(pPrimaryDaughter->GetParticleId()))
        pNeutrinoPfo->SetParticleId(NU_MU);

    pNeutrinoPfo->SetCharge(PdgTable::GetParticleCharge(pNeutrinoPfo->GetParticleId()));
    pNeutrinoPfo->SetMass(PdgTable::GetParticleMass(pNeutrinoPfo->GetParticleId()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoBuildingAlgorithm::GetNTwoDHitsInPfo(const ParticleFlowObject *const pPfo) const
{
    unsigned int nTwoDHits(0);
    const ClusterList &clusterList(pPfo->GetClusterList());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType))
            nTwoDHits += pCluster->GetNCaloHits();
    }

    const PfoList &daughterList(pPfo->GetDaughterPfoList());

    for (PfoList::const_iterator iter = daughterList.begin(), iterEnd = daughterList.end(); iter != iterEnd; ++iter)
    {
        nTwoDHits += this->GetNTwoDHitsInPfo(*iter);
    }

    return nTwoDHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DaughterPfoListNames", m_daughterPfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
