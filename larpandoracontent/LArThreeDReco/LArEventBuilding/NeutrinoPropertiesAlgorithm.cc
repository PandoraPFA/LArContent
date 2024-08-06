/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoPropertiesAlgorithm.cc
 *
 *  @brief  Implementation of the neutrino properties algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoPropertiesAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NeutrinoPropertiesAlgorithm::NeutrinoPropertiesAlgorithm() :
    m_includeIsolatedHits(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoPropertiesAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoPropertiesAlgorithm: unable to find pfo list " << m_neutrinoPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // ATTN Enforces that only one pfo, of neutrino-type, be in the specified input list
    const ParticleFlowObject *const pNeutrinoPfo((1 == pPfoList->size()) ? *(pPfoList->begin()) : nullptr);

    if (!pNeutrinoPfo || !LArPfoHelper::IsNeutrino(pNeutrinoPfo))
        return STATUS_CODE_FAILURE;

    // ATTN At this (maybe unconventional) stage, remove any lone, placeholder neutrinos, with no daughter particles
    if (pNeutrinoPfo->GetDaughterPfoList().empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pNeutrinoPfo, m_neutrinoPfoListName));
    }
    else
    {
        this->SetNeutrinoId(pNeutrinoPfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoPropertiesAlgorithm::SetNeutrinoId(const ParticleFlowObject *const pNeutrinoPfo) const
{
    unsigned int nPrimaryTwoDHits(0);
    const ParticleFlowObject *pPrimaryDaughter(nullptr);

    PfoVector daughterPfoVector(pNeutrinoPfo->GetDaughterPfoList().begin(), pNeutrinoPfo->GetDaughterPfoList().end());
    std::sort(daughterPfoVector.begin(), daughterPfoVector.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfoVector)
    {
        const unsigned int nTwoDHits(this->GetNTwoDHitsInPfoChain(pDaughterPfo));

        if (!pPrimaryDaughter || (nTwoDHits > nPrimaryTwoDHits))
        {
            nPrimaryTwoDHits = nTwoDHits;
            pPrimaryDaughter = pDaughterPfo;
        }
    }

    if (!pPrimaryDaughter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PandoraContentApi::ParticleFlowObject::Metadata metadata;

    if (E_MINUS == std::abs(pPrimaryDaughter->GetParticleId()))
    {
        metadata.m_particleId = NU_E;
    }
    else if (MU_MINUS == std::abs(pPrimaryDaughter->GetParticleId()))
    {
        metadata.m_particleId = NU_MU;
    }

    if (metadata.m_particleId.IsInitialized())
    {
        metadata.m_charge = PdgTable::GetParticleCharge(metadata.m_particleId.Get());
        metadata.m_mass = PdgTable::GetParticleMass(metadata.m_particleId.Get());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pNeutrinoPfo, metadata));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoPropertiesAlgorithm::GetNTwoDHitsInPfoChain(const ParticleFlowObject *const pPfo) const
{
    unsigned int nTwoDHits(0);

    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType))
        {
            nTwoDHits += pCluster->GetNCaloHits();

            if (m_includeIsolatedHits)
                nTwoDHits += pCluster->GetNIsolatedCaloHits();
        }
    }

    for (const ParticleFlowObject *const pDaughterPfo : pPfo->GetDaughterPfoList())
    {
        nTwoDHits += this->GetNTwoDHitsInPfoChain(pDaughterPfo);
    }

    return nTwoDHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoPropertiesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IncludeIsolatedHits", m_includeIsolatedHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
