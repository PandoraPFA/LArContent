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

StatusCode NeutrinoPropertiesAlgorithm::Run()
{
    const PfoList *pPfoList(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoPropertiesAlgorithm: unable to find pfo list " << m_neutrinoPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // ATTN Enforces that only one pfo, of neutrino-type, be in the specified input list
    const ParticleFlowObject *const pNeutrinoPfo((1 == pPfoList->size()) ? *(pPfoList->begin()) : NULL);

    if (!pNeutrinoPfo || !LArPfoHelper::IsNeutrino(pNeutrinoPfo))
        return STATUS_CODE_FAILURE;

    this->SetNeutrinoId(pNeutrinoPfo);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoPropertiesAlgorithm::SetNeutrinoId(const ParticleFlowObject *const pNeutrinoPfo) const
{
    if (pNeutrinoPfo->GetDaughterPfoList().empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    unsigned int nPrimaryTwoDHits(0);
    const ParticleFlowObject *pPrimaryDaughter(NULL);

    for (PfoList::const_iterator dIter = pNeutrinoPfo->GetDaughterPfoList().begin(), dIterEnd = pNeutrinoPfo->GetDaughterPfoList().end();
        dIter != dIterEnd; ++dIter)
    {
        const ParticleFlowObject *const pDaughterPfo(*dIter);
        const unsigned int nTwoDHits(this->GetNTwoDHitsInPfo(pDaughterPfo));

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
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AlterMetadata(*this, pNeutrinoPfo, metadata));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoPropertiesAlgorithm::GetNTwoDHitsInPfo(const ParticleFlowObject *const pPfo) const
{
    unsigned int nTwoDHits(0);

    const ClusterList &pfoClusterList(pPfo->GetClusterList());

    for (ClusterList::const_iterator iter = pfoClusterList.begin(), iterEnd = pfoClusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster(*iter);

        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            continue;

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

StatusCode NeutrinoPropertiesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
