/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/TestBeamParticleCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the test beam particle creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/TestBeamParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

TestBeamParticleCreationAlgorithm::TestBeamParticleCreationAlgorithm() :
    m_pfoListName("")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    PfoList neutrinoPfos;

    for (const Pfo *const pPfo : *pPfoList)
    {
        if (!LArPfoHelper::IsNeutrino(pPfo))
            continue;

        const PfoList &daughterList(pPfo->GetDaughterPfoList());

        const Pfo *pPrimaryPfo(nullptr);
        float caloHitMinZ(std::numeric_limits<float>::max());

        for (const Pfo *const pDaughterPfo : daughterList)
        {
            CaloHitList collectedHits;
            LArPfoHelper::GetCaloHits(pDaughterPfo, TPC_3D, collectedHits);

            for (const CaloHit *const pCaloHit : collectedHits)
            {
                if (pCaloHit->GetPositionVector().GetZ() < caloHitMinZ)
                {
                    caloHitMinZ = pCaloHit->GetPositionVector().GetZ();
                    pPrimaryPfo = pDaughterPfo;
                }
            }
        }

        for (const Pfo *const pPrimaryDaughterPfo : daughterList)
        {
            if (pPrimaryDaughterPfo == pPrimaryPfo)
                continue;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPrimaryPfo, pPrimaryDaughterPfo));
        }

        // ATTN: If the primary pfo is shower like, the target beam particle is most likely an electron/positron.  If the primary pfo is track like, the target 
        // beam particle is most likely a pion as pion interactions are more frequent than proton, kaon and muon interactions in the CERN test beam. 
        if (std::abs(pPrimaryPfo->GetParticleId()) != E_MINUS)
        {
            PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;
            pfoMetadata.m_particleId = PI_PLUS;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPrimaryPfo, pfoMetadata));
        }

        neutrinoPfos.push_back(pPfo);
    }

    for (const Pfo *const pNeutrinoPfo : neutrinoPfos)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Pfo>(*this, pNeutrinoPfo));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
