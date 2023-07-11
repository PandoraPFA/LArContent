/**
 *  @file   larpandoracontent/LArControlFlow/InteractionSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the list moving algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/InteractionSelectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

InteractionSelectionAlgorithm::InteractionSelectionAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode InteractionSelectionAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    std::set<const MCParticle *> mcSet;
    std::map<const MCParticle *, CaloHitList> mcHitMap;
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
        if (!pMCParticle)
            continue;
        const MCParticle *pParent{pMCParticle};
        while (!pParent->GetParentList().empty())
            pParent = pParent->GetParentList().front();
        if (LArMCParticleHelper::IsNeutrino(pParent))
        {
            mcSet.insert(pParent);
            mcHitMap[pParent].emplace_back(pCaloHit);
        }
    }

    std::cout << "Num hits: " << pCaloHitList->size() << " Num MC: " << mcSet.size() << std::endl;
    for (auto &[pParent, caloHits] : mcHitMap)
        std::cout << "MC: " << pParent << " (" << pParent->GetParticleId() << ") " << caloHits.size() << std::endl;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RenameList<CaloHitList>(*this, "CaloHitList2D", "InputCaloHitList2D"));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RenameList<CaloHitList>(*this, "CaloHitListU", "InputCaloHitU"));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RenameList<CaloHitList>(*this, "CaloHitListV", "InputCaloHitV"));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RenameList<CaloHitList>(*this, "CaloHitListW", "InputCaloHitW"));

    int i{0};
    MCParticleList selectedNeutrinos;
    CaloHitList caloHitListU, caloHitListV, caloHitListW, caloHitList2D;
    for (auto &[pParent, caloHits] : mcHitMap)
    {
        if (std::find(m_interactionIds.begin(), m_interactionIds.end(), i) != m_interactionIds.end())
        {
            selectedNeutrinos.emplace_back(pParent);
            for (const CaloHit *pCaloHit : caloHits)
            {
                switch (pCaloHit->GetHitType())
                {
                    case TPC_VIEW_U:
                        caloHitListU.emplace_back(pCaloHit);
                        caloHitList2D.emplace_back(pCaloHit);
                        break;
                    case TPC_VIEW_V:
                        caloHitListV.emplace_back(pCaloHit);
                        caloHitList2D.emplace_back(pCaloHit);
                        break;
                    case TPC_VIEW_W:
                        caloHitListW.emplace_back(pCaloHit);
                        caloHitList2D.emplace_back(pCaloHit);
                        break;
                    default:
                        break;
                }
            }
        }
        std::cout << "MC: " << pParent << " (" << pParent->GetParticleId() << ") " << caloHits.size() << std::endl;
        ++i;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, selectedNeutrinos, "SelectedNeutrinos"));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList2D, "CaloHitList2D"));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitListU, "CaloHitListU"));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitListV, "CaloHitListV"));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitListW, "CaloHitListW"));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode InteractionSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    (void)xmlHandle;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InteractionIds", m_interactionIds));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
