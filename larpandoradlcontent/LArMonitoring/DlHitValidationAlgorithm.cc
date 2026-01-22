/**
 *  @file   larpandoradlcontent/LArMonitoring/DlHitValidationAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArMonitoring/DlHitValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlHitValidationAlgorithm::DlHitValidationAlgorithm() :
    m_confusionU(),
    m_confusionV(),
    m_confusionW()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlHitValidationAlgorithm::~DlHitValidationAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "u_true_shower", m_confusionU[0][0]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "u_false_shower", m_confusionU[1][0]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "u_false_track", m_confusionU[0][1]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "u_true_track", m_confusionU[1][1]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "v_true_shower", m_confusionV[0][0]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "v_false_shower", m_confusionV[1][0]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "v_false_track", m_confusionV[0][1]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "v_true_track", m_confusionV[1][1]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "w_true_shower", m_confusionW[0][0]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "w_false_shower", m_confusionW[1][0]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "w_false_track", m_confusionW[0][1]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "w_true_track", m_confusionW[1][1]));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "confusion_tree"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "confusion_tree", "confusion.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "DlHitValidationAlgorithm: Unable to write confusion_tree to file" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitValidationAlgorithm::Run()
{
    const int SHOWER_IDX{0}, TRACK_IDX{1};
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

        const HitType view{pCaloHitList->front()->GetHitType()};

        if (!(view == TPC_VIEW_U || view == TPC_VIEW_V || view == TPC_VIEW_W))
            return STATUS_CODE_NOT_ALLOWED;

        LArMCParticleHelper::PrimaryParameters parameters;
        // Only care about reconstructability with respect to the current view, so skip good view check
        parameters.m_minHitsForGoodView = 0;
        // Turn off max photo propagation for now, only care about killing off daughters of neutrons
        parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                const int pdg{std::abs(pMCParticle->GetParticleId())};
                const int truth{(pdg == 11 || pdg == 22) ? SHOWER_IDX : TRACK_IDX};
                const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
                if (pLArCaloHit)
                {
                    const float pTrack{pLArCaloHit->GetTrackProbability()};
                    const float pShower{pLArCaloHit->GetShowerProbability()};
                    const int cls{(pShower > pTrack) ? SHOWER_IDX : TRACK_IDX};
                    if (view == TPC_VIEW_U)
                        ++m_confusionU[truth][cls];
                    else if (view == TPC_VIEW_V)
                        ++m_confusionV[truth][cls];
                    else
                        ++m_confusionW[truth][cls];
                }
            }
            catch (const StatusCodeException &)
            {
                continue;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
