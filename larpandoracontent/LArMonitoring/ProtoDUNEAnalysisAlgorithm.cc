/**
 *  @file   larpandoracontent/LArMonitoring/ProtoDUNEAnalysisAlgorithm.cc
 *
 *  @brief  Implementation of the ProtoDUNE data analysis algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/ProtoDUNEAnalysisAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ProtoDUNEAnalysisAlgorithm::ProtoDUNEAnalysisAlgorithm() :
    m_visualDisplay(false),
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ProtoDUNEAnalysisAlgorithm::~ProtoDUNEAnalysisAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoDUNEAnalysisAlgorithm::Run()
{
    m_eventNumber++;

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const PfoList *pPfoList(nullptr);
    (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);

    if (m_visualDisplay)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XY, -1.f, 1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList, "CurrentPfos", AUTO, true, true));
    }

    // Triggered Information
    int isTriggered(pMCParticleList->empty() ? 0 : 1);

    float beamMomentum(std::numeric_limits<float>::max()), beamPositionX(std::numeric_limits<float>::max()), beamPositionY(std::numeric_limits<float>::max()), beamPositionZ(std::numeric_limits<float>::max()),
          beamDirectionX(std::numeric_limits<float>::max()), beamDirectionY(std::numeric_limits<float>::max()), beamDirectionZ(std::numeric_limits<float>::max()), tof(std::numeric_limits<float>::max());
    int ckov0Status(std::numeric_limits<int>::max()), ckov1Status(std::numeric_limits<int>::max());

    if (isTriggered && pMCParticleList->size() == 1)
    {
        // Data Event
        const MCParticle *pMCParticle(pMCParticleList->front());
        beamMomentum = pMCParticle->GetEnergy();
        beamDirectionX = pMCParticle->GetMomentum().GetX();
        beamDirectionY = pMCParticle->GetMomentum().GetY();
        beamDirectionZ = pMCParticle->GetMomentum().GetZ();
        beamPositionX = pMCParticle->GetVertex().GetX();
        beamPositionY = pMCParticle->GetVertex().GetY();
        beamPositionZ = pMCParticle->GetVertex().GetZ();
        tof = pMCParticle->GetEndpoint().GetX();
        ckov0Status = static_cast<int>(pMCParticle->GetEndpoint().GetY());
        ckov1Status = static_cast<int>(pMCParticle->GetEndpoint().GetZ());
    }
    else if (isTriggered && pMCParticleList->size() >  1)
    {
        // MC Event
        const MCParticle *pMCTriggeredParticle(nullptr);

        for (const MCParticle *pMCParticle : *pMCParticleList)
        {
            const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));

            if (LArMCParticleHelper::IsPrimary(pMCParticle) && nuance == 2001)
            {
                if (pMCTriggeredParticle)
                {
                    std::cout << "More than one target triggered mc particle" << std::endl;
                    pMCTriggeredParticle = nullptr;
                    break;
                }
                pMCTriggeredParticle = pMCParticle;
            }
        }

        if (pMCTriggeredParticle)
        {
            beamMomentum = pMCTriggeredParticle->GetMomentum().GetMagnitude();
            beamDirectionX = pMCTriggeredParticle->GetMomentum().GetX();
            beamDirectionY = pMCTriggeredParticle->GetMomentum().GetY();
            beamDirectionZ = pMCTriggeredParticle->GetMomentum().GetZ();
            beamPositionX = pMCTriggeredParticle->GetVertex().GetX();
            beamPositionY = pMCTriggeredParticle->GetVertex().GetY();
            beamPositionZ = pMCTriggeredParticle->GetVertex().GetZ();
            tof = std::numeric_limits<float>::max();
            ckov0Status = pMCTriggeredParticle->GetParticleId();
            ckov1Status = std::numeric_limits<int>::max();
        }
    }

    // Reconstruction Information
    int nBeamPfos(0), nTrkBeamPfos(0), nShwBeamPfos(0);
    IntVector nHitsRecoU, nHitsRecoV, nHitsRecoW, nHitsRecoTotal, recoParticleId;
    IntVector nHitsCosmicU, nHitsCosmicV, nHitsCosmicW, nHitsCosmicTotal;
    FloatVector recoDirectionX, recoDirectionY, recoDirectionZ;;
    int nCosmicRayPfos(0);
    FloatVector cosmicRayX0s;

    for (const Pfo *const pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsTestBeam(pPfo))
        {
            const Vertex *const pVertex = LArPfoHelper::GetVertex(pPfo);
            float directionX(std::numeric_limits<int>::max()), directionY(std::numeric_limits<int>::max()), directionZ(std::numeric_limits<int>::max());
            nBeamPfos++;

            if (std::fabs(pPfo->GetParticleId()) == E_MINUS)
            {
                nShwBeamPfos++;

                const LArShowerPCA showerPCA(LArPfoHelper::GetPrincipalComponents(pPfo, pVertex));
                directionX = showerPCA.GetPrimaryAxis().GetX();
                directionY = showerPCA.GetPrimaryAxis().GetY();
                directionZ = showerPCA.GetPrimaryAxis().GetZ();
            }
            else
            {
                nTrkBeamPfos++;

                // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
                const LArTPC *const pFirstLArTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
                const float layerPitch(pFirstLArTPC->GetWirePitchW());

                // Calculate sliding fit trajectory
                const float slidingFitHalfWindow(20);
                LArTrackStateVector trackStateVector;
                LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, slidingFitHalfWindow, layerPitch, trackStateVector);
                // Check front gives start trajectory? Display?
                if (trackStateVector.size() > 0)
                {
                    const LArTrackState &startTrackState(trackStateVector.back());
                    directionX = startTrackState.GetDirection().GetX();
                    directionY = startTrackState.GetDirection().GetY();
                    directionZ = startTrackState.GetDirection().GetZ();
                }
            }

            CaloHitList uHits, vHits, wHits;
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, uHits);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, vHits);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, wHits);
            nHitsRecoU.push_back(uHits.size());
            nHitsRecoV.push_back(vHits.size());
            nHitsRecoW.push_back(wHits.size());
            nHitsRecoTotal.push_back(uHits.size() + vHits.size() + wHits.size());
            recoParticleId.push_back(pPfo->GetParticleId());
            recoDirectionX.push_back(directionX);
            recoDirectionY.push_back(directionY);
            recoDirectionZ.push_back(directionZ);

            if (m_visualDisplay)
            {
                CartesianVector pfoPosition(pVertex->GetPosition());
                CartesianVector pfoDirection(directionX, directionY, directionZ);
                CartesianVector pfoEndpoint(pfoPosition + pfoDirection * 100.f);
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pfoPosition, &pfoEndpoint, "Beam_PfoInfo_" + std::to_string(nBeamPfos), MAGENTA, 3, 2));

                PfoList testBeamPfo;
                testBeamPfo.push_back(pPfo);
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &testBeamPfo, "TestBeamPfo", MAGENTA, true, true));
            }
        }
        else if (pPfo->GetParticleId() == 13)
        {
            // Cosmic Ray Pfo
            nCosmicRayPfos++;
            const Pfo *const pParentPfo(LArPfoHelper::GetParentPfo(pPfo));

            if (pParentPfo->GetPropertiesMap().count("X0"))
                cosmicRayX0s.push_back(pParentPfo->GetPropertiesMap().at("X0"));

            CaloHitList uHits, vHits, wHits;
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, uHits);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, vHits);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, wHits);
            nHitsCosmicU.push_back(uHits.size());
            nHitsCosmicV.push_back(vHits.size());
            nHitsCosmicW.push_back(wHits.size());
            nHitsCosmicTotal.push_back(uHits.size() + vHits.size() + wHits.size());
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isTriggered", isTriggered));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beamMomentum", beamMomentum));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beamPositionX", beamPositionX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beamPositionY", beamPositionY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beamPositionZ", beamPositionZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beamDirectionX", beamDirectionX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beamDirectionY", beamDirectionY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beamDirectionZ", beamDirectionZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "tof", tof));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ckov0Status", ckov0Status));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ckov1Status", ckov1Status));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nBeamPfos", nBeamPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nShwBeamPfos", nShwBeamPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrkBeamPfos", nTrkBeamPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsRecoU", &nHitsRecoU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsRecoV", &nHitsRecoV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsRecoW", &nHitsRecoW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsRecoTotal", &nHitsRecoTotal));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoParticleId", &recoParticleId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionX", &recoDirectionX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionY", &recoDirectionY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionZ", &recoDirectionZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nCosmicRayPfos", nCosmicRayPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "cosmicRayX0s", &cosmicRayX0s));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsCosmicU", &nHitsCosmicU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsCosmicV", &nHitsCosmicV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsCosmicW", &nHitsCosmicW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsCosmicTotal", &nHitsCosmicTotal));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    if (m_visualDisplay)
    {
        if (beamPositionX < 10000.f)
        {
            CartesianVector beamPosition(beamPositionX, beamPositionY, beamPositionZ);
            CartesianVector beamDirection(beamDirectionX, beamDirectionY, beamDirectionZ);
            CartesianVector endpoint(beamPosition + beamDirection * 100.f);
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &beamPosition, &endpoint, "Beam_TriggerInfo", RED, 3, 1));
        }

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoDUNEAnalysisAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VisualDisplay", m_visualDisplay));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
