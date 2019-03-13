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

    float beamMomentum(std::numeric_limits<float>::max()), beamPositionX(std::numeric_limits<float>::max()), beamPositionY(std::numeric_limits<float>::max());
    float beamPositionZ(std::numeric_limits<float>::max()), beamDirectionX(std::numeric_limits<float>::max()), beamDirectionY(std::numeric_limits<float>::max());
    float beamDirectionZ(std::numeric_limits<float>::max()), tof(std::numeric_limits<float>::max());
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
    FloatVector recoDirectionX, recoDirectionY, recoDirectionZ, recoVertexX, recoVertexY, recoVertexZ;
    FloatVector recoDirectionCRX, recoDirectionCRY, recoDirectionCRZ, recoVertexCRX, recoVertexCRY, recoVertexCRZ;
    FloatVector tbStartPointX, tbStartPointY, tbStartPointZ, tbEndPointX, tbEndPointY, tbEndPointZ;
    FloatVector crStartPointX, crStartPointY, crStartPointZ, crEndPointX, crEndPointY, crEndPointZ;
    FloatVector allTestBeamScores;
    int nCosmicRayPfos(0);
    FloatVector cosmicRayX0s;

    for (const Pfo *const pPfo : *pPfoList)
    {
        allTestBeamScores.push_back((pPfo->GetPropertiesMap().count("TestBeamScore")) ? pPfo->GetPropertiesMap().at("TestBeamScore") : -std::numeric_limits<float>::max());

        const Vertex *pVertex(nullptr);

        try
        {
            pVertex = LArPfoHelper::GetVertex(pPfo);
        }
        catch (...)
        {
            std::cout << "ProtoDUNEAnalysisAlgorithm::Run - PFO found without vertex, skipping" << std::endl;
            continue;
        }

        if (LArPfoHelper::IsTestBeam(pPfo))
        {
            float directionX(std::numeric_limits<int>::max()), directionY(std::numeric_limits<int>::max()), directionZ(std::numeric_limits<int>::max());
            float tbStartX(std::numeric_limits<float>::max()), tbStartY(std::numeric_limits<float>::max());
            float tbStartZ(std::numeric_limits<float>::max()), tbEndX(std::numeric_limits<float>::max());
            float tbEndY(std::numeric_limits<float>::max()), tbEndZ(std::numeric_limits<float>::max());

            nBeamPfos++;

            if (std::fabs(pPfo->GetParticleId()) == E_MINUS)
            {
                nShwBeamPfos++;

                try
                {
                    const LArShowerPCA showerPCA(LArPfoHelper::GetPrincipalComponents(pPfo, pVertex));
                    directionX = showerPCA.GetPrimaryAxis().GetX();
                    directionY = showerPCA.GetPrimaryAxis().GetY();
                    directionZ = showerPCA.GetPrimaryAxis().GetZ();

                    const CartesianVector showerStart(showerPCA.GetCentroid() - (showerPCA.GetPrimaryAxis() * 0.5f * showerPCA.GetPrimaryLength()));
                    tbStartX = showerStart.GetX();
                    tbStartY = showerStart.GetY();
                    tbStartZ = showerStart.GetZ();

                    const CartesianVector endStart(showerPCA.GetCentroid() + (showerPCA.GetPrimaryAxis() *0.5f * showerPCA.GetPrimaryLength()));
                    tbEndX = endStart.GetX();
                    tbEndY = endStart.GetY();
                    tbEndZ = endStart.GetZ();
                }
                catch (...)
                {
                    std::cout << "ProtoDUNEAnalysisAlgorithm::Run - Unable to perform PCA fit on shower like test beam PFO" << std::endl;
                }
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

                try
                {
                    LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, slidingFitHalfWindow, layerPitch, trackStateVector);
                }
                catch (...)
                {
                    trackStateVector.clear();
                    std::cout << "ProtoDUNEAnalysisAlgorithm::Run - Unable to perform sliding linear fit on track like test beam PFO" << std::endl;
                }

                // Check front gives start trajectory? Display?
                if (trackStateVector.size() > 0)
                {
                    const LArTrackState &startTrackState(trackStateVector.back());
                    directionX = startTrackState.GetDirection().GetX();
                    directionY = startTrackState.GetDirection().GetY();
                    directionZ = startTrackState.GetDirection().GetZ();
                }
                // Save start and endpoints
                if (trackStateVector.size() > 1)
                {
                    const LArTrackState &startTrackState(trackStateVector.front());
                    tbStartX = startTrackState.GetPosition().GetX();
                    tbStartY = startTrackState.GetPosition().GetY();
                    tbStartZ = startTrackState.GetPosition().GetZ();

                    const LArTrackState &endTrackState(trackStateVector.back());
                    tbEndX = endTrackState.GetPosition().GetX();
                    tbEndY = endTrackState.GetPosition().GetY();
                    tbEndZ = endTrackState.GetPosition().GetZ();
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
            recoVertexX.push_back(pVertex->GetPosition().GetX());
            recoVertexY.push_back(pVertex->GetPosition().GetY());
            recoVertexZ.push_back(pVertex->GetPosition().GetZ());
            tbStartPointX.push_back(tbStartX);
            tbStartPointY.push_back(tbStartY);
            tbStartPointZ.push_back(tbStartZ);
            tbEndPointX.push_back(tbEndX);
            tbEndPointY.push_back(tbEndY);
            tbEndPointZ.push_back(tbEndZ);

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

            float directionX(std::numeric_limits<float>::max()), directionY(std::numeric_limits<float>::max()), directionZ(std::numeric_limits<float>::max());

            // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
            const LArTPC *const pFirstLArTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
            const float layerPitch(pFirstLArTPC->GetWirePitchW());

            // Calculate sliding fit trajectory
            const float slidingFitHalfWindow(20);
            LArTrackStateVector trackStateVector;

            try
            {
                LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, slidingFitHalfWindow, layerPitch, trackStateVector);
            }
            catch (...)
            {
                trackStateVector.clear();
                std::cout << "ProtoDUNEAnalysisAlgorithm::Run - Unable to perform sliding linear fit on track like cosmic PFO" << std::endl;
            }

            // Check front gives start trajectory? Display?
            if (trackStateVector.size() > 0)
            {
                const LArTrackState &startTrackState(trackStateVector.front());
                directionX = startTrackState.GetDirection().GetX();
                directionY = startTrackState.GetDirection().GetY();
                directionZ = startTrackState.GetDirection().GetZ();
            }

            float crStartX(std::numeric_limits<float>::max()), crStartY(std::numeric_limits<float>::max());
            float crStartZ(std::numeric_limits<float>::max()), crEndX(std::numeric_limits<float>::max());
            float crEndY(std::numeric_limits<float>::max()), crEndZ(std::numeric_limits<float>::max());

            if (trackStateVector.size() > 1)
            {
                const LArTrackState &startTrackState(trackStateVector.front());
                crStartX = startTrackState.GetPosition().GetX();
                crStartY = startTrackState.GetPosition().GetY();
                crStartZ = startTrackState.GetPosition().GetZ();

                const LArTrackState &endTrackState(trackStateVector.back());
                crEndX = endTrackState.GetPosition().GetX();
                crEndY = endTrackState.GetPosition().GetY();
                crEndZ = endTrackState.GetPosition().GetZ();
            }

            CaloHitList uHits, vHits, wHits;
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, uHits);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, vHits);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, wHits);
            nHitsCosmicU.push_back(uHits.size());
            nHitsCosmicV.push_back(vHits.size());
            nHitsCosmicW.push_back(wHits.size());
            nHitsCosmicTotal.push_back(uHits.size() + vHits.size() + wHits.size());
            recoDirectionCRX.push_back(directionX);
            recoDirectionCRY.push_back(directionY);
            recoDirectionCRZ.push_back(directionZ);
            recoVertexCRX.push_back(pVertex->GetPosition().GetX());
            recoVertexCRY.push_back(pVertex->GetPosition().GetY());
            recoVertexCRZ.push_back(pVertex->GetPosition().GetZ());
            crStartPointX.push_back(crStartX);
            crStartPointY.push_back(crStartY);
            crStartPointZ.push_back(crStartZ);
            crEndPointX.push_back(crEndX);
            crEndPointY.push_back(crEndY);
            crEndPointZ.push_back(crEndZ);
        }
    }

    // MC
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

    // Beam
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nBeamPfos", nBeamPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nShwBeamPfos", nShwBeamPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrkBeamPfos", nTrkBeamPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "allTestBeamScores", &allTestBeamScores));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsRecoU", &nHitsRecoU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsRecoV", &nHitsRecoV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsRecoW", &nHitsRecoW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsRecoTotal", &nHitsRecoTotal));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoParticleId", &recoParticleId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionX", &recoDirectionX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionY", &recoDirectionY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionZ", &recoDirectionZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexX", &recoVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexY", &recoVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexZ", &recoVertexZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "tbStartPointX", &tbStartPointX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "tbStartPointY", &tbStartPointY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "tbStartPointZ", &tbStartPointZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "tbEndPointX", &tbEndPointX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "tbEndPointY", &tbEndPointY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "tbEndPointZ", &tbEndPointZ));

    // Cosmics
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nCosmicRayPfos", nCosmicRayPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "cosmicRayX0s", &cosmicRayX0s));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsCosmicU", &nHitsCosmicU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsCosmicV", &nHitsCosmicV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsCosmicW", &nHitsCosmicW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsCosmicTotal", &nHitsCosmicTotal));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionCRX", &recoDirectionCRX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionCRY", &recoDirectionCRY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoDirectionCRZ", &recoDirectionCRZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexCRX", &recoVertexCRX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexCRY", &recoVertexCRY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexCRZ", &recoVertexCRZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "crStartPointX", &crStartPointX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "crStartPointY", &crStartPointY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "crStartPointZ", &crStartPointZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "crEndPointX", &crEndPointX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "crEndPointY", &crEndPointY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "crEndPointZ", &crEndPointZ));

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
