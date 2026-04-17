/**
 *  @file   larpandoracontent/LArMetrics/EventValidationTool.cc
 *
 *  @brief  Implementation of the event validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"

#include "larpandoracontent/LArMetrics/EventValidationTool.h"

using namespace pandora;

namespace lar_content
{

EventValidationTool::EventValidationTool() :
    m_nuVertexPass1ListName("NeutrinoVertices3D_Pass1"),
    m_nuVertexPass2ListName("NeutrinoVertices3D")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const pMCNu, 
    const LArHierarchyHelper::MCMatchesVector &/*mcMatchesVec*/, const MCParticleVector &targetMC, 
    const PfoVector &/*bestRecoMatch*/)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    EventTreeVars eventTreeVars;
    eventTreeVars.m_run = this->GetPandora().GetRun();
    eventTreeVars.m_subrun = this->GetPandora().GetSubrun();
    eventTreeVars.m_event = this->GetPandora().GetEvent();
    eventTreeVars.m_nTargets = targetMC.size();

    this->GetVertexVariables(pAlgorithm, pMCNu, eventTreeVars);
    this->GetInteractionTypeVariables(pMCNu, eventTreeVars);
    this->FillTree(eventTreeVars);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationTool::GetVertexVariables(const Algorithm *const pAlgorithm, const MCParticle *const pMCNu, EventTreeVars &eventTreeVars)
{
    eventTreeVars.m_trueNuVertex = pMCNu->GetVertex();
    eventTreeVars.m_recoNuVertexPass1 = CartesianVector(-9999.f, -9999.f, -9999.f);
    eventTreeVars.m_recoNuVertexPass2 = CartesianVector(-9999.f, -9999.f, -9999.f);
    eventTreeVars.m_recoNuVertexAccPass1 = -1.f;
    eventTreeVars.m_recoNuVertexAccPass2 = -1.f;

    const VertexList *pNuVertexList_pass1(nullptr), *pNuVertexList_pass2(nullptr);
    PandoraContentApi::GetList(*pAlgorithm, m_nuVertexPass1ListName, pNuVertexList_pass1);
    PandoraContentApi::GetList(*pAlgorithm, m_nuVertexPass2ListName, pNuVertexList_pass2);

    if (pNuVertexList_pass1 && !pNuVertexList_pass1->empty())
    {
        eventTreeVars.m_recoNuVertexPass1 = pNuVertexList_pass1->front()->GetPosition();
        eventTreeVars.m_recoNuVertexAccPass1 = (eventTreeVars.m_trueNuVertex - eventTreeVars.m_recoNuVertexPass1).GetMagnitude();
    }

    if (pNuVertexList_pass2 && !pNuVertexList_pass2->empty())
    {
        eventTreeVars.m_recoNuVertexPass2 = pNuVertexList_pass2->front()->GetPosition();
        eventTreeVars.m_recoNuVertexAccPass2 = (eventTreeVars.m_trueNuVertex - eventTreeVars.m_recoNuVertexPass2).GetMagnitude();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationTool::GetInteractionTypeVariables(const MCParticle *const pMCNu, EventTreeVars &eventTreeVars)
{
    eventTreeVars.m_nuPDG = pMCNu->GetParticleId();
    eventTreeVars.m_nuEnergy = pMCNu->GetEnergy(); 
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle *>(pMCNu));
    eventTreeVars.m_nuVisEnergy = pLArMCParticle->GetVisibleEnergy();

    // CC or NC?
    eventTreeVars.m_isCC = 0;
    if (std::abs(pMCNu->GetParticleId()) == 12)
    {
        for (const MCParticle *const pMCParticle : pMCNu->GetDaughterList())
        {
            if (std::abs(pMCParticle->GetParticleId()) == 11)
                eventTreeVars.m_isCC = 1;
        }
    }
    else if (std::abs(pMCNu->GetParticleId()) == 14)
    {
        for (const MCParticle *const pMCParticle : pMCNu->GetDaughterList())
        {
            if (std::abs(pMCParticle->GetParticleId()) == 13)
                eventTreeVars.m_isCC = 1;
        }
    }
    else if (std::abs(pMCNu->GetParticleId()) == 16)
    {
        for (const MCParticle *const pMCParticle : pMCNu->GetDaughterList())
        {
            if (std::abs(pMCParticle->GetParticleId()) == 15)
                eventTreeVars.m_isCC = 1;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationTool::FillTree(EventTreeVars &eventTreeVars)
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "Run", eventTreeVars.m_run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "Subrun", eventTreeVars.m_subrun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "Event", eventTreeVars.m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "MCEvent_NTargets", eventTreeVars.m_nTargets));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "MCNu_PDG", eventTreeVars.m_nuPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "MCNu_Energy", eventTreeVars.m_nuEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "MCNu_VisEnergy", eventTreeVars.m_nuVisEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "MCInt_IsCC", eventTreeVars.m_isCC));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "MCNu_VertexX", eventTreeVars.m_trueNuVertex.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "MCNu_VertexY", eventTreeVars.m_trueNuVertex.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "MCNu_VertexZ", eventTreeVars.m_trueNuVertex.GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "RecoNu_VertexX_Pass1", eventTreeVars.m_recoNuVertexPass1.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "RecoNu_VertexY_Pass1", eventTreeVars.m_recoNuVertexPass1.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "RecoNu_VertexZ_Pass1", eventTreeVars.m_recoNuVertexPass1.GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "RecoNu_VertexAcc_Pass1", eventTreeVars.m_recoNuVertexAccPass1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "RecoNu_VertexX", eventTreeVars.m_recoNuVertexPass2.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "RecoNu_VertexY", eventTreeVars.m_recoNuVertexPass2.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "RecoNu_VertexZ", eventTreeVars.m_recoNuVertexPass2.GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "EventTree", "RecoNu_VertexAcc_Pass2", eventTreeVars.m_recoNuVertexAccPass2));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "EventTree"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexPass1ListName", m_nuVertexPass1ListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexPass2ListName", m_nuVertexPass2ListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
