/**
 *  @file   larpandoracontent/LArMonitoring/ProtoDUNEAnalysisAlgorithm.cc
 *
 *  @brief  Implementation of the ProtoDUNE data analysis algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/ProtoDUNEAnalysisAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ProtoDUNEAnalysisAlgorithm::ProtoDUNEAnalysisAlgorithm() :
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

    int isTriggered(pMCParticleList->empty() ? 0 : 1);

    int nBeamPfos(0);

    for (const Pfo *const pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsTestBeam(pPfo))
            nBeamPfos++;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isTriggered", isTriggered));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nBeamPfos", nBeamPfos));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoDUNEAnalysisAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
