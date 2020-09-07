/**
 *  @file   larpandoracontent/LArMonitoring/HitCheckingAlgorithm.cc
 *
 *  @brief  Implementation of the hit checking algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/HitCheckingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

HitCheckingAlgorithm::HitCheckingAlgorithm() :
    m_writeToTree(false),    
    m_treeName("HitCountTree"),
    m_fileName("HitCount.root"),
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitCheckingAlgorithm::~HitCheckingAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "EventValidationBaseAlgorithm: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode HitCheckingAlgorithm::Run()
{
    ++m_eventNumber;
    
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (pCaloHitList)
    {
        int hitNumber(pCaloHitList->size());

        std::cout << "Number of hits in TPC volume: " << hitNumber << std::endl;

        if (m_writeToTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "hitNumber", hitNumber));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }
    }

    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            //std::cout << "std::fabs(pHitParticle->GetParticleId()): " << std::fabs(pHitParticle->GetParticleId()) << std::endl;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            if (std::fabs(pHitParticle->GetParticleId()) == 11)
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "ELECTRON", RED, 2);

            if (std::fabs(pHitParticle->GetParticleId()) == 13)
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "MUON", BLUE, 2);
        }
        catch (const StatusCodeException&)
        {
            continue;
        }

    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    StatusCode HitCheckingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treeName));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_fileName));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
