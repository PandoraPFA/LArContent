#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArMonitoring/PullDataAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content 
	{
  	PullDataAlgorithm::PullDataAlgorithm() :
    m_eventNumber(4500),
    m_treeName(),
    m_fileName()
  	{
  	}

  PullDataAlgorithm::~PullDataAlgorithm() {
    try {
      PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
    }
    catch (const StatusCodeException &) {
      std::cout << "UNABLE TO WRITE TREE" << std::endl;
    }
  }

  StatusCode PullDataAlgorithm::Run() {

    m_eventNumber++;
    std::cout << m_eventNumber << std::endl;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    
    for(const MCParticle *const pMCParticle : *pMCParticleList) {

      if(pMCParticle->GetParticleId() != 13)
	continue;

      float theta0XZ;
      float theta0YZ;

      GetLArSoftAngles(pMCParticle->GetMomentum(), theta0XZ, theta0YZ);

      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "EventNumber", m_eventNumber));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "ParticleID", pMCParticle->GetParticleId()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Momentum", pMCParticle->GetMomentum().GetMagnitude()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Hierarchy", LArMCParticleHelper::GetHierarchyTier(pMCParticle)));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Energy", pMCParticle->GetEnergy()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "X", pMCParticle->GetVertex().GetX()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Y", pMCParticle->GetVertex().GetY()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Z", pMCParticle->GetVertex().GetZ()));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Theta0XZ", theta0XZ));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Theta0YZ", theta0YZ));

      PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));

    }

   return STATUS_CODE_SUCCESS;

  }

  void PullDataAlgorithm::GetLArSoftAngles(const CartesianVector &vector, float &theta0XZ, float &theta0YZ) {

    theta0YZ = asin(vector.GetY()/vector.GetMagnitude());
    theta0XZ = atan2(vector.GetX(), vector.GetZ());

  }


 StatusCode PullDataAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {

   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));

   return STATUS_CODE_SUCCESS;

 }

} //namespace lar_content

