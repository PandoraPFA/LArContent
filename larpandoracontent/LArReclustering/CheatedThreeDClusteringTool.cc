/**
 *  @file   larpandoracontent/LArReclustering/CheatedThreeDClusteringTool.cc
 *
 *  @brief  Implementation file for the cheated threeD reclustering algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArReclustering/CheatedThreeDClusteringTool.h"


using namespace pandora;

namespace lar_content
{

CheatedThreeDClusteringTool::CheatedThreeDClusteringTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatedThreeDClusteringTool::Run(const Algorithm *const pAlgorithm, std::vector<pandora::CaloHitList*> &newCaloHitListsVector)
{
    if (newCaloHitListsVector.empty() || newCaloHitListsVector.size() != 1)
        return false;

    CaloHitList *initialCaloHitList = newCaloHitListsVector.at(0);
    newCaloHitListsVector.clear();

    std::map<int, CaloHitList*> McIdCaloHitListMap;

    for (const CaloHit *const pCaloHit : *initialCaloHitList)
    {
        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());

        int mainMcParticleIndex = this->GetMainMcParticleIndex(pAlgorithm, pParentCaloHit);

        if (McIdCaloHitListMap.find(mainMcParticleIndex) != McIdCaloHitListMap.end()) 
        {
            McIdCaloHitListMap.at(mainMcParticleIndex)->push_back(pCaloHit);
        }
        else
        {
            CaloHitList* newList = new CaloHitList();
            McIdCaloHitListMap.insert(std::make_pair(mainMcParticleIndex, newList));
        } 

    }
    for (const auto& pair : McIdCaloHitListMap) newCaloHitListsVector.push_back(pair.second);
    return true;
}


//------------------------------------------------------------------------------------------------------------------------------------------

int CheatedThreeDClusteringTool::GetMainMcParticleIndex(const Algorithm *const pAlgorithm, const pandora::CaloHit *const pCaloHit)
{
    if (!pCaloHit) {
        std::cerr << "pCaloHit is null!" << std::endl;
        return -1;
    }
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, m_mcParticleListName, pMCParticleList));
    if (!pMCParticleList) {
        std::cerr << "ERROR in CheatedThreeDClusteringTool: MCParticleList is null!" << std::endl;
        return -1;
    }
    MCParticleVector mcParticleVector(pMCParticleList->begin(),pMCParticleList->end());
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    int iMcPart(0); 
    for (const auto &weightMapEntry : pCaloHit->GetMCParticleWeightMap()) 
    { 
      if(weightMapEntry.second>0.5) 
      { 
        iMcPart=0;  
        for(const MCParticle *const pMCParticle: mcParticleVector) 
        { 
            if(pMCParticle==weightMapEntry.first) { break;} 
            iMcPart++; 
        }
      } 
    }
    return iMcPart;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatedThreeDClusteringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
