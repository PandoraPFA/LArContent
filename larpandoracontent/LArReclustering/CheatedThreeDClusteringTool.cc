/**
 *  @file   larpandoracontent/LArReclustering/CheatedThreeDClusteringTool.cc
 *
 *  @brief  Implementation file for the reclustering algorithm that uses transverse calorimetric profiles.
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
    for (const CaloHit *const pCaloHit : *initialCaloHitList)
    {
        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());

        int mainMcParticleIndex = this->GetMainMcParticleIndex(pAlgorithm, pParentCaloHit);

        // Flag to indicate if a matching list was found
        bool found = false;

        // Search for an existing list with matching MC particle index
        for (CaloHitList *caloHitList : newCaloHitListsVector)
        {
            const CaloHit *const pListParentCaloHit = static_cast<const CaloHit *>(caloHitList->front()->GetParentAddress());
            int listMainMcParticleIndex = this->GetMainMcParticleIndex(pAlgorithm, pListParentCaloHit);
            if (listMainMcParticleIndex == mainMcParticleIndex)
            {
                found = true;
                caloHitList->push_back(pCaloHit);
                break;
            }
        }

        // If a list with matching MC particle hits is not found, create one
        if (!found)
        {
            CaloHitList* newList = new CaloHitList();
            newList->push_back(pCaloHit);
            newCaloHitListsVector.push_back(newList);
        }
    }

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
        std::cerr << "MCParticleList is null!" << std::endl;
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
