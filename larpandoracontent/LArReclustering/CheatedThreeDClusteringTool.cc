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

bool CheatedThreeDClusteringTool::Run(const Algorithm *const /*pAlgorithm*/, std::vector<pandora::CaloHitList*> &newCaloHitListsVector)
{
    if (newCaloHitListsVector.size() != 1)
        return false;

    CaloHitList *initialCaloHitList = newCaloHitListsVector.at(0);
    newCaloHitListsVector.clear();

    std::map<const pandora::MCParticle *, CaloHitList*> McParticleCaloHitListMap;

    for (const CaloHit *const pCaloHit : *initialCaloHitList)
    {
        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());

        const MCParticle *const pMainMCParticle(MCParticleHelper::GetMainMCParticle(pParentCaloHit));

        auto it = McParticleCaloHitListMap.find(pMainMCParticle);

        if (it != McParticleCaloHitListMap.end()) 
        {
            it->second->push_back(pCaloHit);
        }
        else
        {
            CaloHitList* newList = new CaloHitList();
            newList->push_back(pCaloHit);
            McParticleCaloHitListMap.insert(std::make_pair(pMainMCParticle, newList));
        } 

    }
    for (const auto& pair : McParticleCaloHitListMap)
        newCaloHitListsVector.push_back(pair.second);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatedThreeDClusteringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
