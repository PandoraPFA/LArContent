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

bool CheatedThreeDClusteringTool::Run(const Algorithm *const /*pAlgorithm*/, const pandora::CaloHitList &inputCaloHitList,std::vector<pandora::CaloHitList> &outputCaloHitListsVector)
{
    LArMCParticleHelper::MCContributionMap mcParticleCaloHitListMap;

    for (const CaloHit *const pCaloHit : inputCaloHitList)
    {
        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());

        const MCParticle *const pMainMCParticle(MCParticleHelper::GetMainMCParticle(pParentCaloHit));

        auto it = mcParticleCaloHitListMap.find(pMainMCParticle);

        if (it != mcParticleCaloHitListMap.end()) 
        {
            it->second.push_back(pCaloHit);
        }
        else
        {
            CaloHitList newList;
            newList.push_back(pCaloHit);
            mcParticleCaloHitListMap.insert(std::make_pair(pMainMCParticle, newList));
        } 

    }
    for (auto& pair : mcParticleCaloHitListMap)
        outputCaloHitListsVector.push_back(std::ref(pair.second));

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatedThreeDClusteringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
