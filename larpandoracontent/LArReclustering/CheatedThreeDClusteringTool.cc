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

bool CheatedThreeDClusteringTool::Run(const pandora::CaloHitList &inputCaloHitList, std::vector<pandora::CaloHitList> &outputCaloHitListsVector)
{
    LArMCParticleHelper::MCContributionMap mcParticleCaloHitListMap;

    for (const CaloHit *const pCaloHit3D : inputCaloHitList)
    {
        const CaloHit *const pParentCaloHit2D{static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress())};

        const MCParticle *const pMainMCParticle(MCParticleHelper::GetMainMCParticle(pParentCaloHit2D));

        auto it = mcParticleCaloHitListMap.find(pMainMCParticle);

        if (it != mcParticleCaloHitListMap.end())
        {
            it->second.push_back(pCaloHit3D);
        }
        else
        {
            CaloHitList newList;
            newList.push_back(pCaloHit3D);
            mcParticleCaloHitListMap.insert(std::make_pair(pMainMCParticle, newList));
        }
    }
    for (auto &pair : mcParticleCaloHitListMap)
        outputCaloHitListsVector.push_back(std::ref(pair.second));

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatedThreeDClusteringTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
