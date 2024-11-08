/**
 *  @file   larpandoracontent/LArReclustering/CheatedThreeDClusteringTool.h
 *
 *  @brief  Header file for cheated clustering tool class.
 *
 *  $Log: $
 */

#ifndef LAR_CHEATED_THREE_D_CLUSTERING_TOOL_H
#define LAR_CHEATED_THREE_D_CLUSTERING_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

namespace lar_content
{

class ClusteringTool;

class CheatedThreeDClusteringTool : public ClusteringTool
{
public:

    CheatedThreeDClusteringTool();

private:

    bool Run(const pandora::Algorithm *const pAlgorithm, std::vector<pandora::CaloHitList*> &newCaloHitListsVector) override;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Get the MC index for the true particle that contributes most energy to the calo hit 
     *
     *  @param pAlgorithm the algorithm calling the tool 
     *  @param pCaloHit address of the calo hit
     *
     *  @return MCParticle Index
     */

    int GetMainMcParticleIndex(const pandora::Algorithm *const pAlgorithm, const pandora::CaloHit *const pCaloHit);

    std::string m_mcParticleListName; ///< The mc particle list name 
};

} // namespace lar_content

#endif // #endif LAR_CHEATED_THREE_D_CLUSTERING_TOOL_H
