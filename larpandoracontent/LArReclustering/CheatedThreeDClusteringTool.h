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

    std::vector<std::reference_wrapper<pandora::CaloHitList>> Run(const pandora::Algorithm *const pAlgorithm, std::reference_wrapper<pandora::CaloHitList> &inputCaloHitList);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_mcParticleListName; ///< The mc particle list name 
};

} // namespace lar_content

#endif // #endif LAR_CHEATED_THREE_D_CLUSTERING_TOOL_H