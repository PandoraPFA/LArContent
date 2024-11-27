/**
 *  @file   larpandoracontent/LArReclustering/SimplePCAThreeDClusteringTool.h
 *
 *  @brief  Header file for the simple PCA-based clustering tool class.
 *
 *  $Log: $
 */

#ifndef LAR_SIMPLE_PCA_THREE_D_CLUSTERING_TOOL_H
#define LAR_SIMPLE_PCA_THREE_D_CLUSTERING_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

namespace lar_content
{

class SimplePCAThreeDClusteringTool : public ClusteringTool
{
public:

    SimplePCAThreeDClusteringTool();

private:

    std::vector<std::reference_wrapper<pandora::CaloHitList>> Run(const pandora::Algorithm *const pAlgorithm, std::reference_wrapper<pandora::CaloHitList> &inputCaloHitList);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/);
};

} // namespace lar_content

#endif // #endif LAR_SIMPLE_PCA_THREE_D_CLUSTERING_TOOL_H
