/**
 *  @file   larpandoracontent/LArTwoDReco/ClusterSplitting/LayerSplittingAlgorithm.h
 *
 *  @brief  Header file for the layer splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LAYER_SPLITTING_ALGORITHM_H
#define LAR_LAYER_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  LayerSplittingAlgorithm class
 */
class LayerSplittingAlgorithm : public ClusterSplittingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    LayerSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode DivideCaloHits(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList,
        pandora::CaloHitList &secondCaloHitList) const;

    /**
     *  @brief Find the best layer for splitting the cluster
     *
     *  @param pCluster the input cluster
     *  @param splitLayer the best layer
     */
    pandora::StatusCode FindBestSplitLayer(const pandora::Cluster *const pCluster, unsigned int &splitLayer) const;

    /**
     *  @brief  Split the cluster into two fragments at the input layer
     *
     *  @param  pCluster the input cluster
     *  @param  splitLayer the split layer
     *  @param  firstCaloHitList the hits in the first cluster fragment
     *  @param  secondCaloHitList the hits in the second cluster fragment
     */
    pandora::StatusCode DivideCaloHits(const pandora::Cluster *const pCluster, const unsigned int &splitLayer,
        pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    /**
     *  @brief  Calculate rms deviation of cluster centroids between two extremal layers
     *
     *  @param  pCluster the input cluster
     *  @param  firstLayer the first extremal layer
     *  @param  secondLayer the second extremal layer
     */
    float CalculateRms(const pandora::Cluster *const pCluster, const unsigned int &firstLayer, const unsigned int& secondLayer) const;

    unsigned int    m_minClusterLayers;       ///<
    unsigned int    m_layerWindow;            ///<
    float           m_maxScatterRms;          ///<
    float           m_maxScatterCosTheta;     ///<
    float           m_maxSlidingCosTheta;     ///<
};

} // namespace lar_content

#endif // #ifndef LAR_LAYER_SPLITTING_ALGORITHM_H
