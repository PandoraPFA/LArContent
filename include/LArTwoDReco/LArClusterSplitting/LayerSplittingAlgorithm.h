/**
 *  @file   LArContent/include/LArTwoDReco/ClusterSplitting/LayerSplittingAlgorithm.h
 *
 *  @brief  Header file for the layer splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LAYER_SPLITTING_ALGORITHM_H
#define LAR_LAYER_SPLITTING_ALGORITHM_H 1

#include "LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar
{

/**
 *  @brief  LayerSplittingAlgorithm class
 */
class LayerSplittingAlgorithm : public ClusterSplittingAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode SplitCluster(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList,
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
    pandora::StatusCode SplitCluster(const pandora::Cluster *const pCluster, const unsigned int &splitLayer,
        pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    /**
     *  @brief  Calculate rms deviation of cluster centroids between two extremal layers
     *
     *  @param  pCluster the input cluster
     *  @param  firstLayer the first extremal layer
     *  @param  secondLayer the second extremal layer
     */
    float CalculateRms(const pandora::Cluster *const pCluster, const unsigned int &firstLayer, const unsigned int& secondLayer) const;

    unsigned int    m_layerWindow;            ///<
    unsigned int    m_minClusterLayers;       ///<
    float           m_maxScatterRms;          ///<
    float           m_maxScatterCosTheta;     ///<
    float           m_maxSlidingCosTheta;     ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *LayerSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new LayerSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_LAYER_SPLITTING_ALGORITHM_H
