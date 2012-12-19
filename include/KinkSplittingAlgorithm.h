/**
 *  @file   KinkSplittingAlgorithm.h
 * 
 *  @brief  Header file for the kink splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_KINK_SPLITTING_ALGORITHM_H
#define LAR_KINK_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  KinkSplittingAlgorithm class
 */
class KinkSplittingAlgorithm : public pandora::Algorithm
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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Whether a cluster is a kink candidate
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    bool IsPossibleKink(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Split cluster into its two kink fragments
     * 
     *  @param  pCluster address of the cluster
     *  @param  splitLayer the layer at which to perform the split
     */
    pandora::StatusCode SplitCluster(pandora::Cluster *const pCluster, const unsigned int splitLayer) const;

    int         m_minClusterLayers;             ///< Min number of cluster layers for kink identification
    float       m_minScatteringRms;             ///< Min scattering rms for kink identification
    float       m_maxCosScatteringAngle;        ///< Max cosine of scattering angle for kink identification
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *KinkSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new KinkSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_KINK_SPLITTING_ALGORITHM_H
