/**
 *  @file   LArContent/include/LArClusterSplitting/ClusterSplittingAlgorithm.h
 * 
 *  @brief  Header file for the cluster splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <list>

namespace lar
{

/**
 *  @brief  ClusterSplittingAlgorithm class
 */
class ClusterSplittingAlgorithm : public pandora::Algorithm
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


protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
 
    /**
     *  @brief  Split cluster into its two fragments
     * 
     *  @param  pCluster address of the cluster
     *  @param  splitLayer the layer at which to perform the split
     *  @param  daughters the two cluster fragments
     */
    virtual pandora::StatusCode SplitCluster(pandora::Cluster *const pCluster, const unsigned int splitLayer, std::list<pandora::Cluster*>& daughters);

    /**
     *  @brief  Find layer at which cluster should be split     
     * 
     *  @param  pCluster address of the cluster
     *  @param  splitLayer the layer at which to perform the split
     */
    virtual pandora::StatusCode FindBestSplitLayer(const pandora::Cluster* const pCluster, unsigned int& splitLayer);

   /**
     *  @brief  Whether a cluster is a split candidate
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    virtual bool IsPossibleSplit(const pandora::Cluster *const pCluster) const;


private:

    
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterSplittingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar

#endif // #ifndef LAR_CLUSTER_SPLITTING_ALGORITHM_H
