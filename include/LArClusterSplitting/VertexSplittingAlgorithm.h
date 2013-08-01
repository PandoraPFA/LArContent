/**
 *  @file   LArContent/include/LArClusterSplitting/VertexSplittingAlgorithm.h
 * 
 *  @brief  Header file for the vertex splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SPLITTING_ALGORITHM_H
#define LAR_VERTEX_SPLITTING_ALGORITHM_H 1

#include "LArClusterSplitting/ClusterSplittingAlgorithm.h"

#include "Pandora/Algorithm.h"

#include <list>

namespace lar
{

/**
 *  @brief  VertexSplittingAlgorithm class
 */
class VertexSplittingAlgorithm : public ClusterSplittingAlgorithm
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

   
    /**
     *  @brief  Find layer at which cluster should be split     
     * 
     *  @param  pCluster address of the cluster
     *  @param  splitLayer the layer at which to perform the split
     */
    pandora::StatusCode FindBestSplitLayer(const pandora::Cluster* const pCluster, unsigned int& splitLayer);

     /**
     *  @brief  Whether a cluster is a split candidate
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    bool IsPossibleSplit(const pandora::Cluster *const pCluster) const;
    

    float               m_minSplitDisplacement;
    float               m_minSplitDisplacementSquared;  
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_VERTEX_SPLITTING_ALGORITHM_H
