/**
 *  @file   LArContent/include/LArTwoDSeed/VertexSeedFindingAlgorithm.h
 * 
 *  @brief  Header file for the vertex seed finding algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SEED_FINDING_ALGORITHM_H
#define LAR_VERTEX_SEED_FINDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  VertexSeedFindingAlgorithm class
 */
class VertexSeedFindingAlgorithm : public pandora::Algorithm
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
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     * 
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Merge any vertex seed candidates that are unambiguously associated
     * 
     *  @param  eventVertex the event vertex position vector
     *  @param  vertexSeedClusterList the list of vertex seed candidates
     */
    void MakeVertexSeedMerges(const pandora::CartesianVector &eventVertex, pandora::ClusterList &vertexSeedClusterList) const;

    std::string         m_seedClusterListName;          ///< The seed cluster list name
    std::string         m_nonSeedClusterListName;       ///< The non seed cluster list name

    unsigned int        m_minClusterLayers;             ///< The min number of layers for a clean cluster
    float               m_minClusterLengthSquared;      ///< The min length (squared) for a clean cluster
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexSeedFindingAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexSeedFindingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_VERTEX_SEED_FINDING_ALGORITHM_H
