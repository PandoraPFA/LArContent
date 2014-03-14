/**
 *  @file   LArContent/include/LArTwoDReco/LArSeedFinding/VertexSeedFindingAlgorithm.h
 * 
 *  @brief  Header file for the vertex seed finding algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SEED_FINDING_ALGORITHM_H
#define LAR_VERTEX_SEED_FINDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "SeedFindingBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  VertexSeedFindingAlgorithm class
 */
class VertexSeedFindingAlgorithm : public SeedFindingBaseAlgorithm
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
    void GetSeedClusterList(const pandora::ClusterVector &candidateClusters, pandora::ClusterList &seedClusterList) const;

    unsigned int        m_minClusterHitsNode;           ///< The min number of cluster hits for a node 
    unsigned int        m_minClusterHitsEmission;       ///< The min number of cluster hits for an emission
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexSeedFindingAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexSeedFindingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_VERTEX_SEED_FINDING_ALGORITHM_H
