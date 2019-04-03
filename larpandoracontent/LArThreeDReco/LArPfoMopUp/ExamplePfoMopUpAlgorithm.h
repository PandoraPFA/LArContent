/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/ExamplePfoMopUpAlgorithm.h
 * 
 *  @brief  Header file for the example pfo mop up algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_EXAMPLE_PFO_MOP_UP_ALGORITHM_H
#define LAR_EXAMPLE_PFO_MOP_UP_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ExamplePfoMopUpAlgorithm class
 */
class ExamplePfoMopUpAlgorithm : public PfoMopUpBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ExamplePfoMopUpAlgorithm();

    /**
     *  @brief  ClusterMerge class
     */
    class ClusterMerge
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pDaughterCluster the address of the candidate daughter cluster
         *  @param  distance the distance associated with the pfo merge
         */
        ClusterMerge(const pandora::Cluster *const pDaughterCluster, const float distance);

        /**
         *  @brief  Get the address of the candidate daughter cluster
         *
         *  @return the address of the candidate daughter cluster
         */
        const pandora::Cluster *GetDaughterCluster() const;

        /**
         *  @brief  Get the distance associated with the pfo merge
         *
         *  @return the distance associated with the pfo merge
         */
        float GetDistance() const;

        /**
         *  @brief  operator <
         * 
         *  @param  rhs object for comparison
         * 
         *  @return boolean
         */
        bool operator<(const ClusterMerge &rhs) const;

    private:
        const pandora::Cluster *m_pDaughterCluster;     ///< The address of the candidate daughter cluster
        float                   m_distance;             ///< The distance associated with the merge
    };

    typedef std::vector<ClusterMerge> ClusterMergeList;
    typedef std::unordered_map<const pandora::Cluster*, ClusterMergeList> ClusterMergeMap;
    typedef std::unordered_map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    typedef std::unordered_map<const pandora::Cluster*, const pandora::Cluster*> ClusterReplacementMap;

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Get all 3d clusters contained in the input pfo lists and a mapping from clusters to pfos
     * 
     *  @param  clusterList3D to receive the sorted list of 3d clusters
     *  @param  clusterToPfoMap to receive the mapping from 3d cluster to pfo
     */
    void GetThreeDClusters(pandora::ClusterList &clusterList3D, ClusterToPfoMap &clusterToPfoMap) const;

    /**
     *  @brief  Get the cluster merge map describing all potential 3d cluster merges
     * 
     *  @param  clusterList3D the sorted list of 3d clusters
     *  @param  clusterMergeMap to receive the populated cluster merge map
     */
    void GetClusterMergeMap(const pandora::ClusterList &clusterList3D, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Make pfo merges based on the provided cluster merge map
     * 
     *  @param  clusterList3D the sorted list of 3d clusters
     *  @param  clusterToPfoMap the mapping from 3d cluster to pfo
     *  @param  clusterMergeMap the populated cluster merge map
     */
    void MakePfoMerges(const pandora::ClusterList &clusters3D, const ClusterToPfoMap &clusterToPfoMap, const ClusterMergeMap &clusterMergeMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_inputPfoListNames;                ///< The input pfo list names
    float                   m_maxMergeDistance;                 ///< he maximum distance between 3D clusters for merging to occur
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ExamplePfoMopUpAlgorithm::ClusterMerge::ClusterMerge(const pandora::Cluster *const pDaughterCluster, const float distance) :
    m_pDaughterCluster(pDaughterCluster),
    m_distance(distance)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *ExamplePfoMopUpAlgorithm::ClusterMerge::GetDaughterCluster() const
{
    return m_pDaughterCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ExamplePfoMopUpAlgorithm::ClusterMerge::GetDistance() const
{
    return m_distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ExamplePfoMopUpAlgorithm::ClusterMerge::operator<(const ClusterMerge &rhs) const
{
    // ATTN This means that a std::sort will place the cluster merge with the smallest distance first
    return (this->GetDistance() < rhs.GetDistance());
}

} // namespace lar_content

#endif // #ifndef LAR_EXAMPLE_PFO_MOP_UP_ALGORITHM_H
