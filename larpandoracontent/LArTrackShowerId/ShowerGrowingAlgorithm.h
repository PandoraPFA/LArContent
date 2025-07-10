/**
 *  @file   larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the shower growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_GROWING_ALGORITHM_H
#define LAR_SHOWER_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Objects/CaloHit.h"

#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTrackShowerId/BranchGrowingAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ShowerGrowingAlgorithm class
 */
class ShowerGrowingAlgorithm : public BranchGrowingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerGrowingAlgorithm();

    typedef std::unordered_map<const pandora::Cluster *, float> ClusterLengthMap;

protected:
    /**
     *  @brief  Whether a pointing cluster is associated with a provided 2D vertex projection
     *
     *  @param  pointingCluster the pointing cluster
     *  @param  vertexPosition2D the projected vertex position
     *
     *  @return boolean
     */
    bool IsVertexAssociated(const LArPointingCluster &pointingCluster, const pandora::CartesianVector &vertexPosition2D) const;

    /**
     *  @brief  Sorting for clusters to determine order in which seeds are considered
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortClusters(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Cluster sort class, with the ability to use a cluster length cache.
     */
    class ClusterSeedComparator
    {
    public:
        ClusterSeedComparator(const ClusterLengthMap *const pClusterLengthCache) : m_pClusterLengthCache(pClusterLengthCache) {}

        bool operator()(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs) const;

    private:
        const ClusterLengthMap* m_pClusterLengthCache; ///< The cluster length cache
    };

    typedef std::unordered_map<const pandora::Cluster *, LArVertexHelper::ClusterDirection> ClusterDirectionMap;
    mutable ClusterDirectionMap m_clusterDirectionMap; ///< The cluster direction map

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Simple single-pass shower growing mode
     *
     *  @param  pClusterList the list of clusters
     *  @param  clusterListName the cluster list name
     */
    void SimpleModeShowerGrowing(const pandora::ClusterList *const pClusterList, const std::string &clusterListName) const;

    /**
     *  @brief  Get the next seed candidate, using a list of available candidates and a list of those already used
     *
     *  @param  pClusterList the list of available seed candidates
     *  @param  usedClusters the list of candidates already considered
     *  @param  pSeedCluster to receive the address of the next seed candidate
     *
     *  @return whether a seed candidate has been found
     */
    bool GetNextSeedCandidate(const pandora::ClusterList *const pClusterList, const pandora::ClusterSet &usedClusters,
        const pandora::Cluster *&pSeedCluster) const;

    /**
     *  @brief  Get all seed candidates associated with a provided vertex
     *
     *  @param  pClusterList the list of available seed candidates
     *  @param  pVertex the address of the vertex
     *  @param  seedClusters to receive the list of vertex seed candidates
     */
    void GetAllVertexSeedCandidates(
        const pandora::ClusterList *const pClusterList, const pandora::Vertex *const pVertex, pandora::ClusterVector &seedClusters) const;

    /**
     *  @brief  Get the seed association list for a given vector of particle seed candidates
     *
     *  @param  particleSeedVector the particle seed vector
     *  @param  pClusterList the address of the input cluster list
     *  @param  seedAssociationList to receive the populated seed association list
     */
    void GetSeedAssociationList(const pandora::ClusterVector &particleSeedVector, const pandora::ClusterList *const pClusterList,
        SeedAssociationList &seedAssociationList) const;

    /**
     *  @brief  Process the details stored in a specified seed association list
     *
     *  @param  seedAssociationList the seed association list
     *  @param  clusterListName the cluster list name
     *  @param  pfoList the pfo list
     *  @param  usedClusters the list of candidates already considered
     */
    void ProcessSeedAssociationDetails(
        const SeedAssociationList &seedAssociationList, const std::string &clusterListName, pandora::ClusterSet &usedClusters) const;

    /**
     *  @brief  Process the list of branch clusters, merging with specified parent cluster, dealing with any existing pfos as required
     *
     *  @param  pParentCluster the address of the parent cluster
     *  @param  branchClusters the list of branch clusters for the specified seed cluster
     *  @param  listName the cluster list name
     *  @param  pfoList the input pfo list
     */
    void ProcessBranchClusters(
        const pandora::Cluster *const pParentCluster, const pandora::ClusterVector &branchClusters, const std::string &listName) const;

    AssociationType AreClustersAssociated(const pandora::Cluster *const pClusterSeed, const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get a figure of merit representing the consistency of the provided seed associated list
     *
     *  @param  seedAssociationList the seed association list
     *
     *  @return the figure of merit
     */
    float GetFigureOfMerit(const SeedAssociationList &seedAssociationList) const;

    /**
     *  @brief  Get the number of clusters associated with the vertex
     *
     *  @param  vertexPosition2D the projected vertex position
     *  @param  pointingClusterList the list of relevant pointing clusters
     *
     *  @return the number of clusters associated with the vertex
     */
    unsigned int GetNVertexConnections(const pandora::CartesianVector &vertexPosition2D, const LArPointingClusterList &pointingClusterList) const;

    /**
     *  @brief  Pre-compute the lengths of all clusters in the input list, if caching is enabled
     *
     *  @param  pClusterList the address of the input cluster list
     */
    void PreComputeClusterLengths(const pandora::ClusterList *const pClusterList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputClusterListNames; ///< The names of the input cluster lists

    unsigned int m_minCaloHitsPerCluster; ///< The minimum number of calo hits per (seed or branch) cluster
    float m_nearbyClusterDistance;        ///< The nearby cluster distance, used for determining cluster associations
    float m_remoteClusterDistance;        ///< The remote cluster distance, used for determining cluster associations

    float m_directionTanAngle;  ///< Direction determination, look for vertex inside triangle with apex shifted along the cluster length
    float m_directionApexShift; ///< Direction determination, look for vertex inside triangle with apex shifted along the cluster length

    float m_minVertexLongitudinalDistance; ///< Vertex association check: min longitudinal distance cut
    float m_maxVertexLongitudinalDistance; ///< Vertex association check: max longitudinal distance cut
    float m_maxVertexTransverseDistance;   ///< Vertex association check: max transverse distance cut
    float m_vertexAngularAllowance;        ///< Vertex association check: pointing angular allowance in degrees

    mutable ClusterLengthMap m_clusterLengthCache; ///< The cluster length cache
    bool m_useClusterLengthCache;                  ///< Whether to use the cluster length cache
};

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_GROWING_ALGORITHM_H
