/**
 *  @file   larpandoracontent/LArTwoDReco/LArSeedFinding/ShowerGrowingAlgorithm.h
 * 
 *  @brief  Header file for the shower growing algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SHOWER_GROWING_ALGORITHM_H
#define LAR_SHOWER_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArSeedFinding/SeedGrowingAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ShowerGrowingAlgorithm class
 */
class ShowerGrowingAlgorithm : public SeedGrowingAlgorithm
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

    /**
     *  @brief  Default constructor
     */
    ShowerGrowingAlgorithm();

private:
    pandora::StatusCode Run();

    typedef std::unordered_map<const pandora::Cluster*, unsigned int> ClusterInfoMap;

    /**
     *  @brief  Simple single-pass shower growing mode
     * 
     *  @param  pClusterList the list of clusters
     *  @param  clusterListName the cluster list name
     *  @param  pfoList the pfo list
     *  @param  nCaloHitsPerCluster the cluster info map
     *  @param  nBranchesPerCluster the cluster info map
     */
    void SimpleModeShowerGrowing(const pandora::ClusterList *const pClusterList, const std::string &clusterListName, pandora::PfoList &pfoList,
        ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const;

    /**
     *  @brief  Detailed recursive shower-growing and checking mode
     * 
     *  @param  pClusterList the list of clusters
     *  @param  clusterListName the cluster list name
     *  @param  pfoList the pfo list
     *  @param  nCaloHitsPerCluster the cluster info map
     *  @param  nBranchesPerCluster the cluster info map
     */
    void RecursiveModeShowerGrowing(const pandora::ClusterList *const pClusterList, const std::string &clusterListName, pandora::PfoList &pfoList,
        ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const;

    /**
     *  @brief  Get the next seed candidate, using a list of available candidates and a list of those already used
     * 
     *  @param  pClusterList the list of available seed candidates
     *  @param  usedClusters the list of candidates already considered
     *  @param  pSeedCluster to receive the address of the next seed candidate
     * 
     *  @return whether a seed candidate has been found
     */
    bool GetNextSeedCandidate(const pandora::ClusterList *const pClusterList, const pandora::ClusterList &usedClusters,
        const pandora::Cluster *&pSeedCluster) const;

    /**
     *  @brief  Get all seed candidates associated with a provided vertex
     * 
     *  @param  pClusterList the list of available seed candidates
     *  @param  pVertex the address of the vertex
     *  @param  seedClusters to receive the list of vertex seed candidates
     */
    void GetAllVertexSeedCandidates(const pandora::ClusterList *const pClusterList, const pandora::Vertex *const pVertex,
        pandora::ClusterVector &seedClusters) const;

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
     *  @brief  Check a provided seed association list for consistency, making changes as required
     * 
     *  @param  seedIter iterator to an element in the input seed association list
     *  @param  finalSeedAssociationList to receive the output seed association list
     */
    void CheckSeedAssociationList(SeedAssociationList::const_iterator seedIter, SeedAssociationList &finalSeedAssociationList) const;

    /**
     *  @brief  Process the list of branch clusters, merging with specified parent cluster, dealing with any existing pfos as required
     * 
     *  @param  pParentCluster the address of the parent cluster
     *  @param  branchClusters the list of branch clusters for the specified seed cluster
     *  @param  listName the cluster list name
     *  @param  pfoList the input pfo list
     */
    void ProcessBranchClusters(const pandora::Cluster *const pParentCluster, const pandora::ClusterVector &branchClusters, const std::string &listName,
        pandora::PfoList &pfoList) const;

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
     *  @brief  Get a figure of merit using mc information to provide best values technically possible
     *
     *  @param  seedAssociationList the seed association list
     *
     *  @return the figure of merit
     */
    float GetMCFigureOfMerit(const SeedAssociationList &seedAssociationList) const;

    /**
     *  @brief  Get a figure of merit using purely reconstructed quantities
     *
     *  @param  pVertex the address of the reconstructed 3d event vertex
     *  @param  seedAssociationList the seed association list
     *
     *  @return the figure of merit
     */
    float GetRecoFigureOfMerit(const pandora::Vertex *const pVertex, const SeedAssociationList &seedAssociationList) const;

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
     *  @brief  Whether a pointing cluster is assciated with a provided 2D vertex projection
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
     *  @brief  Get the list of all input pfos (possibly from a number of separate named lists)
     *
     *  @param  pfoList to receive the input pfo list
     */
    void GetInputPfoList(pandora::PfoList &pfoList) const;

    /**
     *  @brief  Process the details stored in a specified seed association list
     *
     *  @param  seedAssociationList the seed association list
     *  @param  clusterListName the cluster list name
     *  @param  pfoList the pfo list
     *  @param  usedClusters the list of candidates already considered
     *  @param  nCaloHitsPerCluster the cluster info map
     *  @param  nBranchesPerCluster the cluster info map
     */
    void ProcessSeedAssociationDetails(const SeedAssociationList &seedAssociationList, const std::string &clusterListName,
        pandora::PfoList &pfoList, pandora::ClusterList &usedClusters, ClusterInfoMap &nCaloHitsPerCluster, ClusterInfoMap &nBranchesPerCluster) const;

    /**
     *  @brief  Store the number of calo hits per cluster in a cluster info map
     *
     *  @param  pCluster address of the relevant cluster
     *  @param  clusterInfoMap the cluster info map
     */
    void StoreNCaloHitsPerCluster(const pandora::Cluster *const pCluster, ClusterInfoMap &clusterInfoMap) const;

    /**
     *  @brief  Store the number of branches per cluster in a cluster info map
     *
     *  @param  pCluster address of the relevant cluster
     *  @param  branchList the branch list
     *  @param  clusterInfoMap the cluster info map
     */
    void StoreNBranchesPerCluster(const pandora::Cluster *const pCluster, const pandora::ClusterVector &branchList, ClusterInfoMap &clusterInfoMap) const;

    /**
     *  @brief  Remove pfos containing clusters to which many shower branches have been added
     *
     *  @param  nCaloHitsPerCluster the cluster info map
     *  @param  nBranchesPerCluster the cluster info map
     *  @param  pfoList the input pfo list
     */
    void RemoveShowerPfos(const ClusterInfoMap &nCaloHitsPerCluster, const ClusterInfoMap &nBranchesPerCluster, pandora::PfoList &pfoList) const;

    /**
     *  @brief  Get the address of the pfo containing a specified cluster
     *
     *  @param  pCluster address of the relevant cluster
     *  @param  pfoList the list of all input pfos
     *  @param  pTargetPfo to receive the address of the target pfo
     */
    void FindTargetPfo(const pandora::Cluster *const pCluster, const pandora::PfoList &pfoList, const pandora::Pfo *&pTargetPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::Cluster*, LArVertexHelper::ClusterDirection> ClusterDirectionMap;
    mutable ClusterDirectionMap m_clusterDirectionMap;          ///< The cluster direction map

    pandora::StringVector       m_inputClusterListNames;        ///< The names of the input cluster lists
    pandora::StringVector       m_inputPfoListNames;            ///< The names of the input pfo lists

    unsigned int                m_minCaloHitsPerCluster;        ///< The minimum number of calo hits per (seed or branch) cluster
    float                       m_nearbyTrackDistance;          ///< Prevent track-track associations where the end-to-end separation is smaller than this distance
    float                       m_nearbyClusterDistance;        ///< The nearby cluster distance, used for determining cluster associations
    float                       m_remoteClusterDistance;        ///< The remote cluster distance, used for determining cluster associations

    float                       m_directionTanAngle;            ///< Direction determination, look for vertex inside triangle with apex shifted along the cluster length
    float                       m_directionApexShift;           ///< Direction determination, look for vertex inside triangle with apex shifted along the cluster length

    bool                        m_recursiveMode;                ///< Whether to run in recursive shower-growing and checking mode (rather than a simple one-pass mode)
    bool                        m_useFirstImprovedSeed;         ///< Recursive mode only: use the first daughter seed (from an ordered list) that offers an improved figure of merit
    bool                        m_useMCFigureOfMerit;           ///< Recursive mode only: use a figure of merit based on mc particle information

    bool                        m_shouldRemoveShowerPfos;       ///< Whether to delete any existing pfos to which many shower branches have been added
    unsigned int                m_showerLikeNBranches;          ///< The minimum number of branches before cluster is declared shower like
    float                       m_showerLikeCaloHitRatio;       ///< The minimum ratio of final to original calo hits before cluster is declared shower like

    float                       m_minVertexLongitudinalDistance;///< Vertex association check: min longitudinal distance cut
    float                       m_maxVertexLongitudinalDistance;///< Vertex association check: max longitudinal distance cut
    float                       m_maxVertexTransverseDistance;  ///< Vertex association check: max transverse distance cut
    float                       m_vertexAngularAllowance;       ///< Vertex association check: pointing angular allowance in degrees
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ShowerGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ShowerGrowingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_GROWING_ALGORITHM_H
