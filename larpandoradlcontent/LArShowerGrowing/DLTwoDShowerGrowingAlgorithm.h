/**
 *  @file   larpandoradlcontent/LArShowerGrowing/DLTwoDShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the deep learning shower growing algorithm
 *
 *  $Log: $
 */
#ifndef LAR_DL_TWOD_SHOWER_GROWING_ALGORITHM
#define LAR_DL_TWOD_SHOWER_GROWING_ALGORITHM 1

#include <torch/torch.h>

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_dl_content
{
/**
 *  @brief  DLTwoDShowerGrowingAlgorithm class
 */
class DLTwoDShowerGrowingAlgorithm : public pandora::Algorithm
{
private:
    struct HitFeatures
    {
        HitFeatures();

        float m_xRel;
        float m_zRel;
        float m_rRel;
        float m_cosThetaRel;
        float m_sinThetaRel;
        float m_distToXGap;
        float m_xWidth; 
        float m_energy;
    };

    struct ClusterGroup
    {
        ClusterGroup();

        /**
         *  @brief Add a cluster to the group. If first cluster in the group, use it as the group's representative cluster
         */
        void insert(const pandora::Cluster *pCluster);

        const pandora::Cluster *GetRepresentativeCluster() const { return m_representativeCluster; }
        const pandora::ClusterSet &GetClusters() const { return m_clusters; }
        size_t size() const { return m_clusters.size(); }
        bool empty() const { return m_clusters.empty(); }
        auto begin() const { return m_clusters.begin(); }
        auto end() const { return m_clusters.end(); }

        pandora::ClusterSet m_clusters;                  ///< Container of clusters in the group
        const pandora::Cluster *m_representativeCluster; ///< The address of an arbitrary cluster in the group to use as a key for the group
    };

public:
    /**
     *  @brief Default constructor
     */
    DLTwoDShowerGrowingAlgorithm();

    /**
     *  @brief Default destructor
     */
    virtual ~DLTwoDShowerGrowingAlgorithm();

private:
    typedef std::map<const pandora::Cluster *const, std::map<const pandora::Cluster *const, float>> SimilarityMatrix;
    typedef std::map<const pandora::Cluster *const, std::vector<const pandora::Cluster *>> AdjacencyLists;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Create output TTree for generating the training data. Note that for best results simulation should have EM roll-up
     *         turned off to ensure delta rays are not folded into their parent tracks.
     */ 
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief Do an inference. Includes network similarity prediction and subsequent merges.
     */ 
    pandora::StatusCode Infer();

    /* Start training sample preparation methods */

    /**
     *  @brief Get the folded MCParticle that contributes the most to the CaloHit.
     *
     *  @param[in]     pCaloHit The hit
     *  @param[in,out] mcFoldTo A map for the folding of MCParticles, updates as new MCParticles are seen by the method
     *
     *  @return The folded MCParticle that contributed the most to the CaloHit
     */ 
    const pandora::MCParticle* GetMainMC(
        const pandora::CaloHit *const pCaloHit,
        std::map<const pandora::MCParticle *const, const pandora::MCParticle *const> &mcFoldTo) const;

    /**
     *  @brief Finds the ancestor that a MCParticle should be folded to in order to roll-up an EM shower
     *
     *  @param[in] pMC The MCParticle
     *
     *  @return The ancestor MCParticle, this will be the input pMC an MCParticle that should not be rolled-up
     */
    const pandora::MCParticle* FoldMCTo(const pandora::MCParticle *const pMC) const;

    /**
     *  @brief Recursive function to check descendent particles for signature of a shower (e -> gamma -> e).
     *         Used to identify delta-rays that shower and photons that only undergo compton scatters leaving diffuse hits.
     *
     *  @param[in] pMC                  The MCParticle to check descendents of
     *  @param[in] nDescendentElectrons The number of electrons seen while descending from the original photon MCParticle
     *
     *  @return Flag to indicate if more than 1 electron was seen while descending any of the descendent particle association paths
     */
    bool CausesShower(const pandora::MCParticle *const pMC, int nDescendentElectrons) const;

    /**
     *  @brief Tries to identify and deal with impossible-to-cluster-correctly delta ray hits by assigning the hit to the
     *         parent particle's cluster if the delta-ray is short and does not shower
     *         or the hit has charge contribution from the parent particle.
     *
     *  @param[in] pCaloHit The hit
     *  @param[in] pMC      The main MCParticle of the hit
     *
     *  @return The MCParticle the hit should be assigned to, this will either be the inputted MCParticle or the parent MCParticle
     */
    const pandora::MCParticle* FoldPotentialDeltaRayTo(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pMC) const;

    /* End training sample preparation methods */

    /* Start general helpers */

    /**
     *  @brief Get the 2D reconstructed neutrino vertices for each view.
     *
     *  @return Map from view to the 3D neutrino vertex position projected onto the view
     */
    std::map<pandora::HitType, pandora::CartesianVector> Get2DVertices() const;

    /**
     *  @brief Get clusters from all lists.
     *
     *  @return List of clusters from all list names
     */
    pandora::ClusterList GetAllClusters() const;

    /**
     *  @brief Get clusters from specified list.
     *
     *  @param[in]  clusterListName List name for clusters
     *  @param[out] clusterList     List of clusters 
     */
    pandora::StatusCode GetClusters(const std::string clusterListName, pandora::ClusterList &clusterList) const;

    /**
     *  @brief Calculate all properties of a hit relevant for constructing the hit feature vector.
     *
     *  @param[in]  pCaloHit    The hit
     *  @param[in]  vtxPos      The position of the 2D neutrino vertex associated with the hit
     *  @param[out] hitFeatures Struct to store calculated hit properties
     */
    pandora::StatusCode CalculateHitFeatures(
        const pandora::CaloHit *const pCaloHit, const pandora::CartesianVector vtxPos, HitFeatures &hitFeatures) const;

    /* End general helpers */

    /* Start inference methods */

    /**
     *  @brief Runs inference on the models to predict the similarity matrix for a 2D ClusterList.
     *
     *  @param[in]  clusterList   The list of 2D clusters
     *  @param[in]  view          The view of the 2D clusters
     *  @param[in]  vtxPos        The 2D neutrino vertex of the view
     *  @param[in]  iterationNum  The iteration of the inference (1 is the first)
     *  @param[out] clusterSimMat The similarity matrix encoding all predicted pairwise similarities of the input list of 2D clusters
     */
    pandora::StatusCode PredictClusterSimilarityMatrix(
        const pandora::ClusterList &clusterList,
        const pandora::HitType view,
        const pandora::CartesianVector &vtxPos,
        const unsigned int iterationNum,
        SimilarityMatrix &clusterSimMat);

    /**
     *  @brief Construct the tensor that encodes a cluster.
     *
     *  @param[in]  clusterFeatures Vector of hit properties for the hits of one 2D cluster
     *  @param[in]  view            The view of the 2D clusters
     *  @param[in]  nClusters       The total number of 2D clusters in the view
     *  @param[in]  iterationNum    The iteration of the inference (1 is the first)
     *  @param[out] tensorCluster   The tensor encoding the cluster as a sequence of hit feature vectors
     */
    pandora::StatusCode MakeClusterTensor(
        const std::vector<HitFeatures> &clusterFeatures,
        const pandora::HitType view,
        const size_t nClusters,
        const unsigned int iterationNum,
        torch::Tensor &tensorCluster) const;

    /**
     *  @brief Populates a SimilarityMatrix from the Tensor output of model inference.
     *
     *  @param[in]  tensorSimMat  The predicted similarity matrix tensor, shape [1, N, N] where N is the number of 2D clusters
     *  @param[in]  clusterList   The list of 2D clusters
     *  @param[out] clusterSimMat The cluster similarity matrix object
     */
    pandora::StatusCode PopulateClusterSimilarityMatrix(
        const torch::Tensor &tensorSimMat, const pandora::ClusterList &clusterList, SimilarityMatrix &clusterSimMat) const;

    /**
     *  @brief Uses the predicted similarity matrix to partition the 2D clusters into groups that should be merged. First, core (larger) 
     *         clusters are selected. Their similarities are thresholded to get edges. From these edges, the core clusters are partitioned
     *         by their connected subgraphs. Then, accessory (smaller) clusters are added to the groups of the core cluster partition.
     *         Finally, any remaining accessory clusters are partitioned via the connected subgraphs method.
     *
     *  @param[in]  clusterSimMat       The predicted cluster similarity matrix
     *  @param[in]  similarityThreshold The threshold to convert predicted similarities into edges
     *  @param[out] clusterGroups       The resulting partition of the 2D clusters
     */
    pandora::StatusCode ClusterFromSimilarity(
        const SimilarityMatrix &clusterSimMat, const float similarityThreshold, std::vector<ClusterGroup> &clusterGroups) const;

    /**
     *  @brief Split clusters into accessory and core clusters and convert the similarity into adjacencies (binary edges).
     *
     *  @param[in]  simMat              The similarity matrix
     *  @param[in]  similarityThreshold The threshold to convert similarities into adjacencies
     *  @param[out] coreClusterAdjLists The resulting core cluster adjacencies
     *  @param[out] accClusterAdjLists  The resulting accessory cluster adjacencies
     */
    pandora::StatusCode PopulateAdjacencyLists(
        const SimilarityMatrix &simMat,
        const float similarityThreshold,
        AdjacencyLists &coreClusterAdjLists,
        AdjacencyLists &accClusterAdjLists) const;

    /**
     *  @brief Generate a partition from cluster adjacencies (binary edges) using connected subgraphs.
     *
     *  @param[in]  clusterAdjLists The cluster adjacencies
     *  @param[out] clusterGroups   The resulting 2D cluster partition
     */
    pandora::StatusCode CalculateConnectedGroups(const AdjacencyLists &clusterAdjLists, std::vector<ClusterGroup> &clusterGroups) const;

    /**
     *  @brief Adds clusters to groups of an existing cluster partition. The candidate clusters are added to a group of the partition if
     *         their maximum similarity with any cluster in the group exceeds a threshold. This happens iteratively in case the addition of
     *         a candidate cluster to a group changes the maximum similarity of another candidate cluster to that group. The purpose is to
     *         mop up accessory clusters without allowing them to merge together groups of the current partition.
     *
     *  @param[in]     clusterSimMat            The cluster similarity matrix
     *  @param[in]     ungroupedClusterAdjLists The candidate clusters that are not yet part of the partition
     *  @param[in]     similarityThreshold      The similarity threshold
     *  @param[in,out] clusterGroups            The initial cluster partition that may have clusters added to it
     */
    pandora::StatusCode AddClustersToGroups(
        const SimilarityMatrix &clusterSimMat,
        AdjacencyLists &ungroupedClusterAdjLists,
        const float similarityThreshold,
        std::vector<ClusterGroup> &clusterGroups) const;

    /**
     *  @brief Populates a super-cluster similarity matrix from an initial similarity matrix and the cluster partition. The similarities
     *         between super-clusters (groups of the partition) is the minimum of all pairwise cluster similarities between constituent
     *         clusters of the super-clusters.
     *
     *  @param[in]  clusterGroups      The cluster partition defining the super-clusters
     *  @param[in]  clusterSimMat      The similarity matrix of the individual clusters of the partition
     *  @param[out] superClusterSimMat The super-cluster similarity matrix
     */
    pandora::StatusCode PopulateSuperClusterSimilarityMatrix(
        const std::vector<ClusterGroup> &clusterGroups, const SimilarityMatrix &clusterSimMat, SimilarityMatrix &superClusterSimMat) const;

    /**
     *  @brief Merge clusters according a partition.
     *
     *  @param[in] clusterGroups The cluster partition
     *  @param[in] listName      the list name of the clusters being merged
     */
    pandora::StatusCode MergeGroups(const std::vector<ClusterGroup> &clusterGroups, const std::string &listNmae) const;

    /**
     *  @brief Check if a cluster partition is invalid, valid being disjoint and covers the set of clusters.
     *
     *  @param[in] clusterGroups The cluster partition
     *  @param[in] clusterList   The list clusters
     *
     *  @return Flag for if partition is invalid
     */
    bool IsInvalidPartition(const std::vector<ClusterGroup> &clusterGroups, const pandora::ClusterList  &clusterList) const;

    /**
     *  @brief Check if a cluster partition is all singletons.
     *
     *  @param[in]  clusterGroups The cluster partition
     *
     *  @return  Flag for if partition is all singletons
     */
    bool IsSingletonPartition(const std::vector<ClusterGroup> &clusterGroups) const;

    /* End inference methods */

    /* Start shared mutable members */

    LArDLHelper::TorchModel m_modelEncoder; ///< TorchScript model for encoding hits in a cluster
    LArDLHelper::TorchModel m_modelAttn;    ///< TorchScript model for attention over encoded clusters in a view
    LArDLHelper::TorchModel m_modelSim;     ///< TorchScript model for pairwise similarities over clusters after attention
    float m_polarRScaleFactor;              ///< Scale factor for polar r coordinate input features
    float m_cartesianXScaleFactor;          ///< Scale factor for cartesian x coordinate input features
    float m_cartesianZScaleFactor;          ///< Scale factor for cartesian z coordinate input features
    std::set<double> m_detectorXGaps;       ///< X coordinates where gaps in X direction start/end
    int m_hitFeaturesNHitsIdx;              ///< Hit feature vector index for the optional no. hits in cluster feature
    int m_hitFeaturesNClustersIdx;          ///< Hit feature vector index for the optional no. clusters in event feature
    int m_hitFeaturesIterationNumIdx;       ///< Hit feature vector index for the optional no. clusters in event feature

    /** End shared mutable members ***/

    /** Start hardcoded members ***/

    std::map<pandora::HitType, float> m_deltaRayLengthThresholdSquared; ///< Threshold for defining small delta rays that will be folded to the parent particle
    float m_deltaRayParentWeightThreshold;                              ///< Threshold for weight contribution of parent particle for it take the delta ray's hit
    int m_hitFeatureDim;                                                ///< Feature dimensions of each hit

    /* End hardcoded members */

    /* Start configurable via xml members */

    bool m_trainingMode;                       ///< Training mode
    std::string m_trainingTreeName;            ///< Tree name for training data output
    std::string m_trainingFileName;            ///< File name for training data output
    pandora::StringVector m_clusterListNames;  ///< Names of cluster lists
    std::string m_vertexListName;              ///< Name of vertex list
    float m_similarityThreshold;               ///< Threshold value on similarity for clusters to be connected
    float m_similarityThresholdBeta;           ///< Scaling factor for an optional second clustering pass
    unsigned int m_accessoryClustersMaxHits;   ///< Clusters with this number of less hits are treated as accessory clusters during merging
    unsigned int m_maxIterations;              ///< Max iterative applications of the cluster merging
    bool m_includeHitCardinalityFeatures;      ///< Option to include in hit feature vector the no. hits in cluster and no. clusters in event
    bool m_includeHitNIterationNumFeature;     ///< Option to include in hit feature vector the inference iteration number

    /* End configurable via xml members */
};

} // namespace lar_dl_content

#endif // LAR_DL_TWOD_SHOWER_GROWING_ALGORITHM
