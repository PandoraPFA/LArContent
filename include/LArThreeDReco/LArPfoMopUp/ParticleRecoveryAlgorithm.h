/**
 *  @file   LArContent/include/LArThreeDReco/LArPfoMopUp/ParticleRecoveryAlgorithm.h
 * 
 *  @brief  Header file for the track recovery algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PARTICLE_RECOVERY_ALGORITHM_H
#define LAR_PARTICLE_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ParticleRecoveryAlgorithm class
 */
class ParticleRecoveryAlgorithm : public pandora::Algorithm
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
    ParticleRecoveryAlgorithm();

private:
    /**
     *  @brief  SimpleOverlapTensor class
     */
    class SimpleOverlapTensor
    {
    public:
        /**
         *  @brief  Add an association between two clusters to the simple overlap tensor
         * 
         *  @param  pCluster1 address of cluster 1
         *  @param  pCluster2 address of cluster 2
         */
        void AddAssociation(pandora::Cluster *pCluster1, pandora::Cluster *pCluster2);

        /**
         *  @brief  Get elements connected to a specified cluster
         * 
         *  @param  pCluster address of the cluster
         *  @param  elementList the element list
         *  @param  clusterListU connected u clusters
         *  @param  clusterListV connected v clusters
         *  @param  clusterListW connected w clusters
         */
        void GetConnectedElements(pandora::Cluster *const pCluster, const bool ignoreUnavailable, pandora::ClusterList &clusterListU,
            pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

        /**
         *  @brief  Get the list of key clusters
         *
         *  @return the list of key clusters
         */
        const pandora::ClusterList &GetKeyClusters() const;

    private:
        typedef std::map<pandora::Cluster*, pandora::ClusterList> ClusterNavigationMap;

        pandora::ClusterList    m_keyClusters;                  ///< The list of key clusters
        ClusterNavigationMap    m_clusterNavigationMapUV;       ///< The cluster navigation map U->V
        ClusterNavigationMap    m_clusterNavigationMapVW;       ///< The cluster navigation map V->W
        ClusterNavigationMap    m_clusterNavigationMapWU;       ///< The cluster navigation map W->U
    };

    pandora::StatusCode Run();

    /**
     *  @brief  Get the input cluster lists for processing in this algorithm
     *
     *  @param  inputClusterListU to receive the list of clusters in the u view
     *  @param  inputClusterListU to receive the list of clusters in the v view
     *  @param  inputClusterListU to receive the list of clusters in the w view
     */
    void GetInputClusters(pandora::ClusterList &inputClusterListU, pandora::ClusterList &inputClusterListV, pandora::ClusterList &inputClusterListW) const;

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     *
     *  @param  inputClusterList the input cluster list
     *  @param  selectedClusterList to receive the selected cluster list
     */
    void SelectInputClusters(const pandora::ClusterList &inputClusterList, pandora::ClusterList &selectedClusterList) const;

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     *
     *  @param  inputClusterList the input cluster list
     *  @param  selectedClusterList to receive the selected cluster list
     */
    void StandardClusterSelection(const pandora::ClusterList &inputClusterList, pandora::ClusterList &selectedClusterList) const;

    /**
     *  @brief  Select a subset of input clusters nodally associated with the vertices of existing particles
     *
     *  @param  inputClusterList the input cluster list
     *  @param  selectedClusterList to receive the selected cluster list
     */
    void VertexClusterSelection(const pandora::ClusterList &inputClusterList, pandora::ClusterList &selectedClusterList) const;

    /**
     *  @brief  Find cluster overlaps and record these in the overlap tensor
     *
     *  @param  clusterList1 the first cluster list
     *  @param  clusterList2 the second cluster list
     *  @param  overlapTensor the overlap tensor
     */
    void FindOverlaps(const pandora::ClusterList &clusterList1, const pandora::ClusterList &clusterList2, SimpleOverlapTensor &overlapTensor) const;

    /**
     *  @brief  Whether two clusters overlap convincingly in x
     *
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     */
    bool IsOverlap(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Identify unambiguous cluster overlaps and resolve ambiguous overlaps, creating new track particles
     *
     *  @param  overlapTensor the overlap tensor
     */
    void ExamineTensor(const SimpleOverlapTensor &overlapTensor) const;

    /**
     *  @brief  Whether a trio of clusters are consistent with representing projections of the same 3d trajectory
     *
     *  @param  pClusterU the address of cluster u
     *  @param  pClusterV the address of cluster v
     *  @param  pClusterW the address of cluster w
     * 
     *  @return boolean
     */
    bool CheckConsistency(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW) const;

    /**
     *  @brief  Create and save a track particle containing the provided clusters
     *
     *  @param  clusterList the cluster list
     */
    void CreateTrackParticle(const pandora::ClusterList &clusterList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector       m_inputClusterListNames;        ///< The list of cluster list names
    std::string                 m_outputPfoListName;            ///< The output pfo list name

    bool                        m_includeTracks;                ///< Whether to include fixed tracks in selected cluster list
    bool                        m_includeShowers;               ///< Whether to include clusters not fixed as tracks in selected cluster list

    unsigned int                m_minClusterCaloHits;           ///< The min number of hits in base cluster selection method
    float                       m_minClusterLengthSquared;      ///< The min length (squared) in base cluster selection method
    float                       m_minClusterXSpan;              ///< The min x span required in order to consider a cluster

    bool                        m_vertexClusterMode;            ///< Whether to demand clusters are associated with vertices of existing particles
    float                       m_minVertexLongitudinalDistance;///< Vertex association check: min longitudinal distance cut
    float                       m_maxVertexTransverseDistance;  ///< Vertex association check: max transverse distance cut

    float                       m_minXOverlapFraction;          ///< The min x overlap fraction required in order to id overlapping clusters
    float                       m_pseudoChi2Cut;                ///< The selection cut on the matched chi2
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ParticleRecoveryAlgorithm::Factory::CreateAlgorithm() const
{
    return new ParticleRecoveryAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterList &ParticleRecoveryAlgorithm::SimpleOverlapTensor::GetKeyClusters() const
{
    return m_keyClusters;
}

} // namespace lar_content

#endif // #ifndef LAR_PARTICLE_RECOVERY_ALGORITHM_H
