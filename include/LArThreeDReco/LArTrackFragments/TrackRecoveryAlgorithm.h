/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackFragments/TrackRecoveryAlgorithm.h
 * 
 *  @brief  Header file for the track recovery algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRACK_RECOVERY_ALGORITHM_H
#define LAR_TRACK_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  TrackRecoveryAlgorithm class
 */
class TrackRecoveryAlgorithm : public pandora::Algorithm
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
    TrackRecoveryAlgorithm();

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
     *  @brief  Select a subset of input clusters for processing in this algorithm
     *
     *  @param  inputClusterListName the input cluster list name
     *  @param  selectedClusterList to receive the selected cluster list
     */
    void SelectInputClusters(const std::string &inputClusterListName, pandora::ClusterList &selectedClusterList) const;

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
     *  @brief  Create and save a track particle containing the provided clusters
     *
     *  @param  clusterList the cluster list
     */
    void CreateTrackParticle(const pandora::ClusterList &clusterList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string                 m_outputPfoListName;            ///< The output pfo list name

    unsigned int                m_minClusterCaloHits;           ///< The min number of hits in base cluster selection method
    float                       m_minClusterLengthSquared;      ///< The min length (squared) in base cluster selection method
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TrackRecoveryAlgorithm::Factory::CreateAlgorithm() const
{
    return new TrackRecoveryAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::ClusterList &TrackRecoveryAlgorithm::SimpleOverlapTensor::GetKeyClusters() const
{
    return m_keyClusters;
}

} // namespace lar_content

#endif // #ifndef LAR_TRACK_RECOVERY_ALGORITHM_H
