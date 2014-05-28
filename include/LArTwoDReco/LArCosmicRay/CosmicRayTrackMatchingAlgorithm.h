/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
#define LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

namespace lar
{

/**
 *  @brief  CosmicRayTrackMatchingAlgorithm class
 */
class CosmicRayTrackMatchingAlgorithm : public pandora::Algorithm
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

    typedef std::set<unsigned int>                                  IntList;
    typedef std::map<const pandora::CaloHit*, pandora::Cluster*>    HitToClusterMap;
    typedef std::map<const pandora::Cluster*, pandora::CaloHitList> ClusterToHitMap;
    typedef std::map<const pandora::CaloHit*, IntList>              HitAssociationMap;
    typedef std::map<const unsigned int, pandora::CaloHitList>      ClusterAssociationMap;

    /**
     *  @brief Get a vector of available clusters
     *
     *  @param inputClusterListName the input name of the cluster list
     *  @param clusterVector the output vector of available clusters
     */
    pandora::StatusCode GetAvailableClusters(const std::string inputClusterListName, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief Select a set of clusters judged to be clean
     *
     *  @param inputVector the input vector of all available clusters
     *  @param outputVector the output vector of clean clusters
     */
    void SelectCleanClusters(const pandora::ClusterVector &inputVector, pandora::ClusterVector &outputVector) const;

    /**
     *  @brief Generate a map of sliding linear fit results from a vector of clusters
     *
     *  @param clusterVector the input vector of clusters
     *  @param slidingFitResultMap the output map of sliding linear fit results
     */
    void AddToSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const;

    /**
     *  @brief Match clusters bewteen views
     *
     *  @param slidingFitResultMap the input map of sliding linear fit results
     *  @param clusterVector1 the input vector of clusters for the first view
     *  @param clusterVector2 the input vector of clusters for the second view
     *  @param clusterVector3 the input vector of clusters for the third view
     *  @param hitAssociationMap the output map from associated hits to new clusters
     *  @param clusterAssociationMap the output map from new cluster to associated hits
     */
    void SelectMatchedTracks(const TwoDSlidingFitResultMap &slidingFitResultMap, const pandora::ClusterVector &clusterVector1,
        const pandora::ClusterVector &clusterVector2, const pandora::ClusterVector &clusterVector3,
        HitAssociationMap &hitAssociationMap, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief Match clusters between views
     *
     *  @param clusterID the unique ID number for this match
     *  @param slidingFitResult1 the input sliding linear fit for the first view
     *  @param slidingFitResult2 the input sliding linear fit for the second view
     *  @param availableClusters3 the input vector of available clusters for the third view
     *  @param hitAssociationMap the output map from associated hits to new clusters
     *  @param clusterAssociationMap the output map from new cluster to associated hits
     */
    void SelectMatchedTracks(const unsigned int clusterID, const TwoDSlidingFitResult &slidingFitResult1,
        const TwoDSlidingFitResult &slidingFitResult2, const pandora::ClusterVector &availableClusters3,
        HitAssociationMap &hitAssociationMap, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief Project clusters from two views into a third view using transverse matching
     *
     *  @param slidingFitResult1 the sliding fit result in the first view
     *  @param slidingFitResult2 the sliding fit result in the second view 
     *  @param positionList the projected points in the third view
     */
    void SelectTransverseMatchedPoints(const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2,
        pandora::CartesianPointList &positionList) const;

    /**
     *  @brief Project clusters from two views into a third view using longitudinal matching
     *
     *  @param slidingFitResult1 the sliding fit result in the first view
     *  @param slidingFitResult2 the sliding fit result in the second view 
     *  @param positionList the projected points in the third view
     */
    void SelectLongitudinalMatchedPoints(const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2,
        pandora::CartesianPointList &positionList) const;

    /**
     *  @brief Generate lists of matched hits and positions
     *
     *  @param availableClusters input vector of available clusters
     *  @param projectedPositions input vector of projected points
     *  @param matchedHits output list of matched hits
     *  @param matchedClusters output list of matched clusters
     *  @param matchedPositions output list of matched positions
     */
    void SelectMatchedHits(const pandora::ClusterVector &availableClusters, const pandora::CartesianPointList &projectedPositions, 
        pandora::CaloHitList &matchedHits, pandora::ClusterList &matchedClusters, pandora::CartesianPointList &matchedPositions) const;

    /**
     *  @brief Modify existing clusters and create new clusters
     *
     *  @param inputClusterListName the cluster list name for this view
     *  @param hitAssociationMap the mapping from hits to new clusters
     *  @param clusterAssociationMap the mapping from new clusters to hits
     */
    pandora::StatusCode ModifyClusters(const std::string inputClusterListName, HitAssociationMap &hitAssociationMap,
        ClusterAssociationMap &clusterAssociationMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string    m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string    m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string    m_inputClusterListNameW;        ///< The name of the view W cluster list

    bool           m_useTransverseMode;            ///< Run transverse matching between views
    bool           m_useLongitudinalMode;          ///< Run longitudinal matching between views

    float          m_clusterMinLength;             ///< minimum length of clusters for this algorithm
    float          m_minXOverlap;                  ///< requirement on minimum X overlap for associated clusters
    float          m_minXOverlapFraction;          ///< requirement on minimum X overlap fraction for associated clusters
    float          m_maxPointDisplacement;         ///< maximum distance between projected track and associated hits
    float          m_maxHitDisplacement;           ///< maximum distance between associated hits
    float          m_minMatchedPointFraction;      ///< minimum fraction of matched points

    unsigned int   m_halfWindowLayers;             ///< number of layers to use for half-window of sliding fit 
    unsigned int   m_numSamplingPoints;            ///< The number of matched points to be generated
    unsigned int   m_minMatchedHits;               ///< minimum number of matched hits
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayTrackMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayTrackMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_TRACK_MATCHING_ALGORITHM_H
