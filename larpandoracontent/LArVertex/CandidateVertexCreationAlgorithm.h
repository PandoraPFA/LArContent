/**
 *  @file   larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.h
 *
 *  @brief  Header file for the candidate vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
#define LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CandidateVertexCreationAlgorithm::Algorithm class
 */
class CandidateVertexCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CandidateVertexCreationAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Select a subset of input clusters (contained in the input list names) for processing in this algorithm
     *
     *  @param  clusterVectorU to receive the selected clusters in the u view
     *  @param  clusterVectorV to receive the selected clusters in the v view
     *  @param  clusterVectorW to receive the selected clusters in the w view
     */
    void SelectClusters(pandora::ClusterVector &clusterVectorU, pandora::ClusterVector &clusterVectorV, pandora::ClusterVector &clusterVectorW);

    /**
     *  @brief  Create candidate vertex positions by comparing pairs of cluster end positions
     *
     *  @param  clusterVector1 the clusters in view 1
     *  @param  clusterVector1 the clusters in view 2
     */
    void CreateEndpointCandidates(const pandora::ClusterVector &clusterVector1, const pandora::ClusterVector &clusterVector2) const;

    /**
     *  @brief  Create a candidate vertex position, using an end-point position from one cluster and sliding fit to a second cluster
     *
     *  @param  position1 an end-point position for the first cluster
     *  @param  hitType1 the hit type of the first cluster
     *  @param  fitResult2 the two dimensional sliding fit result for the second cluster
     */
    void CreateEndpointVertex(const pandora::CartesianVector &position1, const pandora::HitType hitType1, const TwoDSlidingFitResult &fitResult2) const;

    /**
     *  @brief  Extrapolate 2D clusters, find where they cross, and match crossing points between views to create vertex candidates
     *
     *  @param  clusterVectorU the clusters in the u view
     *  @param  clusterVectorV the clusters in the v view
     *  @param  clusterVectorW the clusters in the w view
     */
    void CreateCrossingCandidates(const pandora::ClusterVector &clusterVectorU, const pandora::ClusterVector &clusterVectorV,
        const pandora::ClusterVector &clusterVectorW) const;

    /**
     *  @brief  Identify where (extrapolated) clusters plausibly cross in 2D
     *
     *  @param  clusterVector the input clusters
     *  @param  crossingPoints to receive the 2D crossing points
     */
    void FindCrossingPoints(const pandora::ClusterVector &clusterVector, pandora::CartesianPointVector &crossingPoints) const;

    /**
     *  @brief  Get a list of spacepoints representing cluster 2D hit positions and extrapolated positions
     *
     *  @param  pCluster address of the cluster
     *  @param  spacePoints to receive the list of spacepoints
     */
    void GetSpacepoints(const pandora::Cluster *const pCluster, pandora::CartesianPointVector &spacePoints) const;

    /**
     *  @brief  Identify where (extrapolated) clusters plausibly cross in 2D
     *
     *  @param  spacepoints1 space points for cluster 1
     *  @param  spacepoints2 space points for cluster 2
     *  @param  crossingPoints to receive the list of plausible 2D crossing points
     */
    void FindCrossingPoints(const pandora::CartesianPointVector &spacepoints1, const pandora::CartesianPointVector &spacepoints2,
        pandora::CartesianPointVector &crossingPoints) const;

    /**
     *  @brief  Attempt to create candidate vertex positions, using 2D crossing points in 2 views
     *
     *  @param  crossingPoints1 the crossing points in view 1
     *  @param  crossingPoints2 the crossing points in view 2
     *  @param  hitType1 the hit type of crossing points 1
     *  @param  hitType2 the hit type of crossing points 2
     *  @param  nCrossingCandidates to count the number of crossing candidates created
     */
    void CreateCrossingVertices(const pandora::CartesianPointVector &crossingPoints1, const pandora::CartesianPointVector &crossingPoints2,
        const pandora::HitType hitType1, const pandora::HitType hitType2, unsigned int &nCrossingCandidates) const;

    /**
     *  @brief  Add candidate vertices from any input vertices
     */
    void AddInputVertices() const;

    /**
     *  @brief  Creates a 2D sliding fit of a cluster and stores it for later use
     *
     *  @param  pCluster address of the relevant cluster
     */
    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get a sliding fit result from the algorithm cache
     *
     *  @param  pCluster address of the relevant cluster
     */
    const TwoDSlidingFitResult &GetCachedSlidingFitResult(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Clear relevant algorithm member variables between events
     */
    void TidyUp();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::Cluster *, pandora::CartesianPointVector> ClusterToSpacepointsMap;

    pandora::StringVector m_inputClusterListNames; ///< The list of cluster list names
    std::string m_inputVertexListName;             ///< The list name for existing candidate vertices
    std::string m_outputVertexListName;            ///< The name under which to save the output vertex list
    bool m_replaceCurrentVertexList;               ///< Whether to replace the current vertex list with the output list

    unsigned int m_slidingFitWindow;               ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap m_slidingFitResultMap; ///< The sliding fit result map

    unsigned int m_minClusterCaloHits; ///< The min number of hits in base cluster selection method
    float m_minClusterLengthSquared;   ///< The min length (squared) in base cluster selection method
    float m_chiSquaredCut;             ///< The chi squared cut (accept only 3D vertex positions with values below cut)

    bool m_enableEndpointCandidates; ///< Whether to create endpoint-based candidates
    float m_maxEndpointXDiscrepancy; ///< The max cluster endpoint discrepancy

    bool m_enableCrossingCandidates;          ///< Whether to create crossing vertex candidates
    unsigned int m_nMaxCrossingCandidates;    ///< The max number of crossing candidates to create
    float m_maxCrossingXDiscrepancy;          ///< The max cluster endpoint discrepancy
    unsigned int m_extrapolationNSteps;       ///< Number of extrapolation steps, at each end of cluster, of specified size
    float m_extrapolationStepSize;            ///< The extrapolation step size in cm
    float m_maxCrossingSeparationSquared;     ///< The separation (squared) between spacepoints below which a crossing can be identified
    float m_minNearbyCrossingDistanceSquared; ///< The minimum allowed distance between identified crossing positions

    bool m_reducedCandidates;           ///< Whether to reduce the number of candidates
    float m_selectionCutFactorMax;      ///< Maximum factor to multiply the base cluster selection cuts
    float m_nClustersPassingMaxCutsPar; ///< Parameter for number of clusters passing the max base cluster selection cuts
};

} // namespace lar_content

#endif // #ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
