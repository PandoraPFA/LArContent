/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/MissingTrackSegmentTool.h
 * 
 *  @brief  Header file for the missing track segment tool class.
 * 
 *  $Log: $
 */
#ifndef MISSING_TRACK_SEGMENT_TOOL_H
#define MISSING_TRACK_SEGMENT_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  MissingTrackSegmentTool class
 */
class MissingTrackSegmentTool : public TensorManipulationTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

private:
    /**
     *  @brief  Particle class
     */
    class Particle
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  element the tensor element
         */
        Particle(const TensorType::Element &element);

        pandora::Cluster   *m_pShortCluster;            ///< Address of the short cluster
        pandora::Cluster   *m_pCluster1;                ///< Address of long cluster in view 1
        pandora::Cluster   *m_pCluster2;                ///< Address of long cluster in view 2
        pandora::HitType    m_shortHitType;             ///< The hit type of the short cluster
        pandora::HitType    m_hitType1;                 ///< The hit type of the long cluster in view 1
        pandora::HitType    m_hitType2;                 ///< The hit type of the long cluster in view 2
        float               m_shortMinX;                ///< The min x coordinate of the short cluster
        float               m_shortMaxX;                ///< The max x coordinate of the short cluster
        float               m_longMinX;                 ///< The min x coordinate of the long clusters
        float               m_longMaxX;                 ///< The max x coordinate of the long clusters
    };

    /**
     *  @brief  SegmentOverlap class
     */
    class SegmentOverlap
    {
    public:
        /**
         *  @brief  Default constructor
         */
        SegmentOverlap();

        unsigned int        m_nSamplingPoints;          ///< The number of sampling points
        unsigned int        m_nMatchedSamplingPoints;   ///< The number of matched sampling points
        float               m_pseudoChi2Sum;            ///< The pseudo chi2 sum
    };

    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<pandora::Cluster*, TwoDSlidingFitResult> SlidingFitResultMap;
    typedef std::map<pandora::Cluster*, SegmentOverlap> SegmentOverlapMap;
    typedef std::map<pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;

    /**
     *  @brief  Find remaining tracks, hidden by missing track segments (and maybe other ambiguities) in the tensor
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     *  @param  clusterMergeMap to receive the cluster merge map
     */
    void FindTracks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector,
        ClusterMergeMap &clusterMergeMap) const;

    typedef ThreeDTransverseTracksAlgorithm::IteratorList IteratorList;

    /**
     *  @brief  Select a list of the relevant elements from a set of connected tensor elements
     * 
     *  @param  elementList the full list of connected tensor elements
     *  @param  usedClusters the list of clusters already marked as to be added to a pfo
     *  @param  iteratorList to receive a list of iterators to long track-like elements
     */
    void SelectElements(const TensorType::ElementList &elementList, const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Whether a provided tensor element can be used to construct a pfo
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  element the tensor element
     *  @param  overlapTensor the overlap tensor
     *  @param  usedClusters the list of used clusters
     *  @param  clusterMergeMap to receive the cluster merge map
     */
    bool PassesParticleChecks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::Element &element, const TensorType &overlapTensor,
        pandora::ClusterList &usedClusters, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Get a list of candidate clusters, which may represent missing track segments for a provided particle
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  particle the particle
     *  @param  overlapTensor the overlap tensor
     *  @param  usedClusters the list of used clusters
     *  @param  candidateClusters to receive the list of candidate clusters
     */
    void GetCandidateClusters(ThreeDTransverseTracksAlgorithm *pAlgorithm, const Particle &particle, const TensorType &overlapTensor,
        const pandora::ClusterList &usedClusters, pandora::ClusterList &candidateClusters) const;

    /**
     *  @brief  Get a sliding fit result map for the list of candidate clusters
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  candidateClusters the list of candidate clusters
     *  @param  slidingFitResultMap to receive the sliding fit result map
     */
    void GetSlidingFitResultMap(ThreeDTransverseTracksAlgorithm *pAlgorithm, const pandora::ClusterList &candidateClusterList,
        SlidingFitResultMap &slidingFitResultMap) const;

    /**
     *  @brief  Get a segment overlap map, describing overlap between a provided particle and all clusters in a sliding fit result map
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  particle the particle
     *  @param  slidingFitResultMap the sliding fit result map
     *  @param  segmentOverlapMap to receive the segment overlap map
     */
    void GetSegmentOverlapMap(ThreeDTransverseTracksAlgorithm *pAlgorithm, const Particle &particle,
        const SlidingFitResultMap &slidingFitResultMap, SegmentOverlapMap &segmentOverlapMap) const;

    /**
     *  @brief  Make decisions about whether to create a pfo for a provided particle and whether to make cluster merges
     * 
     *  @param  particle the particle
     *  @param  slidingFitResultMap the sliding fit result map
     *  @param  segmentOverlapMap the segment overlap map
     *  @param  usedClusters the list of used clusters
     *  @param  clusterMergeMap to receive details of cluster merges clusterMergeMap
     * 
     *  @return whether to make the particle
     */
    bool MakeDecisions(const Particle &particle, const SlidingFitResultMap &slidingFitResultMap, const SegmentOverlapMap &segmentOverlapMap,
        pandora::ClusterList &usedClusters, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Whether the pointing information is consistent for the short cluster in the tensor and the new cluster, picked-up
     *          via sampling of the two long clusters in the tensor
     * 
     *  @param  pClusterA the address of cluster A
     *  @param  pClusterB the address of cluster B
     * 
     *  @return whether the directions are consistent
     */
    bool AreDirectionsConsistent(pandora::Cluster *const pClusterA, pandora::Cluster *const pClusterB) const;

    /**
     *  @brief  Perform cluster merges described in a provided cluster merge map
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  clusterMergeMap the cluster merge map
     */
    void PerformClusterMerges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ClusterMergeMap &clusterMergeMap) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (between long clusters and short cluster vs. shared overlap)
    unsigned int    m_minMatchedSamplingPointRatio;     ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution

    float           m_minParticleXEndSeparation;        ///< Distance used to ensure that undershoot (due to missing segment) is at one end of x span
    unsigned int    m_minCaloHitsInCandidateCluster;    ///< The min no. of calo hits in a candidate cluster, for matching with long clusters
    float           m_pseudoChi2Cut;                    ///< The pseudo chi2 cut to determine whether a sampling point is matched

    unsigned int    m_makePfoMinSamplingPoints;         ///< The min number of sampling points in order to be able to make pfo
    unsigned int    m_makePfoMinMatchedSamplingPoints;  ///< The min number of matched sampling points in order to be able to make pfo
    float           m_makePfoMinMatchedFraction;        ///< The min matched sampling point fraction in order to be able to make pfo
    float           m_makePfoMaxImpactParameter;        ///< The max transverse impact parameter in order to be able to make pfo

    float           m_mergeMaxChi2PerSamplingPoint;     ///< The max value of chi2 per sampling point in order to merge cluster with parent
    float           m_mergeXContainmentTolerance;       ///< The tolerance in determining whether candidate cluster is contained in x window
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *MissingTrackSegmentTool::Factory::CreateAlgorithmTool() const
{
    return new MissingTrackSegmentTool();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline MissingTrackSegmentTool::SegmentOverlap::SegmentOverlap() :
    m_nSamplingPoints(0),
    m_nMatchedSamplingPoints(0),
    m_pseudoChi2Sum(0.f)
{
}

} // namespace lar

#endif // #ifndef MISSING_TRACK_SEGMENT_TOOL_H
