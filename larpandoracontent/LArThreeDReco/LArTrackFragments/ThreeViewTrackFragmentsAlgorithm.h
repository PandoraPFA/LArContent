/**
 *  @file   larpandoracontent/LArThreeDReco/LArTrackFragments/ThreeViewTrackFragmentsAlgorithm.h
 *
 *  @brief  Header file for the three view fragments algorithm base class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_TRACK_FRAGMENTS_ALGORITHM_H
#define LAR_THREE_VIEW_TRACK_FRAGMENTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewTrackMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"

#include <unordered_map>

namespace lar_content
{

class FragmentTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewTrackFragmentsAlgorithm class
 */
class ThreeViewTrackFragmentsAlgorithm : public NViewTrackMatchingAlgorithm<ThreeViewMatchingControl<FragmentOverlapResult> >
{
public:
    typedef NViewTrackMatchingAlgorithm<ThreeViewMatchingControl<FragmentOverlapResult> > BaseAlgorithm;

    /**
     *  @brief  Default constructor
     */
    ThreeViewTrackFragmentsAlgorithm();

    void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);

    /**
     *  @brief  Rebuild clusters after fragmentation
     *
     *  @param rebuildList the list of clusters containing hits to be rebuilt
     *  @param newClusters the output list of clusters
     */
    void RebuildClusters(const pandora::ClusterList &rebuildList, pandora::ClusterList &newClusters) const;

protected:
    void PerformMainLoop();
    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

    /**
     *  @brief  Calculate overlap result for track fragment candidate consisting of two sliding fit results and a list of available clusters
     *
     *  @param  fitResult1 the first sliding fit result
     *  @param  fitResult2 the second sliding fit result
     *  @param  inputClusterList the input cluster list
     *  @param  pBestMatchedCluster to receive the address of the best matched cluster
     *  @param  fragmentOverlapResult to receive the populated fragment overlap result
     *
     *  @return statusCode, faster than throwing in regular use-cases
     */
    pandora::StatusCode CalculateOverlapResult(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
        const pandora::ClusterList &inputClusterList, const pandora::Cluster *&pBestMatchedCluster, FragmentOverlapResult &fragmentOverlapResult) const;

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;

    /**
     *  @brief  Get the list of projected positions, in the third view, corresponding to a pair of sliding fit results
     *
     *  @param  fitResult1 the first sliding fit result
     *  @param  fitResult2 the second sliding fit result
     *  @param  projectedPositions to receive the list of projected positions
     *
     *  @return statusCode, faster than throwing in regular use-cases
     */
    pandora::StatusCode GetProjectedPositions(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
        pandora::CartesianPointVector &projectedPositions) const;

    /**
     *  @brief  Get the list of hits associated with the projected positions and a useful hit to cluster map
     *
     *  @param  inputClusterList the input cluster list
     *  @param  projectedPositions the list of projected positions
     *  @param  hitToClusterMap to receive the hit to cluster map
     *  @param  matchedCaloHits to receive the list of associated calo hits
     *
     *  @return statusCode, faster than throwing in regular use-cases
     */
    pandora::StatusCode GetMatchedHits(const pandora::ClusterList &inputClusterList, const pandora::CartesianPointVector &projectedPositions,
        HitToClusterMap &hitToClusterMap, pandora::CaloHitList &matchedCaloHits) const;

    /**
     *  @brief  Get the list of the relevant clusters and the address of the single best matched cluster
     *
     *  @param  matchedHits the list of matched calo hits
     *  @param  hitToClusterMap the hit to cluster map
     *  @param  matchedClusters to receive the list of matched clusters
     *  @param  pBestMatchedCluster to receive the address of the single best matched cluster
     *
     *  @return statusCode, faster than throwing in regular use-cases
     */
    pandora::StatusCode GetMatchedClusters(const pandora::CaloHitList &matchedHits, const HitToClusterMap &hitToClusterMap,
        pandora::ClusterList &matchedClusters, const pandora::Cluster *&pBestMatchedCluster) const;

    /**
     *  @brief  Get the populated fragment overlap result
     *
     *  @param  projectedPositions the list of projected positions
     *  @param  matchedHits the list of matched hits
     *  @param  matchedClusters the list of matched clusters
     *  @param  fragmentOverlapResult to receive the populated fragment overlap result
     */
    void GetFragmentOverlapResult(const pandora::CartesianPointVector &projectedPositions, const pandora::CaloHitList &matchedHits,
        const pandora::ClusterList &matchedClusters, FragmentOverlapResult &fragmentOverlapResult) const;

    /**
     *  @brief  Whether the matched clusters are consistent with the projected positions
     *
     *  @param  projectedPositions the list of projected positions
     *  @param  matchedClusters the list of matched clusters
     *
     *  @return boolean
     */
    bool CheckMatchedClusters(const pandora::CartesianPointVector &projectedPositions, const pandora::ClusterList &matchedClusters) const;

    /**
     *  @brief  Whether the matched clusters and hits pass the algorithm quality cuts
     *
     *  @param  fragmentOverlapResult the fragment overlap result
     *
     *  @return boolean
     */
    bool CheckOverlapResult(const FragmentOverlapResult &overlapResult) const;

    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::Cluster*, unsigned int> ClusterToMatchedHitsMap;

    std::string         m_reclusteringAlgorithmName;        ///< Name of daughter algorithm to use for cluster re-building

    typedef std::vector<FragmentTensorTool*> TensorToolVector;
    TensorToolVector    m_algorithmToolVector;              ///< The algorithm tool list

    unsigned int        m_nMaxTensorToolRepeats;            ///< The maximum number of repeat loops over tensor tools

    float               m_minXOverlap;                      ///< requirement on minimum X overlap for associated clusters
    float               m_minXOverlapFraction;              ///< requirement on minimum X overlap fraction for associated clusters
    float               m_maxPointDisplacementSquared;      ///< maximum allowed distance (squared) between projected points and associated hits
    float               m_minMatchedSamplingPointFraction;  ///< minimum fraction of matched sampling points
    unsigned int        m_minMatchedHits;                   ///< minimum number of matched calo hits
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  FragmentTensorTool class
 */
class FragmentTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeViewTrackFragmentsAlgorithm::MatchingType::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_TRACK_FRAGMENTS_ALGORITHM_H
