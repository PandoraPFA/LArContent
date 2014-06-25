/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.h
 *
 *  @brief  Header file for the three dimensional fragments algorithm base class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRACK_FRAGMENTS_ALGORITHM_H
#define LAR_THREE_D_TRACK_FRAGMENTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArObjects/LArTrackOverlapResult.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDTracksBaseAlgorithm.h"

namespace lar
{

class FragmentTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDTrackFragmentsAlgorithm class
 */
class ThreeDTrackFragmentsAlgorithm : public ThreeDTracksBaseAlgorithm<FragmentOverlapResult>
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

    void UpdateForNewCluster(pandora::Cluster *const pNewCluster);

    /**
     *  @brief  Rebuild clusters after fragmentation
     *
     *  @param pCluster the cluster to be rebuilt
     *  @param newClusters the output list of clusters
     */
    void RebuildClusters(pandora::Cluster *pCluster, pandora::ClusterList &newClusters) const;

protected:
    void PerformMainLoop();
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

    /**
     *  @brief  Calculate overlap result for track fragment candidate consisting of two sliding fit results and a list of available clusters
     *
     *  @param  fitResult1 the first sliding fit result
     *  @param  fitResult2 the second sliding fit result
     *  @param  inputClusterList the input cluster list
     *  @param  pBestMatchedCluster to receive the address of the best matched cluster
     *  @param  fragmentOverlapResult to receive the populated fragment overlap result
     */
    void CalculateOverlapResult(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
        const pandora::ClusterList &inputClusterList, pandora::Cluster *&pBestMatchedCluster, FragmentOverlapResult &fragmentOverlapResult) const;

    typedef std::map<const pandora::CaloHit*, pandora::Cluster*> HitToClusterMap;
    typedef std::map<const pandora::CaloHit*, pandora::CaloHitList> HitToHitMap;

    /**
     *  @brief  Get the list of projected positions, in the third view, corresponding to a pair of sliding fit results
     *
     *  @param  fitResult1 the first sliding fit result
     *  @param  fitResult2 the second sliding fit result
     *  @param  projectedPositions to receive the list of projected positions
     */
    void GetProjectedPositions(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
        pandora::CartesianPointList &projectedPositions) const;

    /**
     *  @brief  Get the list of hits associated with the projected positions and a useful hit to cluster map
     *
     *  @param  inputClusterList the input cluster list
     *  @param  projectedPositions the list of projected positions
     *  @param  hitToClusterMap to receive the hit to cluster map
     *  @param  matchedCaloHits to receive the list of associated calo hits
     */
    void GetMatchedHits(const pandora::ClusterList &inputClusterList, const pandora::CartesianPointList &projectedPositions,
        HitToClusterMap &hitToClusterMap, pandora::CaloHitList &matchedCaloHits) const;

    /**
     *  @brief  Get the list of the relevant clusters and the address of the single best matched cluster
     *
     *  @param  matchedHits the list of matched calo hits
     *  @param  hitToClusterMap the hit to cluster map
     *  @param  matchedClusters to receive the list of matched clusters
     *  @param  pBestMatchedCluster to receive the address of the single best matched cluster
     */
    void GetMatchedClusters(const pandora::CaloHitList &matchedHits, const HitToClusterMap &hitToClusterMap,
        pandora::ClusterList &matchedClusters, pandora::Cluster *&pBestMatchedCluster) const;

    /**
     *  @brief  Get the populated fragment overlap result
     *
     *  @param  projectedPositions the list of projected positions
     *  @param  matchedHits the list of matched hits
     *  @param  matchedClusters the list of matched clusters
     *  @param  fragmentOverlapResult to receive the populated fragment overlap result
     */
    void GetFragmentOverlapResult(const pandora::CartesianPointList &projectedPositions, const pandora::CaloHitList &matchedHits,
        const pandora::ClusterList &matchedClusters, FragmentOverlapResult &fragmentOverlapResult) const;

    /**
     *  @brief  Whether the matched clusters are consistent with the projected positions
     *
     *  @param  projectedPositions the list of projected positions
     *  @param  matchedClusters the list of matched clusters
     *
     *  @return boolean
     */
    bool CheckMatchedClusters(const pandora::CartesianPointList &projectedPositions, const pandora::ClusterList &matchedClusters) const;

    /**
     *  @brief  Whether the matched clusters and hits pass the algorithm quality cuts
     *
     *  @param  fragmentOverlapResult the fragment overlap result
     *
     *  @return boolean
     */
    bool CheckOverlapResult(const FragmentOverlapResult &overlapResult) const;

    void ExamineTensor();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<pandora::Cluster*, unsigned int> ClusterToMatchedHitsMap;
    typedef std::vector<FragmentTensorTool*> TensorToolList;

    std::string         m_inputCaloHitListNameU;            ///< Name of the view U calo hit list
    std::string         m_inputCaloHitListNameV;            ///< Name of the view V calo hit list
    std::string         m_inputCaloHitListNameW;            ///< Name of the view W calo hit list
    std::string         m_reclusteringAlgorithmName;        ///< Name of daughter algorithm to use for cluster re-building

    unsigned int        m_nMaxTensorToolRepeats;            ///< The maximum number of repeat loops over tensor tools
    TensorToolList      m_algorithmToolList;                ///< The algorithm tool list

    float               m_minXOverlap;                      ///< requirement on minimum X overlap for associated clusters
    float               m_minXOverlapFraction;              ///< requirement on minimum X overlap fraction for associated clusters
    float               m_maxPointDisplacement;             ///< maximum allowed distance between projected points and associated hits
    float               m_maxPointDisplacementSquared;      ///< maximum allowed distance (squared) between projected points and associated hits
    unsigned int        m_minMatchedSamplingPoints;         ///< minimum number of matched sampling points
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
    typedef ThreeDTrackFragmentsAlgorithm::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeDTrackFragmentsAlgorithm *pAlgorithm, TensorType &overlapTensor) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTrackFragmentsAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDTrackFragmentsAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRACK_FRAGMENTS_ALGORITHM_H
