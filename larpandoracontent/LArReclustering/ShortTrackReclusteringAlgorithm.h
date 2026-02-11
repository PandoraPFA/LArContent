/**
 *  @file   larpandoracontent/LArReclustering/ShortTrackReclusteringAlgorithm.h
 *
 *  @brief  Tries to identify short tracks that are either lost, or under-clustered and recover missing hits.
 *          This algorithm looks at all three views of a PFO to try to identify if one or more views exhibit step changes in ADC values that
 *          might be indicative of a decay or intelastic interaction. If such a change occurs it looks for consistency across views and also
 *          examines nearby unclustered hits, or hits in nearby clusters to see if a more cohenrent clustering can be identified.
 *          If these various conditions are met, then reclustering is performed.
 *
 *  $Log: $
 */

#ifndef LAR_SHORT_TRACK_RECLUSTERING_ALGORITHM_H
#define LAR_SHORT_TRACK_RECLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArReclustering/ThreeDReclusteringFigureOfMeritBaseTool.h"

namespace lar_content
{

/**
  *  @brief  ShortTrackReclusteringAlgorithm class
  */
class ShortTrackReclusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ShortTrackReclusteringAlgorithm();

    /**
    *  @brief  Default destructor
    */
    ~ShortTrackReclusteringAlgorithm() = default;

private:
    typedef std::unordered_map<pandora::HitType, pandora::CaloHitSet> ViewToHitsMap;
    typedef std::unordered_map<pandora::HitType, pandora::ClusterList> ViewToClustersMap;
    typedef std::unordered_map<const pandora::Cluster *, const pandora::Pfo *> ClusterToPfoMap;
    typedef std::unordered_map<const pandora::Cluster *, pandora::FloatVector> ClusterToAdcMap;

    pandora::StatusCode Run();

    /**
     *  @brief  Helper function to get a list from Pandora content, and check that it is valid
     *
     *  @param  listName the name of the list to retrieve
     *  @param  pList the pointer to the list to be retrieved
     *
     *  @return true if the list is successfully retrieved and valid, false otherwise
     */
    template <typename T>
    bool GetList(const std::string &listName, const T *&pList) const;

    /**
     *  @brief  Collects the hits not currently assigned to a PFO, and maps them by view
     *
     *  @param  caloHitList the list of calo hits to consider
     *  @param  viewToUnclusteredHitsMap the map in which to store the unclustered hits, mapped by view
     */
    void CollectUnclusteredHits(const pandora::CaloHitList &caloHitList, ViewToHitsMap &viewToUnclusteredHitsMap) const;

    /**
     *  @brief  Collects the clusters currently assigned to PFOs, and maps them by view, and also maps clusters to their parent PFO
     *
     *  @param  pfoList the list of pfos to consider
     *  @param  viewToClustersMap the map in which to store the clusters, mapped by view
     *  @param  clusterToPfoMap the map in which to store the mapping of clusters to their parent PFO
     */
    void CollectClusters(const pandora::PfoList &pfoList, ViewToClustersMap &viewToClustersMap, ClusterToPfoMap &clusterToPfoMap) const;

    /**
     *  @brief  Gets the hits associated with a cluster, ordered relative to a specified vertex
     *
     *  @param  clusterHits the hits associated with the cluster for which to retrieve the ordered hits
     *  @param  vertex the vertex relative to which to order the hits
     *  @param  orderedHits the vector in which to store the ordered hits
     */
    void OrderHitsRelativeToVertex(const pandora::CaloHitVector &clusterHits, const LArPointingCluster::Vertex &vertex,
        pandora::CaloHitVector &orderedHits) const;

    /**
     *  @brief  Gets the moving average of the ADC values for ordered sets of hits, window is backward looking, so considers the "current"i
     *          hit and the previous (window-1) hits.
     *
     *  @param  adcs an ordered set of ADC values for which to calculate the moving average
     *  @param  movingAdc the vector in which to store the moving average of the ADC values for the hits
     *  @param  movingVariance the vector in which to store the moving variance of the ADC values for the hits
     *  @param  window the size of the window over which to calculate the moving average
     */
    void GetAdcMovingAverage(const pandora::FloatVector &adcs, pandora::FloatVector &movingAdc, pandora::FloatVector &movingVariance,
        const size_t window = 3) const;

    /**
     *  @brief  Gets the median value for a vector. The input vector need not be sorted.
     *
     *  @param  values a set of values for which to calculate the median
     *  @return the median value for the input set of values
     */
    template <typename T>
    double GetMedian(const std::vector<T> &values) const;

    /**
     *  @brief  Gets the indices of hits in an ordered set of hits for which there is a step change in ADC values. This function looks for step changes in
     *          ADC, but with relatively stable ADC either side of the discontinuity, indicative of a decay or inelastic interaction. This method attempts
     *          to filter out Bragg peaks or regions of high volatility.
     *
     *  @param  hits an ordered set of hits for which to identify stable ADC discontinuities
     *  @param  discontinuities the vector in which to store the indices of hits for which there is a stable ADC discontinuity
     *  @param  window the size of the window over which to calculate the moving average and identify step changes
     */
    void GetStableAdcDiscontinuities(const pandora::CaloHitVector &hits, pandora::IntVector &discontinuities, const size_t window = 3) const;

    /**
     *  @brief  Gets the ADC values for a set of hits and normalizes relative to the median ADC
     *
     *  @param  hits a set of hits for which to get the normalized ADC values
     *  @param  normalizedAdc the vector in which to store the normalized ADC values for the hits
     */
    void NormalizeAdc(const pandora::CaloHitVector &hits, pandora::FloatVector &normalizedAdc) const;

    /**
     *  @brief  Identifies whether a set of hits exhibits a Bragg peak. This method considers the consensus view from computation of the linear slope,
     *          quadratic curvature, contrast ratio and monotonicity.
     *
     *  @param  hits the set of hits for which to identify whether there is a Bragg peak
     *  @param  start the index of the first hit in the set of hits to consider when looking for a Bragg peak
     *  @param  end the index of the last hit in the set of hits to consider when looking for a Bragg peak
     *
     *  @return true if there is evidence for a Bragg peak in the specified range of hits, false otherwise
     */
    bool IsBraggPeak(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets a score indicative of the linear slope of the ADC values for a set of hits
     *
     *  @param  hits the set of hits for which to calculate the linear slope score
     *  @param  start the index of the first hit in the set of hits to consider when calculating the linear slope score
     *  @param  end the index of the last hit in the set of hits to consider when calculating the linear slope score
     *
     *  @return a score indicative of the linear slope of the ADC values for the specified range of hits
     */
    float GetLinearSlopeScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets a score indicative of the quadratic curvature of the ADC values for a set of hits
     *
     *  @param  hits the set of hits for which to calculate the quadratic curvature score
     *  @param  start the index of the first hit in the set of hits to consider when calculating the quadratic curvature score
     *  @param  end the index of the last hit in the set of hits to consider when calculating the quadratic curvature score
     *
     *  @return a score indicative of the quadratic curvature of the ADC values for the specified range of hits
     */
    float GetQuadraticCurvatureScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets a score indicative of the contrast in ADC values for a set of hits
     *
     *  @param  hits the set of hits for which to calculate the contrast score
     *  @param  start the index of the first hit in the set of hits to consider when calculating the contrast score
     *  @param  end the index of the last hit in the set of hits to consider when calculating the contrast score
     *
     *  @return a score indicative of the contrast in ADC values for the specified range of hits
     */
    float GetContrastScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    /**
     *  @brief  Gets a score indicative of the monotonicity of the ADC values for a set of hits
     *
     *  @param  hits the set of hits for which to calculate the monotonicity score
     *  @param  start the index of the first hit in the set of hits to consider when calculating the monotonicity score
     *  @param  end the index of the last hit in the set of hits to consider when calculating the monotonicity score
     *
     *  @return a score indicative of the monotonicity of the ADC values for the specified range of hits
     */
    float GetMonotonicityScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< Name of list of calo hits to consider during reclustering
    std::string m_pfoListName; ///< Name of list of track-like pfos to consider for reclustering
};

} // namespace lar_content

#endif // #ifndef LAR_SHORT_TRACK_RECLUSTERING_ALGORITHM_H
