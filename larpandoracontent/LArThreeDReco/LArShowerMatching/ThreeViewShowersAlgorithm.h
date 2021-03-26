/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/ThreeViewShowersAlgorithm.h
 *
 *  @brief  Header file for the three view showers algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_SHOWERS_ALGORITHM_H
#define LAR_THREE_VIEW_SHOWERS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArShowerOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"

namespace lar_content
{

class ShowerTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewShowersAlgorithm class
 */
class ThreeViewShowersAlgorithm : public NViewMatchingAlgorithm<ThreeViewMatchingControl<ShowerOverlapResult>>
{
public:
    typedef NViewMatchingAlgorithm<ThreeViewMatchingControl<ShowerOverlapResult>> BaseAlgorithm;

    /**
     *  @brief  Default constructor
     */
    ThreeViewShowersAlgorithm();

    /**
     *  @brief  Get a sliding shower fit result from the algorithm cache
     *
     *  @param  pCluster address of the relevant cluster
     */
    const TwoDSlidingShowerFitResult &GetCachedSlidingFitResult(const pandora::Cluster *const pCluster) const;

    void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;
    void PrepareInputClusters(pandora::ClusterList &preparedClusterList);

private:
    /**
     *  @brief  XSampling class
     */
    class XSampling
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  fitResultU the sliding fit result for the u view
         *  @param  fitResultV the sliding fit result for the v view
         *  @param  fitResultW the sliding fit result for the w view
         */
        XSampling(const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV, const TwoDSlidingFitResult &fitResultW);

        /**
         *  @brief  Convert an x position into a sampling bin
         *
         *  @param  x the input x coordinate
         *  @param  xBin to receive the x bin
         *
         *  @return status code, faster than throwing in regular use-cases
         */
        pandora::StatusCode GetBin(const float x, int &xBin) const;

        float m_uMinX;        ///< The min x value in the u view
        float m_uMaxX;        ///< The max x value in the u view
        float m_vMinX;        ///< The min x value in the v view
        float m_vMaxX;        ///< The max x value in the v view
        float m_wMinX;        ///< The min x value in the w view
        float m_wMaxX;        ///< The max x value in the w view
        float m_minX;         ///< The min x value of the common x-overlap range
        float m_maxX;         ///< The max x value of the common x-overlap range
        float m_xOverlapSpan; ///< The x-overlap span
        float m_nPoints;      ///< The number of sampling points to be used
    };

    void TidyUp();

    /**
     *  @brief  Add a new sliding fit result, for the specified cluster, to the algorithm cache
     *
     *  @param  pCluster address of the relevant cluster
     */
    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Remova an existing sliding fit result, for the specified cluster, from the algorithm cache
     *
     *  @param  pCluster address of the relevant cluster
     */
    void RemoveFromSlidingFitCache(const pandora::Cluster *const pCluster);

    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

    /**
     *  @brief  Calculate the overlap result for given group of clusters
     *
     *  @param  pClusterU the cluster from the U view
     *  @param  pClusterV the cluster from the V view
     *  @param  pClusterW the cluster from the W view
     *  @param  overlapResult to receive the overlap result
     */
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
        const pandora::Cluster *const pClusterW, ShowerOverlapResult &overlapResult);

    typedef std::pair<ShowerPositionMap, ShowerPositionMap> ShowerPositionMapPair;

    /**
     *  @brief  Get the shower position maps
     *
     *  @param  fitResultU the sliding shower fit result for the u view
     *  @param  fitResultV the sliding shower fit result for the v view
     *  @param  fitResultW the sliding shower fit result for the w view
     *  @param  xSampling the x sampling details
     *  @param  positionMapsU to receive the shower position maps for the u view
     *  @param  positionMapsV to receive the shower position maps for the v view
     *  @param  positionMapsW to receive the shower position maps for the w view
     */
    void GetShowerPositionMaps(const TwoDSlidingShowerFitResult &fitResultU, const TwoDSlidingShowerFitResult &fitResultV,
        const TwoDSlidingShowerFitResult &fitResultW, const XSampling &xSampling, ShowerPositionMapPair &positionMapsU,
        ShowerPositionMapPair &positionMapsV, ShowerPositionMapPair &positionMapsW) const;

    /**
     *  @brief  Get the best fraction of hits, in the common x-overlap range, contained within the provided pair of shower boundaries
     *
     *  @param  pCluster the address of the candidate cluster
     *  @param  xSampling the x sampling details
     *  @param  positionMaps the shower edge position maps
     *  @param  nSampledHits to receive the number of hits in the common x-overlap range
     *  @param  nMatchedHits to receive the number of sampled hits contained within the shower edges
     */
    void GetBestHitOverlapFraction(const pandora::Cluster *const pCluster, const XSampling &xSampling,
        const ShowerPositionMapPair &positionMaps, unsigned int &nSampledHits, unsigned int &nMatchedHits) const;

    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<ShowerTensorTool *> TensorToolVector;
    TensorToolVector m_algorithmToolVector; ///< The algorithm tool vector
    unsigned int m_nMaxTensorToolRepeats;   ///< The maximum number of repeat loops over tensor tools

    unsigned int m_slidingFitWindow;                     ///< The layer window for the sliding linear fits
    TwoDSlidingShowerFitResultMap m_slidingFitResultMap; ///< The sliding shower fit result map

    bool m_ignoreUnavailableClusters;  ///< Whether to ignore (skip-over) unavailable clusters
    unsigned int m_minClusterCaloHits; ///< The min number of hits in base cluster selection method
    float m_minClusterLengthSquared;   ///< The min length (squared) in base cluster selection method

    float m_minShowerMatchedFraction;      ///< The minimum shower matched sampling fraction to allow shower grouping
    unsigned int m_minShowerMatchedPoints; ///< The minimum number of matched shower sampling points to allow shower grouping
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ShowerTensorTool class
 */
class ShowerTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeViewShowersAlgorithm::MatchingType::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeViewShowersAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_SHOWERS_ALGORITHM_H
