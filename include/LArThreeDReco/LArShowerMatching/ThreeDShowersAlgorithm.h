/**
 *  @file   LArContent/include/LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional showers algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_SHOWERS_ALGORITHM_H
#define LAR_THREE_D_SHOWERS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArShowerOverlapResult.h"
#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDBaseAlgorithm.h"

namespace lar
{

class ShowerTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDShowersAlgorithm class
 */
class ThreeDShowersAlgorithm : public ThreeDBaseAlgorithm<ShowerOverlapResult>
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
     *  @brief  SlidingShowerFitResult class
     */
    class SlidingShowerFitResult
    {
    public:
        TwoDSlidingFitResult    m_showerFitResult;              ///< The sliding fit result for the full shower cluster
        TwoDSlidingFitResult    m_negativeEdgeFitResult;        ///< The sliding fit result for the negative shower edge
        TwoDSlidingFitResult    m_positiveEdgeFitResult;        ///< The sliding fit result for the positive shower edge
    };

    typedef std::map<pandora::Cluster*, SlidingShowerFitResult> SlidingShowerFitResultMap;

    /**
     *  @brief  Get a sliding shower fit result from the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    const SlidingShowerFitResult &GetCachedSlidingFitResult(pandora::Cluster *const pCluster) const;

    void UpdateForNewCluster(pandora::Cluster *const pNewCluster);
    void UpdateUponDeletion(pandora::Cluster *const pDeletedCluster);
    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

private:
    /**
     *  @brief  XOverlap class
     */
    class XOverlap
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  fitResultU the sliding fit result for the u view
         *  @param  fitResultV the sliding fit result for the v view
         *  @param  fitResultW the sliding fit result for the w view
         */
        XOverlap(const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV, const TwoDSlidingFitResult &fitResultW);

        float      m_minX;          ///< The min x value of the common x-overlap range
        float      m_maxX;          ///< The max x value of the common x-overlap range
        float      m_xOverlap;      ///< The x-overlap
        float      m_xPitch;        ///< The x sampling pitch to be used
    };

    void PreparationStep();
    void TidyUp();

    /**
     *  @brief  Add a new sliding fit result, for the specified cluster, to the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    void AddToSlidingFitCache(pandora::Cluster *const pCluster);

    /**
     *  @brief  Remova an existing sliding fit result, for the specified cluster, from the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    void RemoveFromSlidingFitCache(pandora::Cluster *const pCluster);

    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

    /**
     *  @brief  Calculate the overlap result for given group of clusters
     *
     *  @param  pClusterU the cluster from the U view
     *  @param  pClusterV the cluster from the V view
     *  @param  pClusterW the cluster from the W view
     *  @param  overlapResult to receive the overlap result
     */
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW,
        ShowerOverlapResult &overlapResult);

    typedef std::map<unsigned int, pandora::CartesianVector> ShowerPositionMap;

    /**
     *  @brief  Get the shower position maps
     * 
     *  @param  fitResultU the sliding fit result for the u view
     *  @param  fitResultV the sliding fit result for the v view
     *  @param  fitResultW the sliding fit result for the w view
     *  @param  xOverlap the common x-overlap details
     *  @param  positionMapU to receive the shower position map for the u view
     *  @param  positionMapV to receive the shower position map for the v view
     *  @param  positionMapW to receive the shower position map for the w view
     */
    void GetShowerPositionMaps(const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV, const TwoDSlidingFitResult &fitResultW,
        const XOverlap &xOverlap, ShowerPositionMap &positionMapU, ShowerPositionMap &positionMapV, ShowerPositionMap &positionMapW) const;

    /**
     *  @brief  Get the fraction of hits, in the common x-overlap range, contained within the provided shower boundaries
     * 
     *  @param  pCluster the address of the candidate cluster
     *  @param  xOverlap the common x-overlap details
     *  @param  positionMap1 the first shower edge position map
     *  @param  positionMap2 the second shower edge position map
     *  @param  nSampledHits to receive the number of hits in the common x-overlap range
     *  @param  nMatchedHits to receive the number of sampled hits contained within the shower edges
     */
    void GetHitOverlapFraction(const pandora::Cluster *const pCluster, const XOverlap &xOverlap, const ShowerPositionMap &positionMap1,
        const ShowerPositionMap &positionMap2, unsigned int &nSampledHits, unsigned int &nMatchedHits) const;

    void ExamineTensor();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<ShowerTensorTool*> TensorToolList;
    TensorToolList              m_algorithmToolList;            ///< The algorithm tool list
    unsigned int                m_nMaxTensorToolRepeats;        ///< The maximum number of repeat loops over tensor tools

    unsigned int                m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    SlidingShowerFitResultMap   m_slidingFitResultMap;          ///< The sliding shower fit result map

    bool                        m_ignoreUnavailableClusters;    ///< Whether to ignore (skip-over) unavailable clusters
    unsigned int                m_minClusterCaloHits;           ///< The min number of hits in base cluster selection method
    float                       m_minClusterLengthSquared;      ///< The min length (squared) in base cluster selection method

    float                       m_minShowerMatchedFraction;     ///< The minimum shower matched sampling fraction to allow shower grouping
    unsigned int                m_minShowerMatchedPoints;       ///< The minimum number of matched shower sampling points to allow shower grouping
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ShowerTensorTool class
 */
class ShowerTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeDShowersAlgorithm::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeDShowersAlgorithm *pAlgorithm, TensorType &overlapTensor) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDShowersAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDShowersAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_SHOWERS_ALGORITHM_H
