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
     *  @brief  Sort tensor elements by number of matched sampling points, using matched fraction then xoverlap span to resolve ties
     *
     *  @param  lhs the first tensor element
     *  @param  rhs the second tensor element
     *
     *  @return boolean
     */
    static bool SortByNMatchedSamplingPoints(const TensorType::Element &lhs, const TensorType::Element &rhs);

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

        float       m_uMinX;         ///< The min x value in the u view
        float       m_uMaxX;         ///< The max x value in the u view
        float       m_vMinX;         ///< The min x value in the v view
        float       m_vMaxX;         ///< The max x value in the v view
        float       m_wMinX;         ///< The min x value in the w view
        float       m_wMaxX;         ///< The max x value in the w view
        float       m_minX;          ///< The min x value of the common x-overlap range
        float       m_maxX;          ///< The max x value of the common x-overlap range
        float       m_xOverlapSpan;  ///< The x-overlap span
        float       m_xPitch;        ///< The x sampling pitch to be used
    };

    /**
     *  @brief  ShowerEdge
     */
    class ShowerEdge
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  xCoordinate the x coordinate
         *  @param  highEdgeZ the shower high edge z coordinate
         *  @param  lowEdgeZ the shower low edge z coordinate
         */
        ShowerEdge(const float xxCoordinate, const float highEdgeZ, const float lowEdgeZ);

        /**
         *  @param  Get the x coordinate
         * 
         *  @return the x coordinate
         */
        float GetXCoordinate() const;

        /**
         *  @param  Get the shower high edge z coordinate
         * 
         *  @return the shower high edge z coordinate
         */
        float GetHighEdgeZ() const;

        /**
         *  @param  Get the shower low edge z coordinate
         * 
         *  @return the shower low edge z coordinate
         */
        float GetLowEdgeZ() const;

    private:
        float       m_xCoordinate;   ///< The x coordinate
        float       m_highEdgeZ;     ///< The shower high edge z coordinate
        float       m_lowEdgeZ;      ///< The shower low edge z coordinate
    };

    typedef std::map<unsigned int, ShowerEdge> ShowerPositionMap;
    typedef std::pair<ShowerPositionMap, ShowerPositionMap> ShowerPositionMapPair;

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
    void GetShowerPositionMaps(const SlidingShowerFitResult &fitResultU, const SlidingShowerFitResult &fitResultV, const SlidingShowerFitResult &fitResultW,
        const XSampling &xSampling, ShowerPositionMapPair &positionMapsU, ShowerPositionMapPair &positionMapsV, ShowerPositionMapPair &positionMapsW) const;

    /**
     *  @brief  Get the most appropriate shower edges at a given x coordinate
     * 
     *  @param  x the x coordinate
     *  @param  fitResult the shower fit result
     *  @param  floatVector to receive the list of intersections of the shower fit at the given x coordinate
     */
    void GetShowerEdges(const float x, const SlidingShowerFitResult &fitResult, pandora::FloatVector &floatVector) const;

    /**
     *  @brief  Get the best fraction of hits, in the common x-overlap range, contained within the provided pair of shower boundaries
     * 
     *  @param  pCluster the address of the candidate cluster
     *  @param  xSampling the x sampling details
     *  @param  positionMaps the shower edge position maps
     *  @param  nSampledHits to receive the number of hits in the common x-overlap range
     *  @param  nMatchedHits to receive the number of sampled hits contained within the shower edges
     */
    void GetBestHitOverlapFraction(const pandora::Cluster *const pCluster, const XSampling &xSampling, const ShowerPositionMapPair &positionMaps,
        unsigned int &nSampledHits, unsigned int &nMatchedHits) const;

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

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDShowersAlgorithm::ShowerEdge::ShowerEdge(const float xCoordinate, const float edge1, const float edge2) :
    m_xCoordinate(xCoordinate),
    m_highEdgeZ(std::max(edge1, edge2)),
    m_lowEdgeZ(std::min(edge1, edge2))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDShowersAlgorithm::ShowerEdge::GetXCoordinate() const
{
    return m_xCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDShowersAlgorithm::ShowerEdge::GetHighEdgeZ() const
{
    return m_highEdgeZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDShowersAlgorithm::ShowerEdge::GetLowEdgeZ() const
{
    return m_lowEdgeZ;
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_SHOWERS_ALGORITHM_H
