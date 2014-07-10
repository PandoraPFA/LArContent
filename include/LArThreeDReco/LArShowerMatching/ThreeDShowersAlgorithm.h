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

    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

private:
    typedef std::map<unsigned int, pandora::CartesianVector> ShowerPositionMap;

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
     *  @param  slidingFitResultU the sliding fit result for the u view
     *  @param  slidingFitResultV the sliding fit result for the v view
     *  @param  slidingFitResultW the sliding fit result for the w view
     *  @param  minX the min x value of the common x-overlap range
     *  @param  maxX the max x value of the common x-overlap range
     *  @param  xPitch the x sampling pitch to be used
     *  @param  showerPositionMapU to receive the shower position map for the u view
     *  @param  showerPositionMapV to receive the shower position map for the v view
     *  @param  showerPositionMapW to receive the shower position map for the w view
     */
    void GetShowerPositionMaps(const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV,
        const TwoDSlidingFitResult &slidingFitResultW, const float minX, const float maxX, const float xPitch, ShowerPositionMap &showerPositionMapU,
        ShowerPositionMap &showerPositionMapV, ShowerPositionMap &showerPositionMapW) const;

    /**
     *  @brief  Get the fraction of hits, in the common x-overlap range, contained within the provided shower boundaries
     * 
     *  @param  pCluster the address of the candidate cluster
     *  @param  minX the min x value of the common x-overlap range
     *  @param  maxX the max x value of the common x-overlap range
     *  @param  xPitch the x sampling pitch used when constructing the shower edge position maps
     *  @param  edgeMap1 the first shower edge position map
     *  @param  edgeMap2 the second shower edge position map
     *  @param  nSampledHits to receive the number of hits in the common x-overlap range
     *  @param  nMatchedHits to receive the number of sampled hits contained within the shower edges
     */
    void GetHitOverlapFraction(const pandora::Cluster *const pCluster, const float minX, const float maxX, const float xPitch,
        const ShowerPositionMap &edgeMap1, const ShowerPositionMap &edgeMap2, unsigned int &nSampledHits, unsigned int &nMatchedHits) const;

    void ExamineTensor();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int                m_nMaxTensorToolRepeats;        ///< The maximum number of repeat loops over tensor tools
    unsigned int                m_slidingFitWindow;             ///< The layer window for the sliding linear fits

    unsigned int                m_minClusterCaloHits;           ///< The min number of hits in base cluster selection method
    float                       m_minClusterLengthSquared;      ///< The min length (squared) in base cluster selection method

    float                       m_minShowerMatchedFraction;     ///< The minimum shower matched sampling fraction to allow shower grouping
    unsigned int                m_minShowerMatchedPoints;       ///< The minimum number of matched shower sampling points to allow shower grouping

    typedef std::vector<ShowerTensorTool*> TensorToolList;
    TensorToolList              m_algorithmToolList;            ///< The algorithm tool list
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
