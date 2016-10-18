/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TracksCrossingGapsTool.h
 * 
 *  @brief  Header file for the long tracks tool class.
 * 
 *  $Log: $
 */
#ifndef TRACKS_CROSSING_GAPS_TOOL_H
#define TRACKS_CROSSING_GAPS_TOOL_H 1

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TracksCrossingGapsTool class
 */
class TracksCrossingGapsTool : public TransverseTensorTool
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

    /**
     *  @brief  Default constructor
     */
    TracksCrossingGapsTool();


    bool Run(ThreeDTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find tracks crossing gaps, with unambiguous connection but poor overlap due to gaps
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindTracks(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Select a list of track-like elements crossing a gap in one or more views from a set of connected tensor elements
     * 
     *  @param  elementList the full list of connected tensor elements
     *  @param  usedClusters the set of clusters already marked as to be added to a pfo
     *  @param  iteratorList to receive a list of iterators to long track-like elements
     */
    void SelectElements(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, 
			const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const;

    /**
     * @brief Calculate the effective overlap fractions given a set of clusters, taking gaps into account
     *
     * @param element the connected tensor element
     * @param xOverlapFractionU,V,W to receive the effective overlap fraction in each view
     */
    pandora::StatusCode CalculateEffectiveOverlapFractions(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const TensorType::Element &element, float &xOverlapFractionU, float &xOverlapFractionV, float &xOverlapFractionW) const;

    /**                                                                                                                                                                
     * @brief Calculate the effective overlap span given a set of clusters, taking gaps into account                                                                
     *                                                                                                                                                                
     * @param pClusterU, pClusterV, pClusterW the three clusters for which to calculate the overlap result                                                     
     * @param xMin/MaxEffU,V,W to receive the effective min and max in each view
     */                                     
    pandora::StatusCode CalculateEffectiveOverlapSpan(ThreeDTransverseTracksAlgorithm *const pAlgorithm,const pandora::Cluster *const pClusterU, float &xMinEffU, float &xMaxEffU, const pandora::Cluster *const pClusterV, float &xMinEffV, float &xMaxEffV, const pandora::Cluster *const pClusterW, float &xMinEffW, float &xMaxEffW) const;     

    /**                                                                                                                                                                
     * @brief Check whether there is any gap in the three U-V-W clusters combination
     *
     * @param xSample, the x coordinate we are checking
     * @param pClusterU, pClusterV, pClusterW the three clusters we are investigating
     * @param gapInU,V,W, booleans to receive whether there is a gap in each view 
     */
    bool PassesGapsChecks(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const float &xSample, const pandora::Cluster *const pClusterU,const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW, bool &gapInU, bool &gapInV, bool &gapInW) const;

    /**                                                                                                                                                                     * @brief Check individually each cluster where a gap might be present
     *
     * @param xSample, the x coordinate we are checking 
     * @param pClusterGap, slidingFitResultGap the cluster where we are investigating the existence of gaps and its sliding fit result
     * @param pSecondCluster, SecondSlidingFitResultGap the second cluster in the connection - same with third
     * @param gapInU,V,W, booleans to receive whether there is a gap in each view 
     */    
    bool CheckXPositionInGap(const float &xSample, const pandora::Cluster *const pClusterGap, const TwoDSlidingFitResult &slidingFitResultGap, const pandora::Cluster *const pSecondCluster, const TwoDSlidingFitResult &secondSlidingFitResult, const pandora::Cluster *const pThirdCluster, const TwoDSlidingFitResult &thirdSlidingFitResult, bool &gapInFirst, bool &gapInSecond, bool &gapInThird)const;

    /**                                                                                                                                                                 
     * @brief Check whether a x position is at the end of the cluster
     * @param xSample, the x coordinate of the point tested 
     * @param pCluster the cluster we are interrogating for its extreme coordinates
     */
    bool IsEndOfCluster(const float &xSample, const pandora::Cluster *const pCluster, const TwoDSlidingFitResult &slidingFitResult) const;

    /**                                                                                                                                                      
     * @brief Check whether the external point at x=xSample is in a gap for a given cluster
     *
     * @param xSample, the x coordinate of the external point tested
     * @param pCluster, the cluster we are investigating 
     * @param slidingFitResult, the TwoDSlidingFitResult for the cluster
     */
    bool CheckGaps(const float &xSample, const pandora::Cluster *const pCluster, const TwoDSlidingFitResult &slidingFitResult) const;


    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (in each view) for particle creation
    unsigned int    m_minMatchedSamplingPointRatio;     ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution
    float           m_maxGapTolerance;                  ///< The max gap tolerance
    float           m_sampleStepSize;                   ///< The sampling step size used in association checks, units cm 
    unsigned int    m_maxAngleRatio;                    ///< The max ratio allowed in the angle
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *TracksCrossingGapsTool::Factory::CreateAlgorithmTool() const
{
    return new TracksCrossingGapsTool();
}

} // namespace lar_content

#endif // #ifndef LONG_TRACKS_TOOL_H
