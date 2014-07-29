/**
 *  @file   LArContent/include/LArObjects/LArTwoDSlidingFitResult.h
 *
 *  @brief  Header file for the lar two dimensional sliding fit result class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_SLIDING_FIT_RESULT_H
#define LAR_TWO_D_SLIDING_FIT_RESULT_H 1

#include "Api/PandoraApi.h"

namespace lar
{

/**
 *  @brief  TwoDSlidingFitResult class
 */
class TwoDSlidingFitResult
{
public:
    /**
     *  @brief  Constructor
     */
    TwoDSlidingFitResult();

    /**
     *  @brief  class LayerFitResult
     */
    class LayerFitResult
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  l the l coordinate
         *  @param  fitT the fitted t coordinate
         *  @param  gradient the fitted gradient dt/dl
         *  @param  rms the rms of the fit residuals
         */
        LayerFitResult(const double l, const double fitT, const double gradient, const double rms);

        /**
         *  @brief  Get the l coordinate
         *
         *  @return the l coordinate
         */
        double GetL() const;

        /**
         *  @brief  Get the fitted t coordinate
         *
         *  @return the fitted t coordinate
         */
        double GetFitT() const;

        /**
         *  @brief  Get the fitted gradient dt/dz
         *
         *  @return the fitted gradient dt/dl
         */
        double GetGradient() const;

        /**
         *  @brief  Get the rms of the fit residuals
         *
         *  @return the rms of the fit residuals
         */
        double GetRms() const;

    private:
        double                  m_l;                            ///< The l coordinate
        double                  m_fitT;                         ///< The fitted t coordinate
        double                  m_gradient;                     ///< The fitted gradient dt/dl
        double                  m_rms;                          ///< The rms of the fit residuals
    };

    typedef std::map<int, LayerFitResult> LayerFitResultMap;

    /**
     *  @brief  LayerFitContribution class
     */
    class LayerFitContribution
    {
    public:
        /**
         *  @brief  Default constructor
         */
        LayerFitContribution();

        /**
         *  @brief  Add point to layer fit
         *
         *  @param  l the longitudinal coordinate
         *  @param  t the transverse coordinate
         */
        void AddPoint(const float l, const float t);

        /**
         *  @brief  Get the sum t
         *
         *  @return the sum t
         */
        double GetSumT() const;

        /**
         *  @brief  Get the sum l
         *
         *  @return the sum l
         */
        double GetSumL() const;

        /**
         *  @brief  Get the sum t * t
         *
         *  @return the sum t * t
         */
        double GetSumTT() const;

        /**
         *  @brief  Get the sum l * t
         *
         *  @return the sum l * t
         */
        double GetSumLT() const;

        /**
         *  @brief  Get the sum l * l
         *
         *  @return the sum z * z
         */
        double GetSumLL() const;

        /**
         *  @brief  Get the number of points used
         *
         *  @return the number of points used
         */
        unsigned int GetNPoints() const;

    private:
        double                  m_sumT;                ///< The sum t
        double                  m_sumL;                ///< The sum l
        double                  m_sumTT;               ///< The sum t * t
        double                  m_sumLT;               ///< The sum l * t
        double                  m_sumLL;               ///< The sum l * l
        unsigned int            m_nPoints;             ///< The number of points used
    };

    typedef std::map<int, LayerFitContribution> LayerFitContributionMap;

    /**
     *  @brief  LayerInterpolation class
     */
    class LayerInterpolation
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param firstayerIter  the iterator for the upstream layer
         *  @param secondLayerIter  the iterator for the downstream layer
         *  @param firstWeight  the weight to be applied to the upstream layer
         *  @param secondWeight  the weight to be applied to the downstream layer
         */
        LayerInterpolation(const LayerFitResultMap::const_iterator &firstLayerIter, const LayerFitResultMap::const_iterator &secondLayerIter,
            const float firstWeight, const float secondWeight);

        /**
         *  @brief  Get the start layer iterator
         *
         *  @return the iterator for the start layer
         */
        LayerFitResultMap::const_iterator GetStartLayerIter() const;

        /**
         *  @brief  Get the end layer iterator
         *
         *  @return the iterator for the end layer
         */
        LayerFitResultMap::const_iterator GetEndLayerIter() const;

        /**
         *  @brief  Get the start layer weight
        *
         *  @return the weight for the start layer
         */
        float GetStartLayerWeight() const;

        /**
         *  @brief  Get the end layer weight
         *
         *  @return the weight for the end layer
         */
        float GetEndLayerWeight() const;

    private:
        LayerFitResultMap::const_iterator m_startLayerIter;     ///< The start layer iterator
        LayerFitResultMap::const_iterator m_endLayerIter;       ///< The end layer iterator
        float                             m_startLayerWeight;   ///< The start layer weight
        float                             m_endLayerWeight;     ///< The end layer weight
    };

    typedef std::vector<LayerInterpolation> LayerInterpolationList;

    /**
     *  @brief  FitSegment class
     */
    class FitSegment
    {
        public:
        /**
         *  @brief Constructor
         *
         *  @param startLayer the start layer
         *  @param endLayer the end layer
         *  @param startX the x position at the start layer
         *  @param endX the x position at the end layer
         */
        FitSegment(const int startLayer, const int endLayer, const float startX, const float endX);

        /**
         *  @brief  Get start layer
         *
         *  @return the start layer
         */
        int GetStartLayer() const;

        /**
         *  @brief  Get end layer
         *
         *  @return the end layer
         */
        int GetEndLayer() const;

        /**
         *  @brief  Get the minimum x value
         *
         *  @return the minimum x value
         */
        float GetMinX() const;

        /**
         *  @brief  Get the maximum x value
         *
         *  @return the maximum x value
         */
        float GetMaxX() const;

        /**
         *  @brief  Whether the x coordinate increases between the start and end layers
         *
         *  @return boolean
         */
        bool IsIncreasingX() const;

    private:
        int         m_startLayer;               ///< The start layer
        int         m_endLayer;                 ///< The end layer
        float       m_minX;                     ///< The minimum x value
        float       m_maxX;                     ///< The maximum x value
        bool        m_isIncreasingX;            ///< Whether the x coordinate increases between the start and end layers
    };

    typedef std::vector<FitSegment> FitSegmentList;

    /**
     *  @brief  Get the address of the cluster
     *
     *  @return the address of the cluster
     */
    const pandora::Cluster *GetCluster() const;

    /**
     *  @brief  Get the layer fit half window
     *
     *  @return the layer fit half window
     */
    unsigned int GetLayerFitHalfWindow() const;

     /**
     *  @brief  Get the layer fit half window length
     *
     *  @return the layer fit half window length
     */
    float GetLayerFitHalfWindowLength() const;

    /**
     *  @brief  Get the axis intercept position
     *
     *  @return the axis intercept position
     */
    const pandora::CartesianVector &GetAxisIntercept() const;

    /**
     *  @brief  Get the axis direction vector
     *
     *  @return the axis direction vector
     */
    const pandora::CartesianVector &GetAxisDirection() const;

    /**
     *  @brief  Get the minimum occupied layer in the sliding fit
     *
     *  @param  the minimum occupied layer in the sliding fit
     */
    int GetMinLayer() const;

    /**
     *  @brief  Get the maximum occupied layer in the sliding fit
     *
     *  @param  the maximum occupied layer in the sliding fit
     */
    int GetMaxLayer() const;

    /**
     *  @brief  Get the minimum and maximum x coordinates associated with the sliding fit
     *
     *  @param  to receive the min x value
     *  @param  to receive the max x value
     */
    void GetMinAndMaxX(float &minX, float &maxX) const;

    /**
     *  @brief  Get the minimum and maximum z coordinates associated with the sliding fit
     *
     *  @param  to receive the min z value
     *  @param  to receive the max z value
     */
    void GetMinAndMaxZ(float &minZ, float &maxZ) const;

    /**
     *  @brief  Get layer number for given sliding linear fit longitudinal coordinate
     *
     *  @param  rL the longitudinal coordinate
     */
    int GetLayer(const float rL) const;

    /**
     *  @brief  Get longitudinal coordinate for a given sliding linear fit layer number
     *
     *  @param  layer the layer number
     */
    float GetL(const int layer) const;

    /**
     *  @brief  Get local sliding fit coordinates for a given global position
     *
     *  @param  position the position cartesian vector
     *  @param  rL to receive the longitudinal coordinate
     *  @param  rT to receive the transverse coordinate
     */
    void GetLocalPosition(const pandora::CartesianVector &position, float &rL, float &rT) const;

    /**
     *  @brief  Get local sliding fit gradient for a given global direction
     *
     *  @param  direction the direction cartesian vector
     *  @param  dTdL to receive the local gradient
     */
    void GetLocalDirection(const pandora::CartesianVector &direction, float &dTdL) const;

    /**
     *  @brief  Get global coordinates for given sliding linear fit coordinates
     *
     *  @param  rL the longitudinal coordinate
     *  @param  rT the transverse coordinate
     *  @param  position to receive the position cartesian vector
     */
    void GetGlobalPosition(const float rL, const float rT, pandora::CartesianVector &position) const;

    /**
     *  @brief  Get global direction coordinates for given sliding linear fit gradient
     *
     *  @param  dTdL the transverse coordinate
     *  @param  direction to receive the direction cartesian vector
     */
    void GetGlobalDirection(const float dTdL, pandora::CartesianVector &direction) const;

    /**
     *  @brief  Get global position corresponding to the fit result in minimum fit layer
     *
     *  @return the position
     */
    pandora::CartesianVector GetGlobalMinLayerPosition() const;

    /**
     *  @brief  Get global position corresponding to the fit result in maximum fit layer
     *
     *  @return the position
     */
    pandora::CartesianVector GetGlobalMaxLayerPosition() const;

    /**
     *  @brief  Get global direction corresponding to the fit result in minimum fit layer
     *
     *  @return the position
     */
    pandora::CartesianVector GetGlobalMinLayerDirection() const;

    /**
     *  @brief  Get global direction corresponding to the fit result in maximum fit layer
     *
     *  @return the position
     */
    pandora::CartesianVector GetGlobalMaxLayerDirection() const;

    /**
     *  @brief  Get rms at minimum layer
     *
     *  @return the rms
     */
    float GetMinLayerRms() const;

    /**
     *  @brief  Get rms at maximum layer
     *
     *  @return the rms
     */
    float GetMaxLayerRms() const;

    /**
     *  @brief  Get fit rms for a given longitudinal coordinate
     *
     *  @param  rL the longitudinal coordinate
     *
     *  @return the fit rms
     */
    float GetFitRms(const float rL) const;

    /**
     *  @brief  Get scattering angle for a given longitudinal coordinate
     *
     *  @param  rL the longitudinal coordinate
     */
    float GetCosScatteringAngle(const float rL) const;

    /**
     *  @brief  Get global fit position for a given longitudinal coordinate
     *
     *  @param  rL the longitudinal coordinate
     *  @param  position the fitted position at these coordinates
     */
    void GetGlobalFitPosition(const float rL, pandora::CartesianVector &position) const;

    /**
     *  @brief  Get global fit direction for a given longitudinal coordinate
     *
     *  @param  rL the longitudinal coordinate
     *  @param  direction the fitted direction at these coordinates
     */
    void GetGlobalFitDirection(const float rL, pandora::CartesianVector &direction) const;

    /**
     *  @brief  Get global fit position for a given input x coordinate
     *
     *  @param  x the input coordinate
     *  @param  position the fitted position at these coordinates
     */
    void GetGlobalFitPositionAtX(const float x, pandora::CartesianVector &position) const;

    /**
     *  @brief  Get global fit direction for a given input x coordinate
     *
     *  @param  x the input coordinate
     *  @param  direction the fitted direction at these coordinates
     */
    void GetGlobalFitDirectionAtX(const float x, pandora::CartesianVector &direction) const;

    /**
     *  @brief  Get projected position on global fit for a given position vector
     *
     *  @param  inputPosition the input coordinate
     *  @param  projectedPosition the projected position on the global fit for these coordinates
     */
    void GetGlobalFitProjection(const pandora::CartesianVector &inputPosition, pandora::CartesianVector &projectedPosition) const;

    /**
     *  @brief Get a list of projected positions for a given input x coordinate
     *
     *  @param x the input x coordinate
     *  @param positionList the output list of positions
     */
    void GetGlobalFitPositionListAtX(const float x, pandora::CartesianPointList &positionList) const;

    /**
     *  @brief Get projected position for a given input x coordinate and fit segment
     *
     *  @param x the input x coordinate
     *  @param fitSegment the portion of sliding linear fit
     *  @param position the output position
     */
    void GetTransverseProjection(const float x, const FitSegment &fitSegment, pandora::CartesianVector &position) const;

    /**
     *  @brief Get projected position and direction for a given input x coordinate and fit segment
     *
     *  @param x the input x coordinate
     *  @param fitSegment the portion of sliding linear fit
     *  @param position the output position
     *  @param position the output direction
     */
    void GetTransverseProjection(const float x, const FitSegment &fitSegment, pandora::CartesianVector &position,
        pandora::CartesianVector &direction) const;

    /**
     *  @brief  Get extrapolated position (beyond span) for a given input x coordinate
     *
     *  @param  x the input coordinate
     *  @param  position the extrapolated position at these coordinates
     */
    void GetExtrapolatedPositionAtX(const float x, pandora::CartesianVector &position) const;

    /**
     *  @brief Get fit segment for a given longitudinal coordinate
     *
     *  @param  rL the longitudinal coordinate
     */
    FitSegment GetFitSegment(const float rL) const;

    /**
     *  @brief  Get the layer fit result map
     *
     *  @return the layer fit result map
     */
    const LayerFitResultMap &GetLayerFitResultMap() const;

    /**
     *  @brief  Get the layer fit contribution map
     *
     *  @return the layer fit contribution map
     */
    const LayerFitContributionMap &GetLayerFitContributionMap() const;

     /**
     *  @brief  Get the fit segment list
     *
     *  @return the fit segment list
     */
    const FitSegmentList &GetFitSegmentList() const;

private:
    /**
     *  @brief  Get the minimum and maximum x or z coordinates associated with the sliding fit
     *
     *  @param  isX whether to provide extremal x or z coordinates
     *  @param  to receive the min coordinate value
     *  @param  to receive the max coordinate value
     */
    void GetMinAndMaxCoordinate(const bool isX, float &min, float &max) const;

    /**
     *  @brief Interpolate a position between two layers
     *
     *  @param layerInterpolation the pair of surrounding layers
     *
     *  @return the interpolated position
     */
    pandora::CartesianVector GetGlobalFitPosition(const LayerInterpolation &layerInterpolation) const;

    /**
     *  @brief Interpolate a direction between two layers
     *
     *  @param layerInterpolation the pair of surrounding layers
     *
     *  @return the interpolated direction
     */
    pandora::CartesianVector GetGlobalFitDirection(const LayerInterpolation &layerInterpolation) const;

    /**
     *  @brief Interpolate a rms between two layers
     *
     *  @param layerInterpolation the pair of surrounding layers
     *
     *  @return the interpolated rms
     */
    float GetFitRms(const LayerInterpolation &layerInterpolation) const;

    /**
     *  @brief  Get the pair of layers surrounding a specified longitudinal position
     *
     *  @param  rL the longitudinal coordinate
     *
     *  @return Layer interpolation object
     */
    LayerInterpolation LongitudinalInterpolation(const float rL) const;

    /**
     *  @brief  Get the surrounding pair of layers for a specified transverse position and fit segment
     *
     *  @param  x the input coordinate
     *  @param  fitSegment the fit segment
     *
     *  @return Layer interpolation object
     */
    LayerInterpolation TransverseInterpolation(const float x, const FitSegment &fitSegment) const;

    /**
     *  @brief  Get the a list of surrounding layer pairs for a specified transverse position
     *
     *  @param  x the input coordinate
     *  @param  layerInterpolationList the output list of layer interpolation objects
     */
    void TransverseInterpolation(const float x, LayerInterpolationList &layerInterpolationList) const;

    /**
     *  @brief  Get iterators for layers surrounding the specified longitudinal position
     *
     *  @param  rL the longitudinal coordinate
     *  @param  firstLayerIter to receive the iterator for the layer just below the input coordinate
     *  @param  secondLayerIter to receive the iterator for the layer just above the input coordinate
     */
    void GetLongitudinalSurroundingLayers(const float rL, LayerFitResultMap::const_iterator &firstLayerIter,
        LayerFitResultMap::const_iterator &secondLayerIter) const;

    /**
     *  @brief  Get iterators for layers surrounding a specified transverse position
     *
     *  @param  x the transverse coordinate
     *  @param  minLayer the minimum allowed layer
     *  @param  maxLayer the maximum allowed layer
     *  @param  firstLayerIter to receive the iterator for the layer just below the input coordinate
     *  @param  secondLayerIter to receive the iterator for the layer just above the input coordinate
     */
    void GetTransverseSurroundingLayers(const float x, const int minLayer, const int maxLayer,
        LayerFitResultMap::const_iterator &firstLayerIter, LayerFitResultMap::const_iterator &secondLayerIter) const;

    /**
     *  @brief  Get interpolation weights for layers surrounding a specified longitudinal position
     *
     *  @param  rL the longitudinal coordinate
     *  @param  firstLayerIter the iterator for the layer below the input coordinate
     *  @param  secondLayerIter the iterator for the layer above the input coordinate
     *  @param  firstWeight the weight assigned to the layer below the input coordinate
     *  @param  secondWeight the weight assigned to the layer above the input coordinate
     */
    void GetLongitudinalInterpolationWeights(const float rL, const LayerFitResultMap::const_iterator &firstLayerIter,
        const LayerFitResultMap::const_iterator &secondLayerIter, float &firstWeight, float &secondWeight) const;

   /**
     *  @brief  Get interpolation weights for layers surrounding a specified transverse position
     *
     *  @param  x the transverse coordinate
     *  @param  firstLayerIter the iterator for the layer below the input coordinate
     *  @param  secondLayerIter the iterator for the layer above the input coordinate
     *  @param  firstWeight the weight assigned to the layer below the input coordinate
     *  @param  firstWeight the weight assigned to the layer above the input coordinate
     */
    void GetTransverseInterpolationWeights(const float x, const LayerFitResultMap::const_iterator &firstLayerIter,
        const LayerFitResultMap::const_iterator &secondLayerIter, float &firstWeight, float &secondWeight) const;

    const pandora::Cluster     *m_pCluster;                 ///< The address of the cluster
    unsigned int                m_layerFitHalfWindow;       ///< The layer fit half window
    pandora::CartesianVector    m_axisIntercept;            ///< The axis intercept position
    pandora::CartesianVector    m_axisDirection;            ///< The axis direction vector
    LayerFitResultMap           m_layerFitResultMap;        ///< The layer fit result map
    LayerFitContributionMap     m_layerFitContributionMap;  ///< The layer fit contribution map
    FitSegmentList              m_fitSegmentList;           ///< The fit segment list

    friend class LArClusterHelper;
};

typedef std::vector<TwoDSlidingFitResult> TwoDSlidingFitResultList;
typedef std::map<pandora::Cluster*, TwoDSlidingFitResult> TwoDSlidingFitResultMap;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *TwoDSlidingFitResult::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TwoDSlidingFitResult::GetLayerFitHalfWindow() const
{
    return m_layerFitHalfWindow;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TwoDSlidingFitResult::GetAxisIntercept() const
{
    return m_axisIntercept;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &TwoDSlidingFitResult::GetAxisDirection() const
{
    return m_axisDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TwoDSlidingFitResult::GetMinLayer() const
{
    if (m_layerFitResultMap.empty())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_layerFitResultMap.begin()->first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TwoDSlidingFitResult::GetMaxLayer() const
{
    if (m_layerFitResultMap.empty())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_layerFitResultMap.rbegin()->first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TwoDSlidingFitResult::GetMinAndMaxX(float &minX, float &maxX) const
{
    return this->GetMinAndMaxCoordinate(true, minX, maxX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TwoDSlidingFitResult::GetMinAndMaxZ(float &minZ, float &maxZ) const
{
    return this->GetMinAndMaxCoordinate(false, minZ, maxZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult::LayerFitResultMap &TwoDSlidingFitResult::GetLayerFitResultMap() const
{
    return m_layerFitResultMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult::LayerFitContributionMap &TwoDSlidingFitResult::GetLayerFitContributionMap() const
{
    return m_layerFitContributionMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult::FitSegmentList &TwoDSlidingFitResult::GetFitSegmentList() const
{
    return m_fitSegmentList;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitResult::GetL() const
{
    return m_l;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitResult::GetFitT() const
{
    return m_fitT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitResult::GetGradient() const
{
    return m_gradient;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitResult::GetRms() const
{
    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitContribution::GetSumT() const
{
    return m_sumT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitContribution::GetSumL() const
{
    return m_sumL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitContribution::GetSumLT() const
{
    return m_sumLT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitContribution::GetSumLL() const
{
    return m_sumLL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TwoDSlidingFitResult::LayerFitContribution::GetSumTT() const
{
    return m_sumTT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TwoDSlidingFitResult::LayerFitContribution::GetNPoints() const
{
    return m_nPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TwoDSlidingFitResult::LayerFitResultMap::const_iterator TwoDSlidingFitResult::LayerInterpolation::GetStartLayerIter() const
{
    return m_startLayerIter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TwoDSlidingFitResult::LayerFitResultMap::const_iterator TwoDSlidingFitResult::LayerInterpolation::GetEndLayerIter() const
{
    return m_endLayerIter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoDSlidingFitResult::LayerInterpolation::GetStartLayerWeight() const
{
    return m_startLayerWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoDSlidingFitResult::LayerInterpolation::GetEndLayerWeight() const
{
    return m_endLayerWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline int TwoDSlidingFitResult::FitSegment::GetStartLayer() const
{
    return m_startLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TwoDSlidingFitResult::FitSegment::GetEndLayer() const
{
    return m_endLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoDSlidingFitResult::FitSegment::GetMinX() const
{
    return m_minX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoDSlidingFitResult::FitSegment::GetMaxX() const
{
    return m_maxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TwoDSlidingFitResult::FitSegment::IsIncreasingX() const
{
    return m_isIncreasingX;
}

} // namespace lar

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_RESULT_H
