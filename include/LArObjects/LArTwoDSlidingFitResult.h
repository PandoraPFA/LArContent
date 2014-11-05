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

#include "LArObjects/LArTwoDSlidingFitObjects.h"

namespace lar_content
{

/**
 *  @brief  TwoDSlidingFitResult class
 */
class TwoDSlidingFitResult
{
public:
    /**
     *  @brief  Constructor using cluster extremal x-z positions to define primary axis
     *
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  layerPitch the layer pitch, units cm
     */
    TwoDSlidingFitResult(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, const float layerPitch);

    /**
     *  @brief  Constructor using specified primary axis
     *
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  layerPitch the layer pitch, units cm
     *  @param  axisIntercept the axis intercept position
     *  @param  axisDirection the axis direction vector
     */
    TwoDSlidingFitResult(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, const float layerPitch,
        const pandora::CartesianVector &axisIntercept, const pandora::CartesianVector &axisDirection);

    /**
     *  @brief  Constructor using specified primary axis and layer fit contribution map. User is responsible for ensuring that
     *          z-pitch, axis intercept and axis direction agree with calculations used to fill the layer fit contribution map.
     *
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  layerPitch the layer pitch, units cm
     *  @param  axisIntercept the axis intercept position
     *  @param  axisDirection the axis direction vector
     *  @param  layerFitContributionMap the layer fit contribution map
     */
    TwoDSlidingFitResult(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, const float layerPitch,
        const pandora::CartesianVector &axisIntercept, const pandora::CartesianVector &axisDirection, const LayerFitContributionMap &layerFitContributionMap);

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
     *  @brief  Get the layer pitch, units cm
     *
     *  @return the layer pitch
     */
    float GetLayerPitch() const;

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

     /**
     *  @brief  Get the layer fit half window length
     *
     *  @return the layer fit half window length
     */
    float GetLayerFitHalfWindowLength() const;

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
    const FitSegment &GetFitSegment(const float rL) const;

private:
    /**
     *  @brief  Fill the layer fit contribution map
     */
    void FillLayerFitContributionMap();

    /**
     *  @brief  Perform the sliding linear fit
     */
    void PerformSlidingLinearFit();

    /**
     *  @brief  Find sliding fit segments; sections with tramsverse direction
     */
    void FindSlidingFitSegments();

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
    float                       m_layerPitch;               ///< The layer pitch, units cm
    pandora::CartesianVector    m_axisIntercept;            ///< The axis intercept position
    pandora::CartesianVector    m_axisDirection;            ///< The axis direction vector
    LayerFitResultMap           m_layerFitResultMap;        ///< The layer fit result map
    LayerFitContributionMap     m_layerFitContributionMap;  ///< The layer fit contribution map
    FitSegmentList              m_fitSegmentList;           ///< The fit segment list
};

typedef std::vector<TwoDSlidingFitResult> TwoDSlidingFitResultList;
typedef std::map<const pandora::Cluster*, TwoDSlidingFitResult> TwoDSlidingFitResultMap;

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

inline float TwoDSlidingFitResult::GetLayerPitch() const
{
    return m_layerPitch;
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

inline const LayerFitResultMap &TwoDSlidingFitResult::GetLayerFitResultMap() const
{
    return m_layerFitResultMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LayerFitContributionMap &TwoDSlidingFitResult::GetLayerFitContributionMap() const
{
    return m_layerFitContributionMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const FitSegmentList &TwoDSlidingFitResult::GetFitSegmentList() const
{
    return m_fitSegmentList;
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

} // namespace lar_content

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_RESULT_H
