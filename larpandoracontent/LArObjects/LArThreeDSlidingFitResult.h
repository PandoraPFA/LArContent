/**
 *  @file   larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h
 *
 *  @brief  Header file for the lar three dimensional sliding fit result class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_SLIDING_FIT_RESULT_H
#define LAR_THREE_D_SLIDING_FIT_RESULT_H 1

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ThreeDSlidingFitResult class
 */
class ThreeDSlidingFitResult
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pT describing the positions to be fitted
     *  @param  slidingFitWindow the sliding fit window
     *  @param  slidingFitLayerPitch the sliding fit z pitch, units cm
     */
    template <typename T>
    ThreeDSlidingFitResult(const T *const pT, const unsigned int slidingFitWindow, const float slidingFitLayerPitch);

    /**
     *  @brief  Get the address of the cluster
     *
     *  @return the address of the cluster
     */
    const pandora::Cluster *GetCluster() const;

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
     *  @brief  Get the first sliding fit result for this cluster
     *
     *  @return the first sliding fit result for this cluster
     */
    const TwoDSlidingFitResult &GetFirstFitResult() const;

    /**
     *  @brief  Get the second sliding fit result for this cluster
     *
     *  @return the second sliding fit result for this cluster
     */
    const TwoDSlidingFitResult &GetSecondFitResult() const;

    /**
     *  @brief  Get global position corresponding to the fit result in minimum fit layer
     *
     *  @return the position
     */
    const pandora::CartesianVector &GetGlobalMinLayerPosition() const;

    /**
     *  @brief  Get global position corresponding to the fit result in maximum fit layer
     *
     *  @return the position
     */
    const pandora::CartesianVector &GetGlobalMaxLayerPosition() const;

    /**
     *  @brief  Get global direction corresponding to the fit result in minimum fit layer
     *
     *  @return the position
     */
    const pandora::CartesianVector &GetGlobalMinLayerDirection() const;

    /**
     *  @brief  Get global direction corresponding to the fit result in maximum fit layer
     *
     *  @return the position
     */
    const pandora::CartesianVector &GetGlobalMaxLayerDirection() const;

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
     *  @brief  Get longitudinal projection onto primary axis
     *
     *  @param  position the input coordinates
     *
     *  @return the longitudinal distance along the axis direction from the axis intercept
     */
    float GetLongitudinalDisplacement(const pandora::CartesianVector &position) const;

    /**
     *  @brief  Get global fit position for a given longitudinal coordinate
     *
     *  @param  rL the longitudinal coordinate
     *  @param  position the fitted position at these coordinates
     *
     *  @return status code, faster than throwing in regular use-cases
     */
    pandora::StatusCode GetGlobalFitPosition(const float rL, pandora::CartesianVector &position) const;

    /**
     *  @brief  Get global fit direction for a given longitudinal coordinate
     *
     *  @param  rL the longitudinal coordinate
     *  @param  direction the fitted direction at these coordinates
     *
     *  @return status code, faster than throwing in regular use-cases
     */
    pandora::StatusCode GetGlobalFitDirection(const float rL, pandora::CartesianVector &direction) const;

private:
    /**
     *  @brief  Get global coordinates for a given pair of sliding linear fit coordinates
     *
     *  @param  rL the longitudinal coordinate
     *  @param  rT1 the first transverse coordinate
     *  @param  rT2 the second transverse coordinate
     *  @param  position to receive the position cartesian vector
     */
    void GetGlobalPosition(const float rL, const float rT1, const float rT2, pandora::CartesianVector &position) const;

    /**
     *  @brief  Get global direction coordinates for a given pair of sliding linear fit gradients
     *
     *  @param  dTdL1 the first transverse coordinate
     *  @param  dTdL2 the second transverse coordinate
     *  @param  direction to receive the direction cartesian vector
     */
    void GetGlobalDirection(const float dTdL1, const float dTdL2, pandora::CartesianVector &direction) const;

    /**
     *  @brief  Calculate the position and direction of the primary axis
     *
     *  @param  pCluster the address of the input cluster
     *  @param  slidingFitLayerPitch the sliding fit z pitch, units cm
     *
     *  @return pandora::TrackState object containing the position and direction
     */
    static pandora::TrackState GetPrimaryAxis(const pandora::Cluster *const pCluster, const float slidingFitLayerPitch);

    /**
     *  @brief  Calculate the position and direction of the primary axis
     *
     *  @param  pPointVector the address of the input point vector
     *  @param  slidingFitLayerPitch the sliding fit z pitch, units cm
     *
     *  @return pandora::TrackState object containing the position and direction
     */
    static pandora::TrackState GetPrimaryAxis(const pandora::CartesianPointVector *const pPointVector, const float slidingFitLayerPitch);

    /**
     *  @brief  Generate a seed vector to be used in calculating the orthogonal axes
     *
     *  @param  axisDirection  the primary axis
     *
     *  @return the seed direction vector
     */
    static pandora::CartesianVector GetSeedDirection(const pandora::CartesianVector &axisDirection);

    const pandora::TrackState         m_primaryAxis;            ///< The primary axis position and direction
    const pandora::CartesianVector    m_axisIntercept;          ///< The axis intercept position
    const pandora::CartesianVector    m_axisDirection;          ///< The axis direction vector
    const pandora::CartesianVector    m_firstOrthoDirection;    ///< The orthogonal direction vector
    const pandora::CartesianVector    m_secondOrthoDirection;   ///< The orthogonal direction vector
    const TwoDSlidingFitResult        m_firstFitResult;         ///< The first sliding fit result
    const TwoDSlidingFitResult        m_secondFitResult;        ///< The second sliding fit result
    const int                         m_minLayer;               ///< The minimum combined layer
    const int                         m_maxLayer;               ///< The maximum combined layer

    pandora::CartesianVector          m_minLayerPosition;       ///< The global position at the minimum combined layer
    pandora::CartesianVector          m_maxLayerPosition;       ///< The global position at the maximum combined layer
    pandora::CartesianVector          m_minLayerDirection;      ///< The global direction at the minimum combined layer
    pandora::CartesianVector          m_maxLayerDirection;      ///< The global direction at the maximum combined layer
};

typedef std::vector<ThreeDSlidingFitResult> ThreeDSlidingFitResultList;
typedef std::unordered_map<const pandora::Cluster*, ThreeDSlidingFitResult> ThreeDSlidingFitResultMap;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ThreeDSlidingFitResult::GetAxisIntercept() const
{
    return m_axisIntercept;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ThreeDSlidingFitResult::GetAxisDirection() const
{
    return m_axisDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult &ThreeDSlidingFitResult::GetFirstFitResult() const
{
    return m_firstFitResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult &ThreeDSlidingFitResult::GetSecondFitResult() const
{
    return m_secondFitResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ThreeDSlidingFitResult::GetGlobalMinLayerPosition() const
{
    return m_minLayerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ThreeDSlidingFitResult::GetGlobalMaxLayerPosition() const
{
    return m_maxLayerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ThreeDSlidingFitResult::GetGlobalMinLayerDirection() const
{
    return m_minLayerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ThreeDSlidingFitResult::GetGlobalMaxLayerDirection() const
{
    return m_maxLayerDirection;
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_SLIDING_FIT_RESULT_H
