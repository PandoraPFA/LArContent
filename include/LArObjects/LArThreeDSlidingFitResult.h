/**
 *  @file   LArContent/include/LArObjects/LArThreeDSlidingFitResult.h
 *
 *  @brief  Header file for the lar three dimensional sliding fit result class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_SLIDING_FIT_RESULT_H
#define LAR_THREE_D_SLIDING_FIT_RESULT_H 1

#include "Api/PandoraApi.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

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
     *  @param  pCluster address of the candidate cluster
     *  @param  slidingFitWindow the sliding fit window
     *  @param  slidingFitLayerPitch the sliding fit z pitch, units cm
     */
    ThreeDSlidingFitResult(const pandora::Cluster *const pCluster, const unsigned int slidingFitWindow, const float slidingFitLayerPitch);

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
     *  @brief  Calculate the position and direction of the primary axis
     *
     *  @param  pCluster the address of the input cluster
     *
     *  @return pandora::TrackState object containing the position and direction
     */
    static pandora::TrackState GetPrimaryAxis(const pandora::Cluster *pCluster);

    /**
     *  @brief  Generate a seed vector to be used in calculating the orthogonal axes
     *
     *  @param  axisDirection  the primary axis
     *
     *  @return the seed direction vector
     */
    static pandora::CartesianVector GetSeedDirection(const pandora::CartesianVector &axisDirection);

private:

    const pandora::Cluster          *m_pCluster;               ///< The address of the cluster
    const pandora::TrackState        m_primaryAxis;            ///< The primary axis position and direction
    const pandora::CartesianVector   m_axisIntercept;          ///< The axis intercept position
    const pandora::CartesianVector   m_axisDirection;          ///< The axis direction vector
    const pandora::CartesianVector   m_firstOrthoDirection;    ///< The orthogonal direction vector
    const pandora::CartesianVector   m_secondOrthoDirection;   ///< The orthogonal direction vector
    const TwoDSlidingFitResult       m_firstFitResult;         ///< The first sliding fit result
    const TwoDSlidingFitResult       m_secondFitResult;        ///< The second sliding fit result
};

typedef std::vector<ThreeDSlidingFitResult> ThreeDSlidingFitResultList;
typedef std::map<pandora::Cluster*, ThreeDSlidingFitResult> ThreeDSlidingFitResultMap;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *ThreeDSlidingFitResult::GetCluster() const
{
    return m_pCluster;
}

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

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_SLIDING_FIT_RESULT_H
