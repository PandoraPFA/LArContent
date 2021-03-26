/**
 *  @file   larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h
 *
 *  @brief  Header file for the lar three dimensional sliding cone fit result class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_SLIDING_CONE_FIT_RESULT_H
#define LAR_THREE_D_SLIDING_CONE_FIT_RESULT_H 1

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <list>
#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ConeSelection enum
 */
enum ConeSelection
{
    CONE_FORWARD_ONLY,
    CONE_BACKWARD_ONLY,
    CONE_BOTH_DIRECTIONS
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SimpleCone class
 */
class SimpleCone
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  coneApex
     *  @param  coneDirection
     *  @param  coneLength
     *  @param  coneTanHalfAngle
     */
    SimpleCone(const pandora::CartesianVector &coneApex, const pandora::CartesianVector &coneDirection, const float coneLength, const float coneTanHalfAngle);

    /**
     *  @brief  Get the cone apex
     *
     *  @return the cone apex
     */
    const pandora::CartesianVector &GetConeApex() const;

    /**
     *  @brief  Get the cone direction
     *
     *  @return the cone direction
     */
    const pandora::CartesianVector &GetConeDirection() const;

    /**
     *  @brief  Get the cone length
     *
     *  @return the cone length
     */
    float GetConeLength() const;

    /**
     *  @brief  Get the tangent of the cone half-angle
     *
     *  @return the tangent of the cone half-angle
     */
    float GetConeTanHalfAngle() const;

    /**
     *  @brief  Get the mean transverse distance to all hits in a cluster (whether contained or not)
     *
     *  @param  pCluster the address of the cluster
     *
     *  @return the mean transverse distance to all hits (whether contained or not)
     */
    float GetMeanRT(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get the fraction of hits in a provided cluster that are bounded within the cone, using fitted cone angle and length
     *
     *  @param  pCluster the address of the cluster
     *
     *  @return the bounded hit fraction
     */
    float GetBoundedHitFraction(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get the fraction of hits in a provided cluster that are bounded within the cone, using provided cone angle and length
     *
     *  @param  pCluster the address of the cluster
     *  @param  coneLength the provided cone length
     *  @param  coneTanHalfAngle the provided tangent of the cone half-angle
     *
     *  @return the bounded hit fraction
     */
    float GetBoundedHitFraction(const pandora::Cluster *const pCluster, const float coneLength, const float coneTanHalfAngle) const;

private:
    pandora::CartesianVector m_coneApex;      ///< The cone apex
    pandora::CartesianVector m_coneDirection; ///< The cone direction
    float m_coneLength;                       ///< The cone length
    float m_coneTanHalfAngle;                 ///< The tangent of the cone half-angle
};

typedef std::vector<SimpleCone> SimpleConeList;

//------------------------------------------------------------------------------------------------------------------------------------------

typedef std::map<int, pandora::TrackState> TrackStateMap;

/**
 *  @brief  ThreeDSlidingConeFitResult class
 */
class ThreeDSlidingConeFitResult
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pT describing the positions to be fitted
     *  @param  slidingFitWindow the sliding fit window
     *  @param  slidingFitLayerPitch the sliding fit pitch, units cm
     */
    template <typename T>
    ThreeDSlidingConeFitResult(const T *const pT, const unsigned int slidingFitWindow, const float slidingFitLayerPitch);

    /**
     *  @brief  Get the sliding fit result for the full cluster
     *
     *  @return the sliding fit result for the full cluster
     */
    const ThreeDSlidingFitResult &GetSlidingFitResult() const;

    /**
     *  @brief  Get the track state map, which caches results from the sliding fit result
     *
     *  @return the track state map
     */
    const TrackStateMap &GetTrackStateMap() const;

    /**
     *  @brief  Get the list of simple cones fitted to the three dimensional cluster
     *
     *  @param  nLayersForConeFit the number of layer to use to extract the cone direction
     *  @param  nCones the number of cones to extract from the cluster (spaced uniformly along the cluster)
     *  @param  coneSelection whether to receive forwards or backwards (or both) cones
     *  @param  simpleConeList to receive the simple cone list
     */
    void GetSimpleConeList(const unsigned int nLayersForConeFit, const unsigned int nCones, const ConeSelection coneSelection,
        SimpleConeList &simpleConeList) const;

private:
    typedef std::list<pandora::TrackState> TrackStateLinkedList; ///< The track state linked list typedef

    const ThreeDSlidingFitResult m_slidingFitResult; ///< The sliding fit result for the full cluster
    TrackStateMap m_trackStateMap;                   ///< The track state map
};

typedef std::vector<ThreeDSlidingConeFitResult> ThreeDSlidingConeFitResultList;
typedef std::unordered_map<const pandora::Cluster *, ThreeDSlidingConeFitResult> ThreeDSlidingConeFitResultMap;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline SimpleCone::SimpleCone(const pandora::CartesianVector &coneApex, const pandora::CartesianVector &coneDirection,
    const float coneLength, const float coneTanHalfAngle) :
    m_coneApex(coneApex),
    m_coneDirection(coneDirection),
    m_coneLength(coneLength),
    m_coneTanHalfAngle(coneTanHalfAngle)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &SimpleCone::GetConeApex() const
{
    return m_coneApex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &SimpleCone::GetConeDirection() const
{
    return m_coneDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SimpleCone::GetConeLength() const
{
    return m_coneLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SimpleCone::GetConeTanHalfAngle() const
{
    return m_coneTanHalfAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SimpleCone::GetBoundedHitFraction(const pandora::Cluster *const pCluster) const
{
    return this->GetBoundedHitFraction(pCluster, this->GetConeLength(), this->GetConeTanHalfAngle());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDSlidingFitResult &ThreeDSlidingConeFitResult::GetSlidingFitResult() const
{
    return m_slidingFitResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TrackStateMap &ThreeDSlidingConeFitResult::GetTrackStateMap() const
{
    return m_trackStateMap;
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_SLIDING_CONE_FIT_RESULT_H
