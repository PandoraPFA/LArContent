/**
 *  @file   larpandoracontent/LArObjects/LArPfoObjects.h
 *
 *  @brief  Header file for lar pfo objects.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_OBJECTS_H
#define LAR_PFO_OBJECTS_H 1

#include "Objects/CartesianVector.h"
#include "Objects/TrackState.h"

#include <vector>

namespace pandora
{
class CaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  LArTrackState class
 */
class LArTrackState : public pandora::TrackState
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  position the position
     *  @param  direction the direction
     *  @param  pCaloHit the address of the associated calo hit
     */
    LArTrackState(const pandora::CartesianVector &position, const pandora::CartesianVector &direction, const pandora::CaloHit *const pCaloHit);

    /**
     *  @brief  Constructor
     *
     *  @param  position the position
     *  @param  direction the direction
     */
    LArTrackState(const pandora::CartesianVector &position, const pandora::CartesianVector &direction);

    /**
     *  @brief  Return direction at this trajectory point
     *
     *  @return the direction at this trajectory point
     */
    const pandora::CartesianVector &GetDirection() const;

    /**
     *  @brief  Return calo hit at this trajectory point
     *
     *  @param  the address of the calo hit at this trajectory point
     */
    const pandora::CaloHit *GetCaloHit() const;

private:
    const pandora::CaloHit *m_pCaloHit;
};

typedef std::vector<LArTrackState> LArTrackStateVector;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArTrackTrajectoryPoint class
 */
class LArTrackTrajectoryPoint : public std::pair<float, LArTrackState>
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  projectedDistance the projected distance
     *  @param  larTrackState the lar track state
     */
    LArTrackTrajectoryPoint(const float projectedDistance, const LArTrackState &larTrackState);

    /**
     *  @brief  Constructor
     *
     *  @param  projectedDistance the projected distance
     *  @param  larTrackState the lar track state
     *  @param  index the index associated with the trajectory point
     */
    LArTrackTrajectoryPoint(const float projectedDistance, const LArTrackState &larTrackState, const int index);

    /**
     *  @brief  Get the index associated with the trajectory point
     *
     *  @return the index associated with the trajectory point
     */
    int GetIndex() const;

private:
    int m_index; ///< The index associated with the trajectory point
};

typedef std::vector<LArTrackTrajectoryPoint> LArTrackTrajectory;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArShowerPCA class
 */
class LArShowerPCA
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  centroid  centroid of shower
     *  @param  primaryAxis  primary axis of shower
     *  @param  secondaryAxis  secondary axis of shower
     *  @param  tertiaryAxis  tertiary axis of shower
     *  @param  axisLengths  ordered vector of shower lengths
     */
    LArShowerPCA(const pandora::CartesianVector &centroid, const pandora::CartesianVector &primaryAxis,
        const pandora::CartesianVector &secondaryAxis, const pandora::CartesianVector &tertiaryAxis, const pandora::CartesianVector &eigenvalues);

    /**
     *  @brief  Return centroid
     *
     *  @return the centroid
     */
    const pandora::CartesianVector &GetCentroid() const;

    /**
     *  @brief  Return primary axis
     *
     *  @return the primary axis
     */
    const pandora::CartesianVector &GetPrimaryAxis() const;

    /**
     *  @brief  Return secondary axis
     *
     *  @return the secondary axis
     */
    const pandora::CartesianVector &GetSecondaryAxis() const;

    /**
     *  @brief  Return tertiary axis
     *
     *  @return the tertiary axis
     */
    const pandora::CartesianVector &GetTertiaryAxis() const;

    /**
     *  @brief  Return vector of eigenvalues
     *
     *  @return the vector of eigenvalues
     */
    const pandora::CartesianVector &GetEigenValues() const;

    /**
     *  @brief  Return vector of lengths
     *
     *  @return the vector of lengths
     */
    const pandora::CartesianVector &GetAxisLengths() const;

    /**
     *  @brief  Return primary length
     *
     *  @return the primary length
     */
    float GetPrimaryLength() const;

    /**
     *  @brief  Return secondary length
     *
     *  @return the secondary length
     */
    float GetSecondaryLength() const;

    /**
     *  @brief  Return tertiary length
     *
     *  @return the tertiary length
     */
    float GetTertiaryLength() const;

private:
    const pandora::CartesianVector m_centroid;      ///< The centroid
    const pandora::CartesianVector m_primaryAxis;   ///< The primary axis
    const pandora::CartesianVector m_secondaryAxis; ///< The secondary axis
    const pandora::CartesianVector m_tertiaryAxis;  ///< The tertiary axis
    const pandora::CartesianVector m_eigenValues;   ///< The vector of eigenvalues
    const pandora::CartesianVector m_axisLengths;   ///< The vector of lengths
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArTrackTrajectoryPoint::LArTrackTrajectoryPoint(const float projectedDistance, const LArTrackState &larTrackState) :
    std::pair<float, LArTrackState>(projectedDistance, larTrackState), m_index(-1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArTrackTrajectoryPoint::LArTrackTrajectoryPoint(const float projectedDistance, const LArTrackState &larTrackState, const int index) :
    std::pair<float, LArTrackState>(projectedDistance, larTrackState), m_index(index)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArTrackTrajectoryPoint::GetIndex() const
{
    return m_index;
}

} // namespace lar_content

#endif // #ifndef LAR_PFO_OBJECTS_H
