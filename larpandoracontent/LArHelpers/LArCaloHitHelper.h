/**
 *  @file   larpandoracontent/LArHelpers/LArCaloHitHelper.h
 *
 *  @brief  Header file for the calo hit helper class.
 *
 *  $Log: $
 */
#ifndef LAR_CALO_HIT_HELPER_H
#define LAR_CALO_HIT_HELPER_H 1

#include "Objects/CaloHit.h"

namespace lar_content
{

/**
 *  @brief  LArCaloHitHelper class
 */
class LArCaloHitHelper
{
public:

    /**
     *  @brief  Get closest distance between a specified position vector and a list of position vectors
     *
     *  @param  position the position vector
     *  @param  positionVector the list of position vectors
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::CartesianVector &position, const pandora::CartesianPointVector &positionVector);

    /**
     *  @brief  Get closest position from a list of positions to a specified input position vector
     *
     *  @param  position the position vector
     *  @param  positionVector the list of position vectors
     *
     *  @return the closest position
     */
    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &position, const pandora::CartesianPointVector &positionVector);

    /**
     *  @brief  Determine whether two hits are in the same TPC volume
     *
     *  @param  pPrevHit the first CaloHit
     *  @param  pCurrentHit the second CaloHit
     *
     *  @return whether two hits are in the same TPC volume
     */
    static bool IsInSameVolume(const pandora::CaloHit *const pPrevHit, const pandora::CaloHit *const pCurrentHit);
};

} // namespace lar_content

#endif // #ifndef LAR_CALO_HIT_HELPER_H
