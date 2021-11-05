/**
 *  @file   larpandoracontent/LArCheating/ThreeDHitRemovalAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_HIT_REMOVAL_ALGORITHM_H
#define LAR_THREE_D_HIT_REMOVAL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  ThreeDHitRemovalAlgorithm class
 */
class ThreeDHitRemovalAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDHitRemovalAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void RemoveSpacePoints(const pandora::ParticleFlowObject *const pCheatedPfo) const;

    pandora::StringVector m_pfoListNames; 
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_HIT_REMOVAL_ALGORITHM_H
