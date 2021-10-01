/**
 *  @file   larpandoracontent/LArCheating/CheatingPfoSubDominantHitRemovalAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_PFO_SUB_DOMINANT_HIT_REMOVAL_ALGORITHM_H
#define LAR_CHEATING_PFO_SUB_DOMINANT_HIT_REMOVAL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  CheatingPfoSubDominantHitRemovalAlgorithm class
 */
class CheatingPfoSubDominantHitRemovalAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingPfoSubDominantHitRemovalAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void FindMatchedMCParticle(const pandora::ParticleFlowObject *const pCheatedPfo, const pandora::MCParticle *&pMatchedMCParticle) const;

    void RemoveSubDominantHits(const pandora::ParticleFlowObject *const pCheatedPfo, const pandora::MCParticle *const pMatchedMCParticle) const;

    void RemoveSpacePoints(const pandora::ParticleFlowObject *const pCheatedPfo) const;

    std::string m_cheatedPfoListName; 
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_PFO_SUB_DOMINANT_HIT_REMOVAL_ALGORITHM_H
