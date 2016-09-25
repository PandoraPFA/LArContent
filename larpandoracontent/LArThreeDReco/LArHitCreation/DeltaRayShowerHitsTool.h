/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.h
 * 
 *  @brief  Header file for the delta ray shower hits tool
 * 
 *  $Log: $
 */
#ifndef DELTA_RAY_SHOWER_HITS_TOOL_H
#define DELTA_RAY_SHOWER_HITS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

namespace lar_content
{

class ThreeDHitCreationAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DeltaRayShowerHitsTool class
 */
class DeltaRayShowerHitsTool : public HitCreationBaseTool
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

    void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitVector &inputTwoDHits,
        pandora::CaloHitVector &newThreeDHits);

private:
     /**
     *  @brief  Create three dimensional hits, using a list of input two dimensional hits and the 3D hits from the parent particle
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  inputTwoDHits the vector of input two dimensional hits
     *  @param  caloHitList3D the vector of 3D hits from the parent particle
     *  @param  newThreeDHits to receive the new three dimensional hits
     */
    void CreateThreeDHits(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::CaloHitVector &inputTwoDHits, const pandora::CaloHitVector &caloHitList3D,
        pandora::CaloHitVector &newThreeDHits) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *DeltaRayShowerHitsTool::Factory::CreateAlgorithmTool() const
{
    return new DeltaRayShowerHitsTool();
}

} // namespace lar_content

#endif // #ifndef DELTA_RAY_SHOWER_HITS_TOOL_H
