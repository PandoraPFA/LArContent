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

/**
 *  @brief  DeltaRayShowerHitsTool class
 */
class DeltaRayShowerHitsTool : public HitCreationBaseTool
{
public:
    virtual void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector);

private:
    /**
     *  @brief  Create three dimensional hits, using a list of input two dimensional hits and the 3D hits from the parent particle
     *
     *  @param  inputTwoDHits the vector of input two dimensional hits
     *  @param  parentHits3D the vector of 3D hits from the parent particle
     *  @param  protoHitVector to receive the new three dimensional proto hits
     */
    void CreateDeltaRayShowerHits3D(
        const pandora::CaloHitVector &inputTwoDHits, const pandora::CaloHitVector &parentHits3D, ProtoHitVector &protoHitVector) const;
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_SHOWER_HITS_TOOL_H
