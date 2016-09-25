/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h
 *
 *  @brief  Header file for the shower hits base tool.
 *
 *  $Log: $
 */
#ifndef SHOWER_HITS_BASE_TOOL_H
#define SHOWER_HITS_BASE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

namespace lar_content
{

class ThreeDHitCreationAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ShowerHitsBaseTool class
 */
class ShowerHitsBaseTool : public HitCreationBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerHitsBaseTool();

    virtual void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitVector &inputTwoDHits,
        pandora::CaloHitVector &newThreeDHits);

protected:
    /**
     *  @brief  Create three dimensional hits, using a list of input two dimensional hits and the hits (contained in the same particle)
     *          from the other two views
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  caloHitVector1 hits in the first alternate view
     *  @param  caloHitVector2 hits in the second alternate view
     *  @param  newThreeDHits to receive the new three dimensional hits
     */
    virtual void CreateThreeDHits(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::CaloHitVector &inputTwoDHits, const pandora::CaloHitVector &caloHitVector1,
        const pandora::CaloHitVector &caloHitVector2, pandora::CaloHitVector &newThreeDHits) const;

    /**
     *  @brief  Filter a list of calo hits to find those within a specified tolerance of a give x position
     *
     *  @param  x the x position
     *  @param  xTolerance the x tolerance
     *  @param  inputCaloHitVector the input calo hit vector
     *  @param  outputCaloHitVector to receive the output calo hit vector
     */
    virtual void FilterCaloHits(const float x, const float xTolerance, const pandora::CaloHitVector &inputCaloHitVector, pandora::CaloHitVector &outputCaloHitVector) const;

    /**
     *  @brief  Get the three dimensional position for to a two dimensional calo hit, using the hit and a list of candidate matched
     *          hits in the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  caloHitVector1 the vector of candidate hits in view 1
     *  @param  caloHitVector2 the vector of candidate hits in view 2
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetThreeDPosition(const pandora::CaloHit *const pCaloHit2D, const pandora::CaloHitVector &caloHitVector1, const pandora::CaloHitVector &caloHitVector2,
        pandora::CartesianVector &position3D, float &chiSquared) const = 0;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float     m_xTolerance;          ///< The x tolerance to use when looking for associated calo hits between views
    float     m_chiSquaredCut;       ///< The chi squared cut (accept only values below the cut value)
};

} // namespace lar_content

#endif // #ifndef SHOWER_HITS_BASE_TOOL_H
