/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h
 *
 *  @brief  Header file for the shower hits base tool.
 *
 *  $Log: $
 */
#ifndef SHOWER_HITS_BASE_TOOL_H
#define SHOWER_HITS_BASE_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

namespace lar
{

class ThreeDHitCreationAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ShowerHitsBaseTool class
 */
class ShowerHitsBaseTool : public HitCreationBaseTool
{
public:
    virtual void Run(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &inputTwoDHits,
        pandora::CaloHitList &newThreeDHits);

protected:
    /**
     *  @brief  Create three dimensional hits, using a list of input two dimensional hits and the hits (contained in the same particle)
     *          from the other two views
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  caloHitList1 the first
     *  @param  caloHitList2 the second
     *  @param  newThreeDHits to receive the new three dimensional hits
     */
    virtual void CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::CaloHitList &inputTwoDHits, const pandora::CaloHitList &caloHitList1,
        const pandora::CaloHitList &caloHitList2, pandora::CaloHitList &newThreeDHits) const;

    /**
     *  @brief  Filter a list of calo hits to find those within a specified tolerance of a give x position
     *
     *  @param  x the x position
     *  @param  xTolerance the x tolerance
     *  @param  inputCaloHitList the input calo hit list
     *  @param  outputCaloHitList to receive the output calo hit list
     */
    virtual void FilterCaloHits(const float x, const float xTolerance, const pandora::CaloHitList &inputCaloHitList, pandora::CaloHitList &outputCaloHitList) const;

    /**
     *  @brief  Get the three dimensional position for to a two dimensional calo hit, using the hit and a list of candidate matched
     *          hits in the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  caloHitList1 the list of candidate hits in view 1
     *  @param  caloHitList2 the list of candidate hits in view 2
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetThreeDPosition(const pandora::CaloHit *const pCaloHit2D, const pandora::CaloHitList &caloHitList1, const pandora::CaloHitList &caloHitList2,
        pandora::CartesianVector &position3D, float &chiSquared) const = 0;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float     m_xTolerance;          ///< The x tolerance to use when looking for associated calo hits between views
    float     m_chiSquaredCut;       ///< The chi squared cut (accept only values below the cut value)
};

} // namespace lar

#endif // #ifndef SHOWER_HITS_BASE_TOOL_H
