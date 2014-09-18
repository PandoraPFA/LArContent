/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/HitCreationBaseTool.cc
 *
 *  @brief  Header file for the hit creation base tool.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_CREATION_BASE_TOOL_H
#define LAR_HIT_CREATION_BASE_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

class ThreeDHitCreationAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  HitCreationBaseTool class
 */
class HitCreationBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo the address of the pfo
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  newThreeDHits to receive the new three dimensional hits
     */
    virtual void Run(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &inputTwoDHits,
        pandora::CaloHitList &newThreeDHits) = 0;

    /**
     *  @brief  Default constructor
     */
    HitCreationBaseTool();

    /**
     *  @brief  Destructor
     */
    virtual ~HitCreationBaseTool();

protected:
    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  hitType1 the hit type identifying the first view
     *  @param  hitType2 the hit type identifying the second view
     *  @param  fitPositionList1 the candidate sliding fit position in the first view
     *  @param  fitPositionList2 the candidate sliding fit position in the second view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetBestPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::HitType hitType1, const pandora::HitType hitType2,
        const pandora::CartesianPointList &fitPositionList1, const pandora::CartesianPointList &fitPositionList2, pandora::CartesianVector &position3D, float &chiSquared) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  hitType1 the hit type identifying the first view
     *  @param  hitType2 the hit type identifying the second view
     *  @param  fitPosition1 the candidate sliding fit position in the first view
     *  @param  fitPosition2 the candidate sliding fit position in the second view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::HitType hitType1, const pandora::HitType hitType2,
        const pandora::CartesianVector &fitPosition1, const pandora::CartesianVector &fitPosition2, pandora::CartesianVector &position3D, float &chiSquared) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and a candidate fit position from another views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  hitType the hit type identifying the other view
     *  @param  fitPosition the candidate sliding fit position in the other view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    virtual void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::HitType hitType, const pandora::CartesianVector &fitPosition,
        pandora::CartesianVector &position3D, float &chiSquared) const;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool    m_useChiSquaredApproach;    ///< Whether to obtain y, z positions via chi2 approach, or projected position approach
    bool    m_useDeltaXCorrection;      ///< Whether to add a term to chi squared accounting for hit combination delta x values
    float   m_sigmaX;                   ///< Resolution in x dimension, used for delta x correction to chi squared
};

} // namespace lar_content

#endif // #ifndef LAR_HIT_CREATION_BASE_TOOL_H
