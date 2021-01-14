/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h
 *
 *  @brief  Header file for the hit creation base tool.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_CREATION_BASE_TOOL_H
#define LAR_HIT_CREATION_BASE_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  HitCreationBaseTool class
 */
class HitCreationBaseTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeDHitCreationAlgorithm::ProtoHit ProtoHit;
    typedef ThreeDHitCreationAlgorithm::ProtoHitVector ProtoHitVector;
    typedef ThreeDHitCreationAlgorithm::TrajectorySample TrajectorySample;

    /**
     *  @brief  Default constructor
     */
    HitCreationBaseTool();

    /**
     *  @brief  Destructor
     */
    virtual ~HitCreationBaseTool();

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo the address of the pfo
     *  @param  inputTwoDHits the vector of input two dimensional hits
     *  @param  protoHitVector to receive the new three dimensional proto hits
     */
    virtual void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector) = 0;

protected:
    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  hitType1 the hit type identifying the first view
     *  @param  hitType2 the hit type identifying the second view
     *  @param  fitPositionList1 the candidate sliding fit position in the first view
     *  @param  fitPositionList2 the candidate sliding fit position in the second view
     *  @param  protoHit to receive the populated proto hit
     */
    virtual void GetBestPosition3D(const pandora::HitType hitType1, const pandora::HitType hitType2, const pandora::CartesianPointVector &fitPositionList1,
        const pandora::CartesianPointVector &fitPositionList2, ProtoHit &protoHit) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  hitType1 the hit type identifying the first view
     *  @param  hitType2 the hit type identifying the second view
     *  @param  fitPosition1 the candidate sliding fit position in the first view
     *  @param  fitPosition2 the candidate sliding fit position in the second view
     *  @param  protoHit to receive the populated proto hit
     */
    virtual void GetBestPosition3D(const pandora::HitType hitType1, const pandora::HitType hitType2, const pandora::CartesianVector &fitPosition1,
        const pandora::CartesianVector &fitPosition2, ProtoHit &protoHit) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and a candidate fit position from another view
     *
     *  @param  hitType the hit type identifying the other view
     *  @param  fitPosition the candidate sliding fit position in the other view
     *  @param  protoHit to receive the populated proto hit
     */
    virtual void GetBestPosition3D(const pandora::HitType hitType, const pandora::CartesianVector &fitPosition, ProtoHit &protoHit) const;

    /**
     *  @brief  Given the current detector volumes and a 3D position, give the distance outside the detector the position is.
     *          Will be 0 for a contained position, otherwise the largest displacement of the Y or Z coord. It is assumed X
     *          is always valid.
     *
     *  @param  larTPCMap a map of LArTPC volumes, to check how contained the given position is.
     *  @param  position3D the current 3D position to check relative to the detector.
     */
    double GetDistanceToDetectorEdge(const pandora::LArTPCMap &larTPCMap, const pandora::CartesianVector &position3D) const;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    double      m_sigmaX2;              ///< The sigmaX squared value, for calculation of chi2 deltaX term
    double      m_sigmaYZ2;             ///< The sigmaYZ squared value, for calculation of chi2 deltaYZ term
    double      m_chiSquaredCut;        ///< The chi squared cut (accept only values below the cut value)
};

} // namespace lar_content

#endif // #ifndef LAR_HIT_CREATION_BASE_TOOL_H
