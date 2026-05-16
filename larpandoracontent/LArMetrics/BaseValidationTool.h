/**
 *  @file   larpandoracontent/LArMetrics/BaseValidationTool.h
 *
 *  @brief  Header file for the base validation tool class.
 *
 *  $Log: $
 */
#ifndef BASE_VALIDATION_TOOL_H
#define BASE_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

namespace lar_content
{

/**
 *  @brief  BaseValidationTool class
 */
class BaseValidationTool : public pandora::AlgorithmTool
{
public:

    /**
     *  @brief  Default constructor
     */
    BaseValidationTool();

    virtual pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch) = 0;

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Pick out the hits of a specified type from an input list of hits
     *
     *  @param[in]   inputList the input hit list
     *  @param[in]   hitType the hit type
     *  @param[out]  outputVector the output vector of hits
     *  @param[out]  totalEnergy the total charge of the filtered hits
     */
    void GetHitsOfType(const pandora::CaloHitList &inputList, const pandora::HitType hitType, pandora::CaloHitVector &outputVector, 
        float &totalEnergy);

    float m_maxMichelSep;      ///< the maximum separation between a michel vertex and parent muon endpoint    
    float m_invalidSmallFloat; ///< small float value used for failure case (visible on plots)
    float m_invalidLargeFloat; ///< large float value used for failure case (used when small value is a valid value)
    int m_invalidInt;          ///< int value used for failure case    
    float m_invalidAngle;      ///< angle value used for failure case
};

} // namespace lar_content

#endif // #ifndef BASE_VALIDATION_TOOL_H
