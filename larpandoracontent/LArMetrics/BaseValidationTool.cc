/**
 *  @file   larpandoracontent/LArMetrics/BaseValidationTool.cc
 *
 *  @brief  Implementation of the base validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

using namespace pandora;

namespace lar_content
{

BaseValidationTool::BaseValidationTool() :
    m_maxMichelSep(3.f),    
    m_invalidSmallFloat(-1.f),    
    m_invalidLargeFloat(-9999.f),
    m_invalidInt(-1),    
    m_invalidAngle(-4.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void BaseValidationTool::GetHitsOfType(const CaloHitList &inputList, const HitType hitType, CaloHitVector &outputVector, float &totalEnergy)
{
    totalEnergy = 0.f;

    for (const CaloHit *const pCaloHit : inputList)
    {
        if (pCaloHit->GetHitType() == hitType)
        {
            totalEnergy += pCaloHit->GetElectromagneticEnergy();
            outputVector.push_back(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BaseValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{    
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxMichelSeparation", m_maxMichelSep));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InvalidSmallFloat", m_invalidSmallFloat));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InvalidLargeFloat", m_invalidLargeFloat));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InvalidInt", m_invalidInt));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InvalidAngle", m_invalidAngle));
    
    return STATUS_CODE_SUCCESS;
}

    
} // namespace lar_content
