/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.cc
 *
 *  @brief  Implementation of the shower hit creation tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerHitsBaseTool::ShowerHitsBaseTool() :
    m_xTolerance(1.f),
    m_chiSquaredCut(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitsBaseTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitVector &inputTwoDHits,
    CaloHitVector &newThreeDHits)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    try
    {
        if (!LArPfoHelper::IsShower(pPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        CaloHitVector caloHitVectorU, caloHitVectorV, caloHitVectorW;
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_U, caloHitVectorU);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_V, caloHitVectorV);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_W, caloHitVectorW);

        this->CreateThreeDHits(pAlgorithm, caloHitVectorU, caloHitVectorV, caloHitVectorW, newThreeDHits);
        this->CreateThreeDHits(pAlgorithm, caloHitVectorV, caloHitVectorU, caloHitVectorW, newThreeDHits);
        this->CreateThreeDHits(pAlgorithm, caloHitVectorW, caloHitVectorU, caloHitVectorV, newThreeDHits);
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitsBaseTool::CreateThreeDHits(ThreeDHitCreationAlgorithm *const pAlgorithm, const CaloHitVector &inputTwoDHits,
    const CaloHitVector &caloHitVector1, const CaloHitVector &caloHitVector2, CaloHitVector &newThreeDHits) const
{
    for (const CaloHit *const pCaloHit2D : inputTwoDHits)
    {
        try
        {
            const float x(pCaloHit2D->GetPositionVector().GetX());

            CaloHitVector filteredHits1, filteredHits2;
            this->FilterCaloHits(x, m_xTolerance, caloHitVector1, filteredHits1);
            this->FilterCaloHits(x, m_xTolerance, caloHitVector2, filteredHits2);

            CartesianVector position3D(0.f, 0.f, 0.f);
            float chiSquared(std::numeric_limits<float>::max());
            this->GetThreeDPosition(pCaloHit2D, filteredHits1, filteredHits2, position3D, chiSquared);

            if (chiSquared > m_chiSquaredCut)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

            const CaloHit *pCaloHit3D(NULL);
            pAlgorithm->CreateThreeDHit(pCaloHit2D, position3D, pCaloHit3D);
            newThreeDHits.push_back(pCaloHit3D);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitsBaseTool::FilterCaloHits(const float x, const float xTolerance, const CaloHitVector &inputCaloHitVector, CaloHitVector &outputCaloHitVector) const
{
    for (const CaloHit *const pCaloHit : inputCaloHitVector)
    {
        const float deltaX(pCaloHit->GetPositionVector().GetX() - x);

        if (std::fabs(deltaX) < xTolerance)
            outputCaloHitVector.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerHitsBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "XTolerance", m_xTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    return HitCreationBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
