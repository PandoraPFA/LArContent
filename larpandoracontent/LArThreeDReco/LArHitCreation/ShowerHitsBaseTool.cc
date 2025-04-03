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
    m_xTolerance(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitsBaseTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfo,
    const CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    try
    {
        if (!LArPfoHelper::IsShower(pPfo))
            return;

        CaloHitVector caloHitVectorU, caloHitVectorV, caloHitVectorW;
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_U, caloHitVectorU);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_V, caloHitVectorV);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_W, caloHitVectorW);

        this->GetShowerHits3D(caloHitVectorU, caloHitVectorV, caloHitVectorW, protoHitVector);
        this->GetShowerHits3D(caloHitVectorV, caloHitVectorU, caloHitVectorW, protoHitVector);
        this->GetShowerHits3D(caloHitVectorW, caloHitVectorU, caloHitVectorV, protoHitVector);
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitsBaseTool::GetShowerHits3D(const CaloHitVector &inputTwoDHits, const CaloHitVector &caloHitVector1,
    const CaloHitVector &caloHitVector2, ProtoHitVector &protoHitVector) const
{
    for (const CaloHit *const pCaloHit2D : inputTwoDHits)
    {
        try
        {
            CaloHitVector filteredHits1, filteredHits2;
            this->FilterCaloHits(pCaloHit2D->GetPositionVector().GetX(), m_xTolerance, caloHitVector1, filteredHits1);
            this->FilterCaloHits(pCaloHit2D->GetPositionVector().GetX(), m_xTolerance, caloHitVector2, filteredHits2);

            ProtoHit protoHit(pCaloHit2D);
            this->GetShowerHit3D(filteredHits1, filteredHits2, protoHit);

            if (protoHit.IsPositionSet() && (protoHit.GetChi2() < m_chiSquaredCut))
                protoHitVector.emplace_back(protoHit);
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
            outputCaloHitVector.emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerHitsBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "XTolerance", m_xTolerance));

    return HitCreationBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
