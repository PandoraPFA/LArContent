/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.cc
 *
 *  @brief  Implementation of the shower hit creation tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPfoHelper.h"

#include "LArThreeDReco/LArHitCreation/ShowerHitsBaseTool.h"
#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerHitsBaseTool::ShowerHitsBaseTool() :
    m_xTolerance(1.f),
    m_chiSquaredCut(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitsBaseTool::Run(ThreeDHitCreationAlgorithm *pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitList &inputTwoDHits,
    CaloHitList &newThreeDHits)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    try
    {
        if (!LArPfoHelper::IsShower(pPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        CaloHitList caloHitListU, caloHitListV, caloHitListW;
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_U, caloHitListU);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_V, caloHitListV);
        pAlgorithm->FilterCaloHitsByType(inputTwoDHits, TPC_VIEW_W, caloHitListW);

        this->CreateThreeDHits(pAlgorithm, caloHitListU, caloHitListV, caloHitListW, newThreeDHits);
        this->CreateThreeDHits(pAlgorithm, caloHitListV, caloHitListU, caloHitListW, newThreeDHits);
        this->CreateThreeDHits(pAlgorithm, caloHitListW, caloHitListU, caloHitListV, newThreeDHits);
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitsBaseTool::CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const CaloHitList &inputTwoDHits,
    const CaloHitList &caloHitList1, const CaloHitList &caloHitList2, CaloHitList &newThreeDHits) const
{
    for (CaloHitList::const_iterator iter = inputTwoDHits.begin(), iterEnd = inputTwoDHits.end(); iter != iterEnd; ++iter)
    {
        try
        {
            CaloHit *pCaloHit2D(*iter);
            const float x(pCaloHit2D->GetPositionVector().GetX());

            CaloHitList filteredList1, filteredList2;
            this->FilterCaloHits(x, m_xTolerance, caloHitList1, filteredList1);
            this->FilterCaloHits(x, m_xTolerance, caloHitList2, filteredList2);

            CartesianVector position3D(0.f, 0.f, 0.f);
            float chiSquared(std::numeric_limits<float>::max());
            this->GetThreeDPosition(pCaloHit2D, filteredList1, filteredList2, position3D, chiSquared);

            if (chiSquared > m_chiSquaredCut)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

            CaloHit *pCaloHit3D(NULL);
            pAlgorithm->CreateThreeDHit(pCaloHit2D, position3D, pCaloHit3D);
            newThreeDHits.insert(pCaloHit3D);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHitsBaseTool::FilterCaloHits(const float x, const float xTolerance, const CaloHitList &inputCaloHitList, CaloHitList &outputCaloHitList) const
{
    outputCaloHitList.clear();

    for (CaloHitList::const_iterator iter = inputCaloHitList.begin(), iterEnd = inputCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;
        const float deltaX(pCaloHit->GetPositionVector().GetX() - x);

        if (std::fabs(deltaX) < xTolerance)
            outputCaloHitList.insert(pCaloHit);
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
