/**
 *  @file   LArContent/src/LArCheating/CheatingEventSlicingTool.cc
 * 
 *  @brief  Implementation of the cheating event slicing tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent/LArCheating/CheatingEventSlicingTool.h"

#include "larpandoracontent/LArContent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoParentAlgorithm::HitTypeToNameMap HitTypeToNameMap;

void CheatingEventSlicingTool::Slice(const NeutrinoParentAlgorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
    const HitTypeToNameMap &/*clusterListNames*/, NeutrinoParentAlgorithm::SliceList &sliceList)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    const MCParticleList *pMCParticleList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, m_mcParticleListName, pMCParticleList));

    MCParticleToSliceMap mcParticleToSliceMap;

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

        if (mcParticleToSliceMap.count(pParentMCParticle))
            continue;

        if (!mcParticleToSliceMap.insert(MCParticleToSliceMap::value_type(pParentMCParticle, NeutrinoParentAlgorithm::Slice())).second)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    this->FillSlices(pAlgorithm, TPC_VIEW_U, caloHitListNames, mcParticleToSliceMap);
    this->FillSlices(pAlgorithm, TPC_VIEW_V, caloHitListNames, mcParticleToSliceMap);
    this->FillSlices(pAlgorithm, TPC_VIEW_W, caloHitListNames, mcParticleToSliceMap);

    for (const MCParticleToSliceMap::value_type &value : mcParticleToSliceMap)
    {
        const NeutrinoParentAlgorithm::Slice &slice(value.second);

        if (!slice.m_caloHitListU.empty() || !slice.m_caloHitListV.empty() || !slice.m_caloHitListW.empty())
            sliceList.push_back(slice);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventSlicingTool::FillSlices(const NeutrinoParentAlgorithm *const pAlgorithm, const HitType hitType, const HitTypeToNameMap &caloHitListNames,
    MCParticleToSliceMap &mcParticleToSliceMap) const
{
    if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, caloHitListNames.at(hitType), pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            if (hitType != pCaloHit->GetHitType())
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const MCParticle *const pMainMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMainMCParticle));

            MCParticleToSliceMap::iterator mapIter = mcParticleToSliceMap.find(pParentMCParticle);

            if (mcParticleToSliceMap.end() == mapIter)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            NeutrinoParentAlgorithm::Slice &slice(mapIter->second);
            CaloHitList &caloHitList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU : (TPC_VIEW_V == hitType) ? slice.m_caloHitListV : slice.m_caloHitListW);
            caloHitList.insert(pCaloHit);
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingEventSlicingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
