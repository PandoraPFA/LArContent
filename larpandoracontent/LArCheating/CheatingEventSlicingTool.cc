/**
 *  @file   larpandoracontent/LArCheating/CheatingEventSlicingTool.cc
 *
 *  @brief  Implementation of the cheating event slicing tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingEventSlicingTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventSlicingTool::RunSlicing(const Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
    const HitTypeToNameMap & /*clusterListNames*/, SliceList &sliceList)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    MCParticleToSliceMap mcParticleToSliceMap;
    this->InitializeMCParticleToSliceMap(pAlgorithm, caloHitListNames, mcParticleToSliceMap);

    this->FillSlices(pAlgorithm, TPC_VIEW_U, caloHitListNames, mcParticleToSliceMap);
    this->FillSlices(pAlgorithm, TPC_VIEW_V, caloHitListNames, mcParticleToSliceMap);
    this->FillSlices(pAlgorithm, TPC_VIEW_W, caloHitListNames, mcParticleToSliceMap);

    MCParticleVector mcParticleVector;
    for (const auto &mapEntry : mcParticleToSliceMap)
        mcParticleVector.push_back(mapEntry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleVector)
    {
        const Slice &slice(mcParticleToSliceMap.at(pMCParticle));

        if (!slice.m_caloHitListU.empty() || !slice.m_caloHitListV.empty() || !slice.m_caloHitListW.empty())
            sliceList.push_back(slice);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventSlicingTool::InitializeMCParticleToSliceMap(
    const Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames, MCParticleToSliceMap &mcParticleToSliceMap) const
{
    for (const auto &mapEntry : caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, mapEntry.second, pCaloHitList));

        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            MCParticleVector mcParticleVector;
            for (const auto &weightMapEntry : pCaloHit->GetMCParticleWeightMap())
                mcParticleVector.push_back(weightMapEntry.first);
            std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

            for (const MCParticle *const pMCParticle : mcParticleVector)
            {
                const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

                if (mcParticleToSliceMap.count(pParentMCParticle))
                    continue;

                if (!mcParticleToSliceMap.insert(MCParticleToSliceMap::value_type(pParentMCParticle, Slice())).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventSlicingTool::FillSlices(const Algorithm *const pAlgorithm, const HitType hitType,
    const HitTypeToNameMap &caloHitListNames, MCParticleToSliceMap &mcParticleToSliceMap) const
{
    if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, caloHitListNames.at(hitType), pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMainMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMainMCParticle));

            MCParticleToSliceMap::iterator mapIter = mcParticleToSliceMap.find(pParentMCParticle);

            if (mcParticleToSliceMap.end() == mapIter)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            Slice &slice(mapIter->second);
            CaloHitList &caloHitList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU
                    : (TPC_VIEW_V == hitType)                ? slice.m_caloHitListV
                                                             : slice.m_caloHitListW);
            caloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingEventSlicingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
