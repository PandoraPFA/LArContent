/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/RecursivePfoMopUpAlgorithm.cc
 *
 *  @brief  Implementation file for the recursive pfo mop up algorithm that runs other algs.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/RecursivePfoMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode RecursivePfoMopUpAlgorithm::Run()
{
    PfoMergeStatsList mergeStatsListBefore(this->GetPfoMergeStats());

    for (unsigned int iter = 0; iter < m_maxIterations; ++iter)
    {
        for (auto const &mopUpAlg : m_mopUpAlgorithms)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, mopUpAlg));

        PfoMergeStatsList mergeStatsListAfter(this->GetPfoMergeStats());

        if (std::equal(mergeStatsListBefore.cbegin(), mergeStatsListBefore.cend(), mergeStatsListAfter.cbegin(), mergeStatsListAfter.cend(), PfoMergeStatsComp))
            break;

        mergeStatsListBefore = std::move(mergeStatsListAfter);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

RecursivePfoMopUpAlgorithm::PfoMergeStatsList RecursivePfoMopUpAlgorithm::GetPfoMergeStats() const
{
    PfoMergeStatsList pfoMergeStatsList;

    for (auto const &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
            continue;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            ClusterNumHitsList pfoHits;
            ClusterList clusterList;
            LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);

            for (auto const &cluster : clusterList)
                pfoHits.push_back(cluster->GetNCaloHits());

            const PropertiesMap &pfoMeta(pPfo->GetPropertiesMap());
            const auto &trackScoreIter(pfoMeta.find("TrackScore"));
            const float trackScore(trackScoreIter != pfoMeta.end() ? trackScoreIter->second : -1.f);

            pfoMergeStatsList.emplace_back(PfoMergeStats{pfoHits, trackScore});
        }
    }
    return pfoMergeStatsList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RecursivePfoMopUpAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "MopUpAlgorithms", m_mopUpAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MaxIterations", m_maxIterations));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
