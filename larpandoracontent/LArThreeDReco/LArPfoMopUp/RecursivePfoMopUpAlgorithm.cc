/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/RecursivePfoMopUpAlgorithm.h
 *
 *  @brief  Implementation file for the recursive pfo mop up algorithm that runs other algs.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/RecursivePfoMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content {

StatusCode RecursivePfoMopUpAlgorithm::Run()
{
  bool unchanged(false);
  std::vector<RecursivePfoMopUpAlgorithm::pfoMergeStats> mergeStatsVecBefore(GetPfoMergeStats());

  while (!unchanged) {

    for (auto const& mopUpAlg : m_mopUpAlgorithms) {
      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, mopUpAlg));
    }

    std::vector<RecursivePfoMopUpAlgorithm::pfoMergeStats> mergeStatsVecAfter(GetPfoMergeStats());

    unchanged = std::equal(mergeStatsVecBefore.cbegin(), mergeStatsVecBefore.cend(), mergeStatsVecAfter.cbegin(),
        mergeStatsVecAfter.cend(), RecursivePfoMopUpAlgorithm::pfoMergeStatsComp);

    mergeStatsVecBefore = std::move(mergeStatsVecAfter);
  } // while !unchanged

  return STATUS_CODE_SUCCESS;
} // Run

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<RecursivePfoMopUpAlgorithm::pfoMergeStats> RecursivePfoMopUpAlgorithm::GetPfoMergeStats() const
{
  std::vector<RecursivePfoMopUpAlgorithm::pfoMergeStats> pfoMergeStatsVec;

  for (auto const& pfoListName : m_pfoListNames) {

    const PfoList* pPfoList = NULL;
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
      continue;

    for (const ParticleFlowObject* const pPfo : *pPfoList) {

      std::vector<unsigned int> pfoHits;
      ClusterList clusterList;
      LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
      for (auto const& cluster : clusterList) {
        pfoHits.push_back(cluster->GetNCaloHits());
      }

      const PropertiesMap pfoMeta(pPfo->GetPropertiesMap());
      const auto& trackScoreIter = pfoMeta.find("TrackScore");
      const float trackScore(trackScoreIter != pfoMeta.end() ? trackScoreIter->second : -1.f);

      pfoMergeStatsVec.push_back(RecursivePfoMopUpAlgorithm::pfoMergeStats{ pfoHits, trackScore });
    } // pPfo : pPfoList
  }   // pfoListName : m_pfoListNames
  return pfoMergeStatsVec;
} // GetPfoMergeStats

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RecursivePfoMopUpAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
      XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "MopUpAlgorithms", m_mopUpAlgorithms));

  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
      XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

  return STATUS_CODE_SUCCESS;
} // ReadSettings

} // namespace lar_content
