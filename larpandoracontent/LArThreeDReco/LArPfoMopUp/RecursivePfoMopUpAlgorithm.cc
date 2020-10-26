/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/RecursivePfoMopUpAlgorithm.h
 *
 *  @brief  Implementation file for the mega cluster mop up algorithm that runs other algs.
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

    unchanged = std::equal(mergeStatsVecBefore.cbegin(), mergeStatsVecBefore.cend(), mergeStatsVecAfter.cbegin(), mergeStatsVecAfter.cend(), RecursivePfoMopUpAlgorithm::pfoMergeStatsComp);

    mergeStatsVecBefore = std::move(mergeStatsVecAfter);
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<RecursivePfoMopUpAlgorithm::pfoMergeStats> RecursivePfoMopUpAlgorithm::GetPfoMergeStats()
{
  std::vector<RecursivePfoMopUpAlgorithm::pfoMergeStats> pfoMergeStatsVec;

  for (StringVector::const_iterator sIter = m_pfoListNames.begin(), sIterEnd = m_pfoListNames.end(); sIter != sIterEnd; ++sIter) {
    const PfoList* pPfoList = NULL;
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *sIter, pPfoList))
      continue;

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter) {
      const ParticleFlowObject* const pPfo = *pIter;

      unsigned int pfoHits(0);
      ClusterList clusterList;
      LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
      for (auto const& cluster : clusterList) {
        pfoHits += cluster->GetNCaloHits();
      }
      // TODO get TrackScore Metadata
      float trackScore(-1.f);

      pfoMergeStatsVec.push_back(RecursivePfoMopUpAlgorithm::pfoMergeStats{ pfoHits, clusterList.size(), trackScore });
    }
  }
  return pfoMergeStatsVec;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RecursivePfoMopUpAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "MopUpAlgorithms", m_mopUpAlgorithms));

  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

  return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
