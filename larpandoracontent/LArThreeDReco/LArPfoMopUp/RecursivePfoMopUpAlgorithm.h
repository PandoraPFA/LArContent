/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/RecursivePfoMopUpAlgorithm.h
 *
 *  @brief  Header file for the mega cluster mop up algorithm that runs other algs.
 *
 *  $Log: $
 */
#ifndef LAR_RECURSIVE_PFO_MOP_UP_ALGORITHM_H
#define LAR_RECURSIVE_PFO_MOP_UP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content {

/**
 *  @brief  RecursivePfoMopUpAlgorithm class
 */
class RecursivePfoMopUpAlgorithm : public pandora::Algorithm {
  private:
  // Struct to determine if a pfo has been merged
  struct pfoMergeStats {
    const std::vector<unsigned int> mNumHits;
    const float mTrackScore;
  };

  static bool pfoMergeStatsComp(const pfoMergeStats& lhs, const pfoMergeStats& rhs)
  {
    return ((lhs.mNumHits == rhs.mNumHits) && ((lhs.mTrackScore - rhs.mTrackScore) < std::numeric_limits<float>::epsilon()));
  };

  pandora::StatusCode Run();
  std::vector<pfoMergeStats> GetPfoMergeStats();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  pandora::StringVector m_pfoListNames;    ///< The list of pfo list names
  pandora::StringVector m_mopUpAlgorithms; ///< Ordered list of mop up algorithms to run
};

} // namespace lar_content

#endif // #ifndef LAR_RECURSIVE_PFO_MOP_UP_ALGORITHM_H
