/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/RecursivePfoMopUpAlgorithm.h
 *
 *  @brief  Header file for the recursive pfo mop up algorithm that runs other algs:
 *          Recusively loop over a series of algortihms until no more changes are made
 *
 *  $Log: $
 */
#ifndef LAR_RECURSIVE_PFO_MOP_UP_ALGORITHM_H
#define LAR_RECURSIVE_PFO_MOP_UP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  RecursivePfoMopUpAlgorithm class
 */
class RecursivePfoMopUpAlgorithm : public pandora::Algorithm
{
private:
    struct pfoMergeStats
    {                                             ///< Object to compare PFO before/after merging algs have run to see if anything changed
        const std::vector<unsigned int> mNumHits; ///< Vector filled with number of hits in each of the PFO's clusters
        const float mTrackScore;                  ///< MVA "Track Score" for the PFO
    };

    static bool pfoMergeStatsComp(const pfoMergeStats &lhs, const pfoMergeStats &rhs) ///< Comparitor  for pfoMergeStats
    {
        return ((lhs.mNumHits == rhs.mNumHits) && (std::abs(lhs.mTrackScore - rhs.mTrackScore) < std::numeric_limits<float>::epsilon()));
    };

    pandora::StatusCode Run();
    std::vector<pfoMergeStats> GetPfoMergeStats() const; ///< Vector filled with pfoMergeStats for each PFP in m_pfoListNames
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_maxIterations;            ///< Maximum number of iterations
    pandora::StringVector m_pfoListNames;    ///< The list of pfo list names
    pandora::StringVector m_mopUpAlgorithms; ///< Ordered list of mop up algorithms to run
};

} // namespace lar_content

#endif // #ifndef LAR_RECURSIVE_PFO_MOP_UP_ALGORITHM_H
