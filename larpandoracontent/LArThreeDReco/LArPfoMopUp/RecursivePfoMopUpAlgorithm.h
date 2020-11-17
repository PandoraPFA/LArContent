/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/RecursivePfoMopUpAlgorithm.h
 *
 *  @brief  Header file for the recursive pfo mop up algorithm that runs other algs:
 *          Recusively loop over a series of algorithms until no more changes are made
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
public:
    RecursivePfoMopUpAlgorithm() : m_maxIterations(10){};

private:
    typedef std::vector<unsigned int> ClusterNumHitsList;
    struct PfoMergeStats ///< Object to compare PFO before/after merging algs have run to see if anything changed
    {
        const ClusterNumHitsList m_numClusterHits; ///< Vector filled with number of hits in each of the PFO's clusters
        const float m_trackScore;                  ///< MVA "Track Score" for the PFO
    };
    typedef std::vector<PfoMergeStats> PfoMergeStatsList;

    static bool PfoMergeStatsComp(const PfoMergeStats &lhs, const PfoMergeStats &rhs); ///< Comparator  for PfoMergeStats

    PfoMergeStatsList GetPfoMergeStats() const; ///< Vector filled with PfoMergeStats for each PFP in m_pfoListNames

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_maxIterations;            ///< Maximum number of iterations
    pandora::StringVector m_pfoListNames;    ///< The list of pfo list names
    pandora::StringVector m_mopUpAlgorithms; ///< Ordered list of mop up algorithms to run
};

inline bool RecursivePfoMopUpAlgorithm::PfoMergeStatsComp(
    const RecursivePfoMopUpAlgorithm::PfoMergeStats &lhs, const RecursivePfoMopUpAlgorithm::PfoMergeStats &rhs)
{
    return ((lhs.m_numClusterHits == rhs.m_numClusterHits) && (std::abs(lhs.m_trackScore - rhs.m_trackScore) < std::numeric_limits<float>::epsilon()));
}

} // namespace lar_content

#endif // #ifndef LAR_RECURSIVE_PFO_MOP_UP_ALGORITHM_H
