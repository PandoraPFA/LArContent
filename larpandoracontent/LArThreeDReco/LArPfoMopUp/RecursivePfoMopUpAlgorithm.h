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
    /**
     *  @brief  Default constructor
     */
    RecursivePfoMopUpAlgorithm();

private:
    typedef std::vector<unsigned int> ClusterNumHitsList;

    /**
     *  @brief  PfoMergeStats class: Object to compare PFO before/after merging algs have run to see if anything changed
     */
    struct PfoMergeStats
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  Vector filled with number of hits in each of the PFO's clusters
         *  @param  MVA "Track Score" for the PFO
         */
        PfoMergeStats(const ClusterNumHitsList &numClusterHits, const float trackScore);

        const ClusterNumHitsList m_numClusterHits; ///< Vector filled with number of hits in each of the PFO's clusters
        const float m_trackScore;                  ///< MVA "Track Score" for the PFO
    };

    typedef std::vector<PfoMergeStats> PfoMergeStatsList;

    /**
     *  @brief  Equality comparator for two PfoMergeStats
     *
     *  @param  lhs PfoMergeStats for comparison
     *  @param  rhs PfoMergeStats for comparison
     *
     *  @return boolean if the lhs and rhs are the same
     */
    static bool PfoMergeStatsComp(const PfoMergeStats &lhs, const PfoMergeStats &rhs);

    /**
     *  @brief  Get the PfoMergeStats for all of the particles in the event from m_pfoListNames
     *
     *  @return List of PfoMergeStats for each Pfo
     */
    PfoMergeStatsList GetPfoMergeStats() const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_maxIterations;            ///< Maximum number of iterations
    pandora::StringVector m_pfoListNames;    ///< The list of pfo list names
    pandora::StringVector m_mopUpAlgorithms; ///< Ordered list of mop up algorithms to run
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline RecursivePfoMopUpAlgorithm::RecursivePfoMopUpAlgorithm() : m_maxIterations(10)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline RecursivePfoMopUpAlgorithm::PfoMergeStats::PfoMergeStats(const ClusterNumHitsList &numClusterHits, const float trackScore) :
    m_numClusterHits(numClusterHits), m_trackScore(trackScore)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool RecursivePfoMopUpAlgorithm::PfoMergeStatsComp(
    const RecursivePfoMopUpAlgorithm::PfoMergeStats &lhs, const RecursivePfoMopUpAlgorithm::PfoMergeStats &rhs)
{
    return ((lhs.m_numClusterHits == rhs.m_numClusterHits) && (std::abs(lhs.m_trackScore - rhs.m_trackScore) < std::numeric_limits<float>::epsilon()));
}

} // namespace lar_content

#endif // #ifndef LAR_RECURSIVE_PFO_MOP_UP_ALGORITHM_H
