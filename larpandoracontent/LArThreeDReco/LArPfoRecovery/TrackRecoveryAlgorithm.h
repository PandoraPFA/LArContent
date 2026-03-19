/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.h
 *
 *  @brief  Header file for the track recovery algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_RECOVERY_ALGORITHM_H
#define LAR_TRACK_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  TrackRecoveryAlgorithm class
 */
class TrackRecoveryAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackRecoveryAlgorithm();

private:
    typedef std::unordered_map<pandora::HitType, const pandora::Cluster *> ViewToClusterMap;
    typedef std::unordered_map<const pandora::Pfo *, ViewToClusterMap> PfoToViewClusterMap;
    typedef std::unordered_map<const pandora::Cluster *, pandora::CaloHitList> ClusterToHitMap;
    typedef std::unordered_map<const pandora::CaloHit *, const pandora::Cluster *> HitToClusterToMap;
    typedef std::unordered_map<const pandora::Cluster *, const pandora::Pfo *> ClusterToPfoMap;

    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Performs a sliding linear fit to get an ordered set of hits along a cluster trajectory.
     *
     *  @param  pCluster the cluster for which hits should be ordered
     *  @param  orderedHits the list in which to store the ordered hits
     */
    void FitAndOrderCluster(const pandora::Cluster *const pCluster, pandora::CaloHitList &orderedHits) const;

    std::string m_inputPfoListName; ///< The name of the input PFO list containing track-like PFOs
    pandora::StringVector m_inputClusterListNames; ///< The list of cluster list names to consider when looking for clusters to recover
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_RECOVERY_ALGORITHM_H
