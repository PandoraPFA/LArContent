/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/NViewTrackMatchingAlgorithm.h
 *
 *  @brief  Header file for the n view track matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_N_VIEW_TRACK_MATCHING_ALGORITHM_H
#define LAR_N_VIEW_TRACK_MATCHING_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"

namespace lar_content
{

typedef std::unordered_map<const pandora::Cluster*, pandora::CartesianPointVector> SplitPositionMap;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  NViewTrackMatchingAlgorithm class
 */
template<typename T>
class NViewTrackMatchingAlgorithm : public NViewMatchingAlgorithm<T>
{
public:
    /**
     *  @brief  Default constructor
     */
    NViewTrackMatchingAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~NViewTrackMatchingAlgorithm();

    /**
     *  @brief  Get a sliding fit result from the algorithm cache
     *
     *  @param  pCluster address of the relevant cluster
     */
    const TwoDSlidingFitResult &GetCachedSlidingFitResult(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get the layer window for the sliding linear fits
     *
     *  @return the layer window for the sliding linear fits
     */
    unsigned int GetSlidingFitWindow() const;

    /**
     *  @brief  Make cluster splits
     *
     *  @param  splitPositionMap the split position map
     *
     *  @return whether changes to the overlap container have been made
     */
    virtual bool MakeClusterSplits(const SplitPositionMap &splitPositionMap);

    /**
     *  @brief  Make a cluster split
     *
     *  @param  splitPosition the split position
     *  @param  pCurrentCluster the cluster to split
     *  @param  pLowXCluster to receive the low x cluster
     *  @param  pHighXCluster to receive the high x cluster
     *
     *  @return whether a cluster split occurred
     */
    virtual bool MakeClusterSplit(const pandora::CartesianVector &splitPosition, const pandora::Cluster *&pCurrentCluster,
        const pandora::Cluster *&pLowXCluster, const pandora::Cluster *&pHighXCluster) const;

    /**
     *  @brief  Sort split position cartesian vectors by increasing x coordinate
     *
     *  @param  lhs the first cartesian vector
     *  @param  rhs the second cartesian vector
     */
    static bool SortSplitPositions(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs);

    virtual void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);
    virtual void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    virtual void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;
    virtual void PrepareInputClusters(pandora::ClusterList &preparedClusterList);
    virtual void SetPfoParticleId(PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const;

protected:
    /**
     *  @brief  Add a new sliding fit result, for the specified cluster, to the algorithm cache
     *
     *  @param  pCluster address of the relevant cluster
     */
    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Remova an existing sliding fit result, for the specified cluster, from the algorithm cache
     *
     *  @param  pCluster address of the relevant cluster
     */
    void RemoveFromSlidingFitCache(const pandora::Cluster *const pCluster);

    virtual void TidyUp();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    unsigned int                m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap     m_slidingFitResultMap;          ///< The sliding fit result map

    unsigned int                m_minClusterCaloHits;           ///< The min number of hits in base cluster selection method
    float                       m_minClusterLengthSquared;      ///< The min length (squared) in base cluster selection method
};

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline unsigned int NViewTrackMatchingAlgorithm<T>::GetSlidingFitWindow() const
{
    return m_slidingFitWindow;
}

} // namespace lar_content

#endif // #ifndef LAR_N_VIEW_TRACK_MATCHING_ALGORITHM_H
