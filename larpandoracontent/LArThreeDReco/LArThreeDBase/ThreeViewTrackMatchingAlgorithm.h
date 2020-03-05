/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewTrackMatchingAlgorithm.h
 *
 *  @brief  Header file for the three view track matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_TRACK_MATCHING_ALGORITHM_H
#define LAR_THREE_VIEW_TRACK_MATCHING_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ThreeDTransverseTracksAlgorithm class
 */
template<typename T>
class ThreeViewTrackMatchingAlgorithm : public ThreeViewMatchingAlgorithm<T>
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeViewTrackMatchingAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~ThreeViewTrackMatchingAlgorithm();

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

    virtual void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);
    virtual void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    virtual void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;
    virtual void SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const;

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

    /**
     *  @brief  Preparation step for a specific cluster list
     *
     *  @param  clusterList the cluster list
     */
    virtual void PreparationStep(pandora::ClusterList &clusterList);

    virtual void PreparationStep();
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
inline unsigned int ThreeViewTrackMatchingAlgorithm<T>::GetSlidingFitWindow() const
{
    return m_slidingFitWindow;
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_TRACK_MATCHING_ALGORITHM_H
