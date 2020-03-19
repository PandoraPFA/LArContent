/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingContainer.h
 *
 *  @brief  Header file for the matching container class.
 *
 *  $Log: $
 */
#ifndef LAR_MATCHING_CONTAINER_H
#define LAR_MATCHING_CONTAINER_H 1

namespace lar_content
{

class MatchingBaseAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MatchingContainer class
 */
class MatchingContainer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pAlgorithm address of the matching base algorithm
     */
    MatchingContainer(MatchingBaseAlgorithm *const pAlgorithm);

    /**
     *  @brief  Destructor
     */
    virtual ~MatchingContainer();
// TODO look at removing repeated interface
    /**
     *  @brief  Update to reflect addition of a new cluster to the problem space
     *
     *  @param  pNewCluster address of the new cluster
     */
    virtual void UpdateForNewCluster(const pandora::Cluster *const pNewCluster) = 0;

    /**
     *  @brief  Update to reflect cluster deletion
     *
     *  @param  pDeletedCluster address of the deleted cluster
     */
    virtual void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster) = 0;

    /**
     *  @brief  Get the cluster list name corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the cluster list name
     */
    virtual const std::string &GetClusterListName(const pandora::HitType hitType) const = 0;

    /**
     *  @brief  Get the input cluster list corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the input cluster list
     */
    virtual const pandora::ClusterList &GetInputClusterList(const pandora::HitType hitType) const = 0;

    /**
     *  @brief  Get the selected cluster list corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the selected cluster list
     */
    virtual const pandora::ClusterList &GetSelectedClusterList(const pandora::HitType hitType) const = 0;

    /**
     *  @brief  Get the list of all selected clusters for all hit types
     *
     *  @param  clusterList to receive the selected cluster list
     */
    virtual void GetAllSelectedClusters(pandora::ClusterList &clusterList) const = 0;

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     */
    virtual void SelectAllInputClusters() = 0;

    /**
     *  @brief  Main loop over cluster combinations in order to populate the overlap container. Responsible for calling CalculateOverlapResult.
     */
    virtual void PerformMainLoop() = 0;

    /**
     *  @brief  Tidy member variables in derived class
     */
    virtual void TidyUp() = 0;

    /**
     *  @brief  Read settings from xml
     *
     *  @param  xmlHandle the xml handle
     */
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) = 0;

private:
    MatchingBaseAlgorithm  *m_pAlgorithm;   ///< The address of the matching base algorithm
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline MatchingContainer::MatchingContainer(MatchingBaseAlgorithm *const pAlgorithm) :
    m_pAlgorithm(pAlgorithm)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline MatchingContainer::~MatchingContainer()
{
}

} // namespace lar_content

#endif // #ifndef LAR_MATCHING_CONTAINER_H
