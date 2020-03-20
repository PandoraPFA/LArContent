/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h
 *
 *  @brief  Header file for the n view matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_N_VIEW_MATCHING_ALGORITHM_H
#define LAR_N_VIEW_MATCHING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  NViewMatchingAlgorithm class
 */
template<typename T>
class NViewMatchingAlgorithm : public MatchingBaseAlgorithm
{
public:
    typedef T ContainerType;

    /**
     *  @brief  Default constructor
     */
    NViewMatchingAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~NViewMatchingAlgorithm();

    void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    const std::string &GetClusterListName(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetInputClusterList(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetSelectedClusterList(const pandora::HitType hitType) const;

protected:
    /**
     *  @brief  Get the matching container
     */
    ContainerType &GetMatchingContainer();

    virtual void SelectAllInputClusters();
    virtual void PerformMainLoop();
    virtual void TidyUp();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    ContainerType   m_matchingContainer;     ///< The matching container
};

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline T &NViewMatchingAlgorithm<T>::GetMatchingContainer()
{
    return m_matchingContainer;
}

} // namespace lar_content

#endif // #ifndef LAR_N_VIEW_MATCHING_ALGORITHM_H
