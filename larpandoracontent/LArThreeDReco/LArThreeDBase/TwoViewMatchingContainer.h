/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingContainer.h
 *
 *  @brief  Header file for the two view matching container class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_MATCHING_CONTAINER_H
#define LAR_TWO_VIEW_MATCHING_CONTAINER_H 1

#include "larpandoracontent/LArObjects/LArOverlapMatrix.h"

#include <unordered_map>

namespace lar_content
{

class MatchingBaseAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TwoViewMatchingContainer class
 */
template<typename T>
class TwoViewMatchingContainer
{
public:
    typedef OverlapMatrix<T> MatrixType;

    /**
     *  @brief  Constructor
     *
     *  @param  pAlgorithm address of the matching base algorithm
     */
    TwoViewMatchingContainer(MatchingBaseAlgorithm *const pAlgorithm);

    /**
     *  @brief  Destructor
     */
    virtual ~TwoViewMatchingContainer();

    /**
     *  @brief  Get the overlap matrix
     *
     *  @return the overlap matrix
     */
    MatrixType &GetOverlapMatrix();

    /**
     *  @brief  Update to reflect addition of a new cluster to the problem space
     *
     *  @param  pNewCluster address of the new cluster
     */
    void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);

    /**
     *  @brief  Update to reflect cluster deletion
     *
     *  @param  pDeletedCluster address of the deleted cluster
     */
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);

    /**
     *  @brief  Get the cluster list name corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the cluster list name
     */
    const std::string &GetClusterListName(const pandora::HitType hitType) const;

    /**
     *  @brief  Get the input cluster list corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the input cluster list
     */
    const pandora::ClusterList &GetInputClusterList(const pandora::HitType hitType) const;

    /**
     *  @brief  Get the selected cluster list corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the selected cluster list
     */
    const pandora::ClusterList &GetSelectedClusterList(const pandora::HitType hitType) const;

    /**
     *  @brief  Get all selected clusters including all hit types
     *
     *  @param  clusterList to receive the list of selected clusters
     */
    void GetAllSelectedClusters(pandora::ClusterList &clusterList) const;

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     */
    void SelectAllInputClusters();

    /**
     *  @brief  Main loop over cluster combinations in order to populate the overlap container. Responsible for calling CalculateOverlapResult.
     */
    void PerformMainLoop();

    /**
     *  @brief  Tidy member variables
     */
    void TidyUp();

    /**
     *  @brief  Read the settings
     *
     *  @param  xmlHandle the relevant xml handle
     */
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    MatchingBaseAlgorithm      *m_pAlgorithm;                   ///< The address of the matching base algorithm

    const pandora::ClusterList *m_pInputClusterList1;           ///< Address of the input cluster list 1
    const pandora::ClusterList *m_pInputClusterList2;           ///< Address of the input cluster list 2

    pandora::ClusterList        m_clusterList1;                 ///< The selected modified cluster list 1
    pandora::ClusterList        m_clusterList2;                 ///< The selected modified cluster list 2

    MatrixType                  m_overlapMatrix;                ///< The overlap matrix

    typedef std::unordered_map<pandora::HitType, unsigned int, std::hash<int> > HitTypeToIndexMap;
    HitTypeToIndexMap           m_hitTypeToIndexMap;            ///< The hit type to index map

    std::string                 m_inputClusterListName1;        ///< The name of the view 1 cluster list
    std::string                 m_inputClusterListName2;        ///< The name of the view 2 cluster list
};

} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_MATCHING_CONTAINER_H
