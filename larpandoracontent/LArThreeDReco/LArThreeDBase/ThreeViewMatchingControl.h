/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h
 *
 *  @brief  Header file for the three view matching control class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_MATCHING_CONTROL_H
#define LAR_THREE_VIEW_MATCHING_CONTROL_H 1

#include "larpandoracontent/LArObjects/LArOverlapTensor.h"

namespace lar_content
{

class MatchingBaseAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewMatchingControl class
 */
template<typename T>
class ThreeViewMatchingControl
{
public:
    typedef OverlapTensor<T> TensorType;

    /**
     *  @brief  Constructor
     *
     *  @param  pAlgorithm address of the matching base algorithm
     */
    ThreeViewMatchingControl(MatchingBaseAlgorithm *const pAlgorithm);

    /**
     *  @brief  Destructor
     */
    virtual ~ThreeViewMatchingControl();

    /**
     *  @brief  Get the overlap tensor
     *
     *  @return the overlap tensor
     */
    TensorType &GetOverlapTensor();

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

    const pandora::ClusterList *m_pInputClusterListU;           ///< Address of the input cluster list U
    const pandora::ClusterList *m_pInputClusterListV;           ///< Address of the input cluster list V
    const pandora::ClusterList *m_pInputClusterListW;           ///< Address of the input cluster list W

    pandora::ClusterList        m_clusterListU;                 ///< The selected modified cluster list U
    pandora::ClusterList        m_clusterListV;                 ///< The selected modified cluster list V
    pandora::ClusterList        m_clusterListW;                 ///< The selected modified cluster list W

    TensorType                  m_overlapTensor;                ///< The overlap tensor

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_MATCHING_CONTROL_H
