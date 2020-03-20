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

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingControl.h"

namespace lar_content
{

/**
 *  @brief  ThreeViewMatchingControl class
 */
template<typename T>
class ThreeViewMatchingControl : public NViewMatchingControl
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

private:
    void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    const std::string &GetClusterListName(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetInputClusterList(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetSelectedClusterList(const pandora::HitType hitType) const;
    void GetAllSelectedClusters(pandora::ClusterList &clusterList) const;
    void SelectAllInputClusters();
    void PerformMainLoop();
    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

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

    friend class ThreeViewTrackFragmentsAlgorithm;              ///< ATTN This is for legacy purposes only

    template <typename U>
    friend class NViewMatchingAlgorithm;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_MATCHING_CONTROL_H
