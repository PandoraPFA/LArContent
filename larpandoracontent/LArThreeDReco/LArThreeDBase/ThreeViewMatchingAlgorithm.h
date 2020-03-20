/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingAlgorithm.h
 *
 *  @brief  Header file for the three view matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_MATCHING_ALGORITHM_H
#define LAR_THREE_VIEW_MATCHING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.h"

#include "larpandoracontent/LArObjects/LArOverlapTensor.h"

namespace lar_content
{

/**
 *  @brief  ThreeViewMatchingAlgorithm class
 */
template<typename T>
class ThreeViewMatchingAlgorithm : public MatchingBaseAlgorithm
{
public:
    typedef OverlapTensor<T> TensorType;

    /**
     *  @brief  Default constructor
     */
    ThreeViewMatchingAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~ThreeViewMatchingAlgorithm();

    void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    const std::string &GetClusterListName(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetInputClusterList(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetSelectedClusterList(const pandora::HitType hitType) const;
    
protected:
    /**
     *  @brief  Calculate cluster overlap result and store in tensor
     *
     *  @param  pClusterU address of U view cluster
     *  @param  pClusterV address of V view cluster
     *  @param  pClusterW address of W view cluster
     */
    virtual void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW) = 0;

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     *
     *  @param  pInputClusterList address of an input cluster list
     *  @param  selectedClusterList to receive the selected cluster list
     */
    virtual void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const = 0;

    virtual void SelectAllInputClusters();
    virtual void PerformMainLoop();
    virtual void TidyUp();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const pandora::ClusterList *m_pInputClusterListU;           ///< Address of the input cluster list U
    const pandora::ClusterList *m_pInputClusterListV;           ///< Address of the input cluster list V
    const pandora::ClusterList *m_pInputClusterListW;           ///< Address of the input cluster list W

    pandora::ClusterList        m_clusterListU;                 ///< The selected modified cluster list U
    pandora::ClusterList        m_clusterListV;                 ///< The selected modified cluster list V
    pandora::ClusterList        m_clusterListW;                 ///< The selected modified cluster list W

    TensorType                  m_overlapTensor;                ///< The overlap tensor

private:
    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_MATCHING_ALGORITHM_H
