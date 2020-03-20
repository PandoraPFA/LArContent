/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingAlgorithm.h
 *
 *  @brief  Header file for the two view matching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_MATCHING_ALGORITHM_H
#define LAR_TWO_VIEW_MATCHING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.h"

#include "larpandoracontent/LArObjects/LArOverlapMatrix.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  TwoViewMatchingAlgorithm class
 */
template<typename T>
class TwoViewMatchingAlgorithm : public MatchingBaseAlgorithm
{
public:
    typedef OverlapMatrix<T> MatrixType;

    /**
     *  @brief  Default constructor
     */
    TwoViewMatchingAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~TwoViewMatchingAlgorithm();

    void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    const std::string &GetClusterListName(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetInputClusterList(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetSelectedClusterList(const pandora::HitType hitType) const;
    
protected:
    /**
     *  @brief  Calculate cluster overlap result and store in matrix
     *
     *  @param  pClusterU address of view 1 cluster
     *  @param  pClusterV address of view 2 cluster
     */
    virtual void CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) = 0;

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

    const pandora::ClusterList *m_pInputClusterList1;           ///< Address of the input cluster list 1
    const pandora::ClusterList *m_pInputClusterList2;           ///< Address of the input cluster list 2

    pandora::ClusterList        m_clusterList1;                 ///< The selected modified cluster list 1
    pandora::ClusterList        m_clusterList2;                 ///< The selected modified cluster list 2

    MatrixType                  m_overlapMatrix;                ///< The overlap matrix

private:
    typedef std::unordered_map<pandora::HitType, unsigned int, std::hash<int> > HitTypeToIndexMap;
    HitTypeToIndexMap           m_hitTypeToIndexMap;            ///< The hit type to index map

    std::string                 m_inputClusterListName1;        ///< The name of the view 1 cluster list
    std::string                 m_inputClusterListName2;        ///< The name of the view 2 cluster list
};

} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_MATCHING_ALGORITHM_H
