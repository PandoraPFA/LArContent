/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingContainer.h
 *
 *  @brief  Header file for the two view matching container class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_MATCHING_CONTAINER_H
#define LAR_TWO_VIEW_MATCHING_CONTAINER_H 1

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingContainer.h"

#include "larpandoracontent/LArObjects/LArOverlapMatrix.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  TwoViewMatchingContainer class
 */
template<typename T>
class TwoViewMatchingContainer : public MatchingContainer
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

private:
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
