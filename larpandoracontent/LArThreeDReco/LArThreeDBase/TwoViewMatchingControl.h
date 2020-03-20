/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h
 *
 *  @brief  Header file for the two view matching control class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_MATCHING_CONTROL_H
#define LAR_TWO_VIEW_MATCHING_CONTROL_H 1

#include "larpandoracontent/LArObjects/LArOverlapMatrix.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingControl.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  TwoViewMatchingControl class
 */
template<typename T>
class TwoViewMatchingControl : public NViewMatchingControl
{
public:
    typedef OverlapMatrix<T> MatrixType;

    /**
     *  @brief  Constructor
     *
     *  @param  pAlgorithm address of the matching base algorithm
     */
    TwoViewMatchingControl(MatchingBaseAlgorithm *const pAlgorithm);

    /**
     *  @brief  Destructor
     */
    virtual ~TwoViewMatchingControl();

    /**
     *  @brief  Get the overlap matrix
     *
     *  @return the overlap matrix
     */
    MatrixType &GetOverlapMatrix();

private:
    void UpdateForNewCluster(const pandora::Cluster *const pNewCluster);
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    const std::string &GetClusterListName(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetInputClusterList(const pandora::HitType hitType) const;
    const pandora::ClusterList &GetSelectedClusterList(const pandora::HitType hitType) const;
    void SelectAllInputClusters();
    void PrepareAllInputClusters();
    void PerformMainLoop();
    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const pandora::ClusterList *m_pInputClusterList1;           ///< Address of the input cluster list 1
    const pandora::ClusterList *m_pInputClusterList2;           ///< Address of the input cluster list 2

    pandora::ClusterList        m_clusterList1;                 ///< The selected modified cluster list 1
    pandora::ClusterList        m_clusterList2;                 ///< The selected modified cluster list 2

    MatrixType                  m_overlapMatrix;                ///< The overlap matrix

    typedef std::unordered_map<pandora::HitType, unsigned int, std::hash<int> > HitTypeToIndexMap;
    HitTypeToIndexMap           m_hitTypeToIndexMap;            ///< The hit type to index map

    std::string                 m_inputClusterListName1;        ///< The name of the view 1 cluster list
    std::string                 m_inputClusterListName2;        ///< The name of the view 2 cluster list

    template <typename U>
    friend class NViewMatchingAlgorithm;
};

} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_MATCHING_CONTROL_H
