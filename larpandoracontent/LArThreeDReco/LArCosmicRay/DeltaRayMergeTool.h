/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMergeTool.h
 *
 *  @brief  Header file for the delta ray merge tool class
 *
 *  $Log: $
 */
#ifndef DELTA_RAY_MERGE_TOOL_H
#define DELTA_RAY_MERGE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  DeltaRayMergeTool class
 */
class DeltaRayMergeTool : public DeltaRayTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayMergeTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Identify ambiguous matches (e.g. 3:2:1) and attempt to merge clusters together
     *
     *  @param  overlapTensor the overlap tensor
     *
     *  @return  whether any merges have been made
     */
    bool ExamineConnectedElements(TensorType &overlapTensor) const;

    /**
     *  @brief  Search for two matches with two common clusters and attempt to merge the clusters in the third view together
     *
     *  @param  elementList the tensor element list
     *
     *  @return  whether a merge was made
     */
    bool MakeTwoCommonViewMerges(const TensorType::ElementList &elementList) const;

    /**
     *  @brief  Create a list of the shared common muon pfos of two elements
     *
     *  @param  commonMuonPfoList1 the common muon pfo list of the first element
     *  @param  commonMuonPfoList2 the common muon pfo list of the second element
     *  @param  commonMuonPfoList the output common muon pfo list
     */
    void CombineCommonMuonPfoLists(
        const pandora::PfoList &commonMuonPfoList1, const pandora::PfoList &commonMuonPfoList2, pandora::PfoList &commonMuonPfoList) const;

    /**
     *  @brief  Determine, from a topological point of view, whether two delta ray clusters should be merged together
     *
     *  @param  element1 the first tensor element
     *  @param  element2 the second tensor element
     *  @param  mergeHitType the hit type of the view in which to assess the merge
     *
     *  @return  whether the clusters are topologically associated
     */
    bool AreAssociated(const TensorType::Element &element1, const TensorType::Element &element2, const pandora::HitType &mergeHitType) const;

    /**
     *  @brief  Return the list of muon pfos that a specified delta ray cluster is directly connected to
     *
     *  @param  pDeltaRayCluster the address of the input delta ray cluster
     *  @param  commonMuonPfoList the common muon pfo list of the element to which the DR cluster belongs
     *  @param  connectedMuonPfoList the output list of connected muon pfos
     */
    void GetConnectedMuons(const pandora::Cluster *const pDeltaRayCluster, const pandora::PfoList &commonMuonPfoList,
        pandora::PfoList &connectedMuonPfoList) const;

    /**
     *  @brief  Determine whether a given cluster is connected to a cosmic ray pfo
     *
     *  @param  pCluster the address of the input cluster
     *  @param  pCommonMuonPfo the address of the cosmic ray pfo
     *
     *  @return  whether the cluster is connected to the cosmic ray pfo
     */
    bool IsConnected(const pandora::Cluster *const pCluster, const pandora::Pfo *const pCommonMuonPfo) const;

    /**
     *  @brief  Determine whether two delta ray clusters have been split
     *
     *  @param  pClusterToEnlarge the address of one delta ray cluster
     *  @param  pClusterToDelete the address of the other delta ray cluster
     *
     *  @return  whether the clusters have been split
     */
    bool IsBrokenCluster(const pandora::Cluster *const pClusterToEnlarge, const pandora::Cluster *const pClusterToDelete) const;

    /**
     *  @brief  Determine whether two delta ray clusters are actually a single cluster that is hidden behind a cosmic ray track
     *
     *  @param  pMuonPfo the address of the cosmic ray pfo
     *  @param  pCluster1 the address of one delta ray cluster
     *  @param  pCluster2 the address of the other delta ray cluster
     *
     *  @return  whether the delta ray clusters are one delta ray cluster, hidden behind a cosmic ray track
     */
    bool IsHiddenByTrack(const pandora::ParticleFlowObject *const pMuonPfo, const pandora::Cluster *const pCluster1,
        const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Find all connection points of a delta ray cluster and a cosmic ray pfo
     *
     *  @param  pCommonMuonPfo the address of the cosmic ray pfo
     *  @param  pCluster the address of the delta ray cluster
     *  @param  vertexList the output list of connection points
     */
    void FindVertices(const pandora::Pfo *const pCommonMuonPfo, const pandora::Cluster *const pCluster, pandora::CaloHitList &vertexList) const;

    /**
     *  @brief  Search for two matches with a single common cluster and attempt to merge the clusters in the other two views together
     *
     *  @param  elementList the tensor element list
     *
     *  @return  whether a merge was made
     */
    bool MakeOneCommonViewMerges(const TensorType::ElementList &elementList) const;

    float m_maxDRSeparationFromTrack; ///< The maximum distance of a connected delta ray from a cosmic ray track
    float m_maxClusterSeparation;     ///< The maximum separation of two broken clusters that should be merged
    float m_maxVertexSeparation; ///< The maximum separation of the connection points of two delta ray clusters that are hidden by a CR track and should be merged
    float m_maxGoodMatchReducedChiSquared; ///< The threshold reduced chi squared value for a potential two view merge to go ahead
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_MERGE_TOOL_H
