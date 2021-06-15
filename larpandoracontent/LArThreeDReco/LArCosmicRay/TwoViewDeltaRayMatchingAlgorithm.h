/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatching.h
 *
 *  @brief  Header file for the two view delta ray matching class
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/NViewDeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{

class DeltaRayMatrixTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TwoViewDeltaRayMatchingAlgorithm class
 */
class TwoViewDeltaRayMatchingAlgorithm : public NViewDeltaRayMatchingAlgorithm<TwoViewMatchingControl<TwoViewDeltaRayOverlapResult>>
{
public:
    typedef NViewDeltaRayMatchingAlgorithm<TwoViewMatchingControl<TwoViewDeltaRayOverlapResult>> BaseAlgorithm;
    typedef TwoViewDeltaRayMatchingAlgorithm::MatchingType::MatrixType MatrixType;

    /**
     *  @brief  Default constructor
     */
    TwoViewDeltaRayMatchingAlgorithm();

    /**
     *  @brief  Get the name of the third view clusters
     *
     *  @return the third view cluster list name
     */
    const std::string &GetThirdViewClusterListName() const;

    /**
     *  @brief  Get the name of the clustering algorithm to be used to recluster created delta ray remnants
     *
     *  @return the clustering algorithm name
     */
    const std::string &GetClusteringAlgName() const;

    /**
     *  @brief  Obtain the HitTypeVector of input views
     *
     *  @return  the HitTypeVector of input views
     */
    HitTypeVector GetHitTypeVector();

    /**
     *  @brief  Get the address of the given hit type cluster
     *
     *  @param  hitType hit type of the required cluster
     *
     *  @return address of the required cluster
     */
    const pandora::Cluster *GetCluster(const MatrixType::Element &element, const pandora::HitType hitType);

    /**
     *  @brief  Create delta ray pfos out of a given element, merging the third view clusters together and adding in any stray clusters
     *
     *  @param  element the matrix element
     *
     *  @return  whether any pfos were created
     */
    bool CreatePfo(const MatrixType::Element &element);

    /**
     *  @brief  Update the matrix after a third view cluster modification - remove delta ray clusters and reassess the matching of cosmic ray clusters
     *
     *  @param  pModifiedCluster the address of the modified cluster
     *  @param  isMuon whether the modified cluster belongs to a cosmic ray pfo
     */
    void UpdateForThirdViewClusterModification(const pandora::Cluster *const pModifiedCluster, const bool isMuon);

private:
    void CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3);

    /**
     *  @brief  To check whether a given cluster meets the requirements to be added into the matching container (tensor/matrix)
     *
     *  @param  pCluster the address of the input cluster
     *
     *  @return  whether the checks were met
     */
    virtual bool DoesClusterPassTensorThreshold(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Calculate the overlap result for given pair of clusters
     *
     *  @param  pCluster1 the cluster from the first input view
     *  @param  pCluster2 the cluster from the second input view
     *  @param  overlapResult to receive the overlap result
     *
     *  @return statusCode, faster than throwing in regular use-cases
     */
    pandora::StatusCode CalculateOverlapResult(
        const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, TwoViewDeltaRayOverlapResult &overlapResult) const;

    /**
     *  @brief  Find the cosmic ray pfos that, in each view, lie close to the clusters of the matrix element
     *
     *  @param  pCluster1 the cluster from the first input view
     *  @param  pCluster2 the cluster from the second input view
     *  @param  commonMuonPfoList the output list of common cosmic ray pfos
     */
    void FindCommonMuonParents(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, pandora::PfoList &commonMuonPfoList) const;

    /**
     *  @brief  Collect the available and unavailable third view clusters that lie close to the projected delta ray hits
     *
     *  @param  pCluster1 the cluster from the first input view
     *  @param  pCluster2 the cluster from the second input view
     *  @param  projectedPositions the projected positions of the matched cluster pair
     *  @param  matchedClusters the output list of collected clusters
     */
    void CollectThirdViewClusters(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const pandora::CartesianPointVector &projectedPositions, pandora::ClusterList &matchedClusters) const;

    /**
     *  @brief  Determine the best matched third view cluster and calculate the reduced chi-squared value of the three cluster match
     *
     *  @param  pCluster1 the cluster from the first input view
     *  @param  pCluster2 the cluster from the second input view
     *  @param  commonMuonPfoList the list of common cosmic ray pfos
     *  @param  matchedClusters the list of third view matched clusters
     *  @param  reducedChiSquared to receive the calculated reduced chi-squared value
     *
     *  @return  the address of the best matched third view cluster
     */
    const pandora::Cluster *GetBestMatchedCluster(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const pandora::PfoList &commonMuonPfoList, const pandora::ClusterList &matchedClusters, float &reducedChiSquared) const;

    /**
     *  @brief  Form the third view cluster by removing hits from cosmic ray clusters and merging the matched clusters where appropriate
     *
     *  @param  element the matrix element
     *  @param  protoParticle the output proto particle
     */
    void FormThirdViewCluster(const MatrixType::Element &element, ProtoParticle &protoParticle);

    /**
     *  @brief  Starting with an input seed cluster, sequentially merge in matched clusters that retain a good reduced chi-squared
     *
     *  @param  element the matrix element
     *  @param  pSeedCluster the address of the input seed cluster
     */
    void MergeThirdView(const MatrixType::Element &element, const pandora::Cluster *const pSeedCluster);

    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<DeltaRayMatrixTool *> MatrixToolVector;

    std::string m_inputClusterListName;      ///< The name of the cluster list in the view in which to project into
    std::string m_reclusteringAlgorithmName; ///< The name of the clustering algorithm to be used to recluster created delta ray remnants
    MatrixToolVector m_algorithmToolVector;  ///< The algorithm tool vector
    unsigned int m_nMaxMatrixToolRepeats;    ///< The maximum number of repeat loops over matrix tools
    unsigned int m_minClusterCaloHits;       ///< The threshold number of hits for a cluster to be considered
    float m_maxDistanceFromPrediction;       ///< The maximum distance of a matched cluster from the third view projection points
    float m_maxGoodMatchReducedChiSquared;   ///< The maximum reduced chi squared value of a good 1:1:1 match
    float m_minDistanceFromMuon;             ///< The minimum distance of a hit from the cosmic ray track required for removal
    float m_maxDistanceToCollected;          ///< The maximim distance of a hit from the projected delta ray hits required for removal
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &TwoViewDeltaRayMatchingAlgorithm::GetThirdViewClusterListName() const
{
    return m_inputClusterListName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &TwoViewDeltaRayMatchingAlgorithm::GetClusteringAlgName() const
{
    return m_reclusteringAlgorithmName;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DeltaRayTensorTool class
 */
class DeltaRayMatrixTool : public pandora::AlgorithmTool
{
public:
    typedef TwoViewDeltaRayMatchingAlgorithm::MatchingType::MatrixType MatrixType;
    typedef std::vector<MatrixType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &matrixTensor) = 0;

    TwoViewDeltaRayMatchingAlgorithm *m_pParentAlgorithm; ///< Address of the parent matching algorithm
};

} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
