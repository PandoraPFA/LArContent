/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatching.h
 *
 *  @brief  Header file for the three view delta ray matching class
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/NViewDeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"

namespace lar_content
{

class DeltaRayTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewDeltaRayMatchingAlgorithm class
 */
class ThreeViewDeltaRayMatchingAlgorithm : public NViewDeltaRayMatchingAlgorithm<ThreeViewMatchingControl<DeltaRayOverlapResult>>
{
public:
    typedef NViewDeltaRayMatchingAlgorithm<ThreeViewMatchingControl<DeltaRayOverlapResult>> BaseAlgorithm;
    typedef ThreeViewDeltaRayMatchingAlgorithm::MatchingType::TensorType TensorType;

    /**
     *  @brief  Default constructor
     */
    ThreeViewDeltaRayMatchingAlgorithm();

    /**
     *  @brief  Get the name of the clustering algorithm to be used to recluster created delta ray remnants
     *
     *  @return the clustering algorithm name
     */
    std::string GetClusteringAlgName() const;

private:
    typedef std::vector<DeltaRayTensorTool *> TensorToolVector;

    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);
    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  To check whether a given cluster meets the requirements to be added into the matching container (tensor/matrix)
     *
     *  @param  pCluster the address of the input cluster
     *
     *  @return  whether the checks were met
     */
    virtual bool DoesClusterPassTensorThreshold(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Calculate the overlap result for given group of clusters
     *
     *  @param  pClusterU the cluster from the U view
     *  @param  pClusterV the cluster from the V view
     *  @param  pClusterW the cluster from the W view
     *  @param  overlapResult to receive the overlap result
     *
     *  @return statusCode, faster than throwing in regular use-cases
     */
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
        const pandora::Cluster *const pClusterW, DeltaRayOverlapResult &overlapResult) const;

    /**
     *  @brief  Find the cosmic ray pfos that, in each view, lie close to the clusters of the tensor element   
     *
     *  @param  pClusterU the cluster from the U view
     *  @param  pClusterV the cluster from the V view
     *  @param  pClusterW the cluster from the W view
     *  @param  commonMuonPfoList the output list of common cosmic ray pfos 
     */
    void FindCommonMuonParents(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
        const pandora::Cluster *const pClusterW, pandora::PfoList &commonMuonPfoList) const;

    TensorToolVector m_algorithmToolVector;  ///< The algorithm tool vector
    std::string m_reclusteringAlgorithmName; ///< The name of the clustering algorithm to be used to recluster created delta ray remnants
    unsigned int m_minClusterCaloHits;       ///< The threshold number of hits for a cluster to be considered
    unsigned int m_nMaxTensorToolRepeats;    ///< The maximum number of repeat loops over tensor tools
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DeltaRayTensorTool class
 */
class DeltaRayTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeViewDeltaRayMatchingAlgorithm::MatchingType::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;

    ThreeViewDeltaRayMatchingAlgorithm *m_pParentAlgorithm; ///< Address of the parent matching algorithm
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string ThreeViewDeltaRayMatchingAlgorithm::GetClusteringAlgName() const
{
    return m_reclusteringAlgorithmName;
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
