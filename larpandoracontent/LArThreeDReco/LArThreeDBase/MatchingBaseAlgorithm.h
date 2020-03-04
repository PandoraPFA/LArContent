/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.h
 *
 *  @brief  Header file for the three dimension algorithm base class.
 *
 *  $Log: $
 */
#ifndef LAR_MATCHING_BASE_ALGORITHM_H
#define LAR_MATCHING_BASE_ALGORITHM_H 1

#include "Api/PandoraContentApi.h"

#include "Pandora/Algorithm.h"

#include <unordered_map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  ProtoParticle
 */
class ProtoParticle
{
public:
    pandora::ClusterList    m_clusterListU;                 ///< List of 2D U clusters in a 3D proto particle
    pandora::ClusterList    m_clusterListV;                 ///< List of 2D V clusters in a 3D proto particle
    pandora::ClusterList    m_clusterListW;                 ///< List of 2D W clusters in a 3D proto particle
};

typedef std::vector<ProtoParticle> ProtoParticleVector;
typedef std::unordered_map<const pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;

//------------------------------------------------------------------------------------------------------------------------------------------

typedef std::unordered_map<const pandora::Cluster*, pandora::CartesianPointVector> SplitPositionMap;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MatchingBaseAlgorithm class
 */

class MatchingBaseAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MatchingBaseAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~MatchingBaseAlgorithm();

    /**
     *  @brief  Update to reflect addition of a new cluster to the problem space
     *
     *  @param  pNewCluster address of the new cluster
     */
    virtual void UpdateForNewCluster(const pandora::Cluster *const pNewCluster) = 0;

    /**
     *  @brief  Update to reflect cluster deletion
     *
     *  @param  pDeletedCluster address of the deleted cluster
     */
    virtual void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster) = 0;

    /**
     *  @brief  Get the cluster list name corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the cluster list name
     */
    virtual const std::string &GetClusterListName(const pandora::HitType hitType) const = 0;

    /**
     *  @brief  Get the input cluster list corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the input cluster list
     */
    virtual const pandora::ClusterList &GetInputClusterList(const pandora::HitType hitType) const = 0;

    /**
     *  @brief  Get the selected cluster list corresponding to a specified hit type
     *
     *  @param  hitType the hit type
     *
     *  @return the selected cluster list
     */
    virtual const pandora::ClusterList &GetSelectedClusterList(const pandora::HitType hitType) const = 0;

    /**
     *  @brief  Create particles using findings from recent algorithm processing
     *
     *  @param  protoParticleVector the proto particle vector
     *
     *  @param  whether particles were created
     */
    virtual bool CreateThreeDParticles(const ProtoParticleVector &protoParticleVector);

    /**
     *  @brief  Calculate Pfo properties from proto particle
     *
     *  @param  protoParticle the input proto particle
     *  @param  pfoParameters the output pfo parameters
     */
    virtual void SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const;

    /**
     *  @brief  Merge clusters together
     *
     *  @param  clusterMergeMap the cluster merge map
     *
     *  @return whether changes to the tensor have been made
     */
    virtual bool MakeClusterMerges(const ClusterMergeMap &clusterMergeMap);

    /**
     *  @brief  Make cluster splits
     *
     *  @param  splitPositionMap the split position map
     *
     *  @return whether changes to the tensor have been made
     */
    virtual bool MakeClusterSplits(const SplitPositionMap &splitPositionMap);

    /**
     *  @brief  Make a cluster split
     *
     *  @param  splitPosition the split position
     *  @param  pCurrentCluster the cluster to split
     *  @param  pLowXCluster to receive the low x cluster
     *  @param  pHighXCluster to receive the high x cluster
     *
     *  @return whether a cluster split occurred
     */
    virtual bool MakeClusterSplit(const pandora::CartesianVector &splitPosition, const pandora::Cluster *&pCurrentCluster,
        const pandora::Cluster *&pLowXCluster, const pandora::Cluster *&pHighXCluster) const;

    /**
     *  @brief  Sort split position cartesian vectors by increasing x coordinate
     *
     *  @param  lhs the first cartesian vector
     *  @param  rhs the second cartesian vector
     */
    static bool SortSplitPositions(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs);

protected:
    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     */
    virtual void SelectAllInputClusters() = 0;

    /**
     *  @brief  Perform any preparatory steps required, e.g. caching expensive fit results for clusters
     */
    virtual void PreparationStep();

    /**
     *  @brief  Main loop over cluster combinations in order to populate the tensor. Responsible for calling CalculateOverlapResult.
     */
    virtual void PerformMainLoop() = 0;

    /**
     *  @brief  Examine contents of tensor, collect together best-matching 2D particles and modify clusters as required
     */
    virtual void ExamineTensor() = 0;

    /**
     *  @brief  Tidy member variables in derived class
     */
    virtual void TidyUp();

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    pandora::StatusCode Run();

    std::string                 m_outputPfoListName;            ///< The output pfo list name
};

} // namespace lar_content

#endif // #ifndef LAR_MATCHING_BASE_ALGORITHM_H
