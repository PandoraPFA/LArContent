/**
 *  @file   larpandoracontent/LArVertex/DirectionFlowProbabilityTool.h
 *
 *  @brief  Header file for the svm vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DIRECTION_FLOW_PROBABILITY_TOOL_H
#define LAR_DIRECTION_FLOW_PROBABILITY_TOOL_H 1

#include "Api/PandoraContentApi.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArControlFlow/TrackDirectionTool.h"

#include<functional>

namespace lar_content
{
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DirectionFlowProbabilityTool class
 */
class DirectionFlowProbabilityTool : public pandora::AlgorithmTool
{
public:

    /**
     *  @brief  Default constructor
     */
    DirectionFlowProbabilityTool();
    
    float GetDirectionFlowProbability(std::function< TrackDirectionTool::DirectionFitObject(const pandora::Cluster* const pCluster) > &directionToolLambda, const pandora::CartesianVector &vertexPosition, pandora::ClusterList &inputClusterVector) const;

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    typedef std::unordered_map<const pandora::Cluster*, pandora::CartesianPointVector> ClusterToSpacepointsMap;
    typedef std::map<const pandora::Cluster*, TrackDirectionTool::DirectionFitObject> DirectionFitMap;
    
    void SelectClusters(std::function< TrackDirectionTool::DirectionFitObject(const pandora::Cluster* const pCluster) > &directionToolLambda, pandora::ClusterList &clusterList, pandora::ClusterVector &selectedClusterVector) const;
    
    pandora::ClusterVector GetPrimaryClusters(const pandora::CartesianVector &positionVector, const pandora::ClusterVector &inputClusterVector) const;
    
    pandora::ClusterVector GetOrderedDaughters(const pandora::CartesianVector &positionVector, const pandora::Cluster* const pParentCluster, const pandora::ClusterVector &inputClusterVector, pandora::ClusterVector &primaryClusters) const;
    
    bool ClusterPointsToPosition(const pandora::Cluster *const pCluster, const pandora::CartesianVector &positionVector, ClusterToSpacepointsMap &clusterToSpacepointsMap) const;
    
    ClusterToSpacepointsMap FillClusterToSpacepointsMap(const pandora::ClusterVector &clusterVector) const;
    
    void GetSpacepoints(const pandora::Cluster *const pCluster, pandora::CartesianPointVector &spacePoints) const;
    
    void AddToSlidingFitCache(std::function< TrackDirectionTool::DirectionFitObject(const pandora::Cluster* const pCluster) > &directionToolLambda, const pandora::Cluster *const pCluster) const;
    
    const TwoDSlidingFitResult &GetCachedSlidingFit(const pandora::Cluster *const pCluster) const;
    
    TrackDirectionTool::DirectionFitObject GetCachedDirectionFit(const pandora::Cluster *const pCluster) const;
    
    unsigned int                    m_slidingFitWindow;                 ///< The layer window for the sliding linear fits
    mutable TwoDSlidingFitResultMap m_slidingFitResultMap;              ///< The sliding fit result map
    mutable DirectionFitMap         m_directionFitMap;                  ///< The direction fit map
    float                           m_impactRadius;                     ///< The impact radius determining whether a sliding fit extrapolation points to a position
    unsigned int                    m_extrapolationNSteps;              ///< The number of steps used in the sliding fit extrapolation method
    float                           m_extrapolationStepSize;            ///< The extrapolation step size.
    float                           m_minimumClusterLength;             ///< The minimum length a cluster must be in order to be considered
};

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_DIRECTION_FLOW_PROBABILITY_TOOL_H
