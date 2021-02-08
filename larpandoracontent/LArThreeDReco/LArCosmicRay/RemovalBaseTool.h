/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/RemovalBase.h
 *
 *  @brief  Header file for the removal base tool class
 *
 *  $Log: $
 */
#ifndef REMOVAL_BASE_TOOL_H
#define REMOVAL_BASE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  RemovalBaseTool class
 */
class RemovalBaseTool : public DeltaRayTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    RemovalBaseTool();

protected:
    virtual bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;    
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) = 0;

    /**
     *  @brief  Return
     *
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  elementList the tensor element list
     *  @param  changesMade receive boolean indicating whether clusters in element have been modified
     */    
    pandora::StatusCode GetMuonCluster(const TensorType::Element &element, const pandora::HitType &hitType, const pandora::Cluster *&pMuonCluster) const;    

    /**
     *  @brief  Determine whether the input element is the best to use to modify the contaminated cluster (best is defined by the total hit count)
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  elementList the tensor element list
     *
     *  @return  whether the input element is the best element
     */
    bool IsBestElement(const TensorType::Element &element, const pandora::HitType &hitType, const TensorType::ElementList &elementList) const;

    
    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart, const pandora::CartesianVector &lineEnd,
        const float distanceToLine) const;     

    void FindExtrapolatedHits(const pandora::Cluster *const pCluster, const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary,
	    pandora::CaloHitList &collectedHits) const;

    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary, const pandora::CartesianVector &point) const;

    pandora::StatusCode ProjectDeltaRayPositions(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
        const pandora::HitType &hitType, pandora::CartesianPointVector &projectedPositions) const;

    float m_distanceToLine;
};

} // namespace lar_content

#endif // #ifndef REMOVAL_BASE_TOOL_H
