/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPBaseHierarchyTool.h
 *
 *  @brief  Header file for the MLP base hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_MLP_BASE_HIERARCHY_TOOL_H
#define LAR_MLP_BASE_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace lar_dl_content
{

/**
 *   @brief  MLPBaseHierarchyTool to calculate variables related to the initial shower region
 */
class MLPBaseHierarchyTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    MLPBaseHierarchyTool();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
  
    /**
     *  @brief  Set the detector boundaries 
     */
    void SetDetectorBoundaries();

    /**
     *  @brief  Return whether and input point is within the bounds of the detector
     *
     *  @param  position the input CartesianVector
     *
     *  @return Whether the input position is within the detector
     */
    bool IsInFV(const pandora::CartesianVector &position) const;

    /**
     *  @brief  Get the number of 3D hits owned by a pfo
     *
     *  @param  hierarchyPfo the HierarchyPfo of the corresponding pfo
     *
     *  @return The number of 3D hits
     */
    float GetNSpacepoints(const HierarchyPfo &hierarchyPfo) const;

    /**
     *  @brief  Return the number of 3D hits and the number of corresponding 
     *          pfos of a given pfo about a point
     *
     *  @param  pAlgorithm a pointer to the pandora algorithm calling the tool
     *  @param  pPfo a pointer to the pfo
     *  @param  pointOfInterest the input position
     *
     *  @return std::pair of the number of 3D hits and the number of corresponding pfos
     */
    std::pair<float, float> GetParticleInfoAboutPfoPosition(const pandora::Algorithm *const pAlgorithm, 
        const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &pointOfInterest) const;

    /**
     *  @brief  Shift and normalise a network parameter with respect to an input range
     *
     *  @param  minLimit the minimum allowed value of the variable
     *  @param  maxLimit the maximum allowed value of the variable
     *  @param  networkParam the input network parameter value
     */
    void NormaliseNetworkParam(const float minLimit, const float maxLimit, float &networkParam) const;

    /**
     *  @brief  Whether an input vector is not equal to a default value, 
     *          and therefore has been correctly set
     *
     *  @param  vector the input vector
     *
     *  @return whether the input vector is not equal to a default value
     */
    bool IsVectorSet(const pandora::CartesianVector &vector) const;

    /**
     *  @brief  Add a vector of network parameter values to a torch vector 
     *
     *  @param  startIndex the index from which to begin the insert
     *  @param  paramVector the vector of network parameter values
     *  @param  modelInput the torch vector to which to add
     */
    int AddToInput(const int startIndex, const pandora::FloatVector &paramVector, LArDLHelper::TorchInput &modelInput) const;

    float m_bogusFloat;                   ///< a default float value
    float m_vertexRegionRadius;           ///< the radius in which to search for particle hits
    pandora::StringVector m_pfoListNames; ///< the input pfo list name vector
    float m_detectorMinX;                 ///< the minimum x detector boundary
    float m_detectorMaxX;                 ///< the maximum x detector boundary
    float m_detectorMinY;                 ///< the minimum y detector boundary
    float m_detectorMaxY;                 ///< the maximum y detector boundary
    float m_detectorMinZ;                 ///< the minimum z detector boundary
    float m_detectorMaxZ;                 ///< the maximum z detector boundary
    bool areBoundariesSet;                ///< whether the detector boundaries have been set
};

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_BASE_HIERARCHY_TOOL_H
