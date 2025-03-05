/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLBaseHierarchyTool.h
 *
 *  @brief  Header file for the DL base hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_DL_BASE_HIERARCHY_TOOL_H
#define LAR_DL_BASE_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

#include <torch/script.h>
#include <torch/torch.h>

namespace lar_dl_content
{

/**
 *   @brief  DLBaseHierarchyTool to calculate variables related to the initial shower region
 */
class DLBaseHierarchyTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLBaseHierarchyTool();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Set the detector boundaries 
     */
    void SetDetectorBoundaries();

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

    float m_vertexRegionRadiusSq;                               ///< the radius (squared) in which to search for particle hits
    pandora::StringVector m_pfoListNames;                       ///< the input pfo list name vector
    LArGeometryHelper::DetectorBoundaries m_detectorBoundaries; ///< the detector boundaries
    bool m_areBoundariesSet;                                    ///< whether the detector boundaries have been set
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_BASE_HIERARCHY_TOOL_H
