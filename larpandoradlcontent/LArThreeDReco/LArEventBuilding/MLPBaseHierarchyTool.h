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

    void SetDetectorBoundaries();

    bool IsInFV(const pandora::CartesianVector &position) const;

    float GetNSpacepoints(const HierarchyPfo &hierarchyPfo);

    std::pair<float, float> GetParticleInfoAboutPfoPosition(const pandora::Algorithm *const pAlgorithm, 
        const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &pointOfInterest) const;

    void NormaliseNetworkParam(const float minLimit, const float maxLimit, float &networkParam) const;

    bool IsVectorSet(const pandora::CartesianVector &vector);

    int AddToInput(const int startIndex, const pandora::FloatVector &paramVector, LArDLHelper::TorchInput &modelInput);

    float m_detectorMinX;
    float m_detectorMaxX;
    float m_detectorMinY;
    float m_detectorMaxY;
    float m_detectorMinZ;
    float m_detectorMaxZ;
    bool areBoundariesSet;
    float m_vertexRegionRadius;
    pandora::StringVector m_pfoListNames;
};

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_BASE_HIERARCHY_TOOL_H
