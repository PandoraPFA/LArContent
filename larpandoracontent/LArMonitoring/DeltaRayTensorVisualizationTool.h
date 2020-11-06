/**
 *  @file   larpandoracontent/LArMonitoring/DeltaRayTensorVisualizationTool.h
 *
 *  @brief  Header file for the transverse tensor visualization tool class.
 *
 *  $Log: $
 */
#ifndef DELTA_RAY_TENSOR_VISUALIZATION_TOOL_H
#define DELTA_RAY_TENSOR_VISUALIZATION_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayTensorVisualizationTool class
 */
class DeltaRayTensorVisualizationTool : public DeltaRayTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayTensorVisualizationTool();

    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);

    typedef std::map<const pandora::MCParticle*, unsigned int> MCParticleToIDMap;
    typedef std::map<unsigned int, pandora::CaloHitList> IDToHitMap;

    void FillMCParticleIDMap(const pandora::Cluster *const pCluster, MCParticleToIDMap &mcParticleToIDMap, IDToHitMap &idToHitMap);
    void PrintClusterHitOwnershipMap(IDToHitMap &idToHitMap);
    const pandora::MCParticle *GetLeadingParticle(const pandora::MCParticle *const pMCParticle);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_minClusterConnections;        ///< The minimum number of cluster connections for display
    bool            m_ignoreUnavailableClusters;    ///< Whether to ignore (skip-over) unavailable clusters in the tensor
    bool            m_showEachIndividualElement;    ///< Whether to draw each individual tensor element
    bool            m_showContext;                  ///< Whether to show input cluster lists to add context to tensor elements
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_TENSOR_VISUALIZATION_TOOL_H
