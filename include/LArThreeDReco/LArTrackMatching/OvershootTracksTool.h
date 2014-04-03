/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/OvershootTracksTool.h
 * 
 *  @brief  Header file for the overshoot tracks tool class.
 * 
 *  $Log: $
 */
#ifndef OVERSHOOT_TRACKS_TOOL_H
#define OVERSHOOT_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDKinkBaseTool.h"

namespace lar
{

/**
 *  @brief  OvershootTracksTool class
 */
class OvershootTracksTool : public ThreeDKinkBaseTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    OvershootTracksTool();

private:
    /**
     *  @brief  Particle class
     */
    class Particle
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  elementA the tensor element A
         *  @param  elementB the tensor element B
         */
        Particle(const TensorType::Element &elementA, const TensorType::Element &elementB);

        pandora::Cluster           *m_pCommonCluster;       ///< Address of the common cluster
        pandora::Cluster           *m_pClusterA1;           ///< Address of cluster in element A, view 1
        pandora::Cluster           *m_pClusterA2;           ///< Address of cluster in element A, view 2
        pandora::Cluster           *m_pClusterB1;           ///< Address of cluster in element B, view 1
        pandora::Cluster           *m_pClusterB2;           ///< Address of cluster in element B, view 2
        pandora::CartesianVector    m_splitPosition;        ///< The candidate split position
    };

    void GetIteratorListModifications(ThreeDTransverseTracksAlgorithm *pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Whether a pair of vertices pass longitudinal projection cuts
     * 
     *  @param  vertexA vertex from cluster in tensor element a
     *  @param  vertexB vertex from cluster in tensor element b
     */
    bool PassesVertexCuts(const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB) const;

    /**
     *  @brief  Set split position for a provided particle
     * 
     *  @param  vertexA1 vertex for tensor element a in view 1
     *  @param  vertexA2 vertex for tensor element a in view 2
     *  @param  vertexB1 vertex for tensor element b in view 1
     *  @param  vertexB2 vertex for tensor element b in view 2
     *  @param  particle the particle
     */
    void SetSplitPosition(const LArPointingCluster::Vertex &vertexA1, const LArPointingCluster::Vertex &vertexA2,
        const LArPointingCluster::Vertex &vertexB1, const LArPointingCluster::Vertex &vertexB2, Particle &particle) const;

    bool    m_splitMode;                    ///< Whether to run in cluster splitting mode, as opposed to cluster merging mode
    float   m_maxVertexXSeparation;         ///< The max separation between accompanying clusters vertex x positions to make split
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *OvershootTracksTool::Factory::CreateAlgorithmTool() const
{
    return new OvershootTracksTool();
}

} // namespace lar

#endif // #ifndef OVERSHOOT_TRACKS_TOOL_H
