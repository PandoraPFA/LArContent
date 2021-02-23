/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewOvershootTracksTool.h
 *
 *  @brief  Header file for the overshoot tracks tool class.
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_OVERSHOOT_TRACKS_TOOL_H
#define TWO_VIEW_OVERSHOOT_TRACKS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewThreeDKinkBaseTool.h"

namespace lar_content
{

/**
 *  @brief  TwoViewOvershootTracksTool class
 */
class TwoViewOvershootTracksTool : public TwoViewThreeDKinkBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewOvershootTracksTool();

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
        Particle(const MatrixType::Element &elementA, const MatrixType::Element &elementB);

        const pandora::Cluster     *m_pCommonCluster;       ///< Address of the common cluster
        const pandora::Cluster     *m_pClusterA;           ///< Address of cluster in element A, view 1
        //const pandora::Cluster     *m_pClusterA2;           ///< Address of cluster in element A, view 2
        const pandora::Cluster     *m_pClusterB;           ///< Address of cluster in element B, view 1
        //const pandora::Cluster     *m_pClusterB2;           ///< Address of cluster in element B, view 2
        pandora::CartesianVector    m_splitPositionAB;       ///< The candidate split position in view 1
        pandora::CartesianVector    m_splitPositionCommon;        ///< The candidate split position for the common cluster
        //pandora::CartesianVector    m_splitPosition2;       ///< The candidate split position in view 2
    };

    void GetIteratorListModifications(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const;
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
    void SetSplitPosition(const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB, Particle &particle) const;

    /**
     *  @brief  Whether the provided particle is consistent with being a kink, when examined in three dimensions at the split position
     *
     *  @param  pAlgorithm the calling algorithm
     *  @param  particle the particle
     *  @param  isA1LowestInX whether cluster associated with tensor element a extends to lowest x positions in view 1
     *  @param  isA2LowestInX whether cluster associated with tensor element a extends to lowest x positions in view 2
     *
     *  @return boolean
     */
    bool IsThreeDKink(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const Particle &particle, const bool isALowestInX) const;

    bool            m_splitMode;                        ///< Whether to run in cluster splitting mode, as opposed to cluster merging mode
//    float           m_maxVertexXSeparation;             ///< The max separation between accompanying clusters vertex x positions to make split
    float           m_cosThetaCutForKinkSearch;         ///< The cos theta cut used for the kink search in three dimensions
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_OVERSHOOT_TRACKS_TOOL_H
