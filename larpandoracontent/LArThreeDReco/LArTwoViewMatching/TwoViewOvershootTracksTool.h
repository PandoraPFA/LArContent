/**
 *  @file
 * larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewOvershootTracksTool.h
 *
 *  @brief  Header file for the two view overshoot tracks tool class.
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
     *  @param  elementA the matrix element A
     *  @param  elementB the matrix element B
     */
        Particle(const MatrixType::Element &elementA, const MatrixType::Element &elementB);

        const pandora::Cluster *m_pCommonCluster;       ///< Address of the common cluster
        const pandora::Cluster *m_pClusterA;            ///< Address of cluster in element A
        const pandora::Cluster *m_pClusterB;            ///< Address of cluster in element B
        pandora::CartesianVector m_splitPositionAB;     ///< The candidate split position in view with A and
                                                        ///< B
        pandora::CartesianVector m_splitPositionCommon; ///< The candidate split position for the common
                                                        ///< cluster
    };

    void GetIteratorListModifications(
        TwoViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
   *  @brief  Whether a pair of vertices pass longitudinal projection cuts
   *
   *  @param  vertexA vertex from cluster in matrix element a
   *  @param  vertexB vertex from cluster in matrix element b
   */
    bool PassesVertexCuts(const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB) const;

    /**
   *  @brief  Set split position for a provided particle
   *
   *  @param  vertexA vertex for matrix element a
   *  @param  vertexB vertex for matrix element b
   *  @param  particle the particle
   */
    void SetSplitPosition(const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB, Particle &particle) const;

    /**
   *  @brief  Whether the provided particle is consistent with being a kink,
   * when examined in three dimensions at the split position
   *
   *  @param  pAlgorithm the calling algorithm
   *  @param  particle the particle
   *  @param  isALowestInX whether cluster associated with matrix element a
   * extends to lowest x positions
   *
   *  @return boolean
   */
    bool IsThreeDKink(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const Particle &particle, const bool isALowestInX) const;

    bool m_splitMode;                 ///< Whether to run in cluster splitting mode, as opposed to
                                      ///< cluster merging mode
    float m_cosThetaCutForKinkSearch; ///< The cos theta cut used for the kink
                                      ///< search in three dimensions
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_OVERSHOOT_TRACKS_TOOL_H
