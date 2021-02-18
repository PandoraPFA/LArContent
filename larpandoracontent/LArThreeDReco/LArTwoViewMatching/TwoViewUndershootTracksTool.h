/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewUndershootTracksTool.h
 *
 *  @brief  Header file for the undershoot tracks tool class.
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_UNDERSHOOT_TRACKS_TOOL_H
#define TWO_VIEW_UNDERSHOOT_TRACKS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewThreeDKinkBaseTool.h"

namespace lar_content
{

/**
 *  @brief  TwoViewUndershootTracksTool class
 */
class TwoViewUndershootTracksTool : public TwoViewThreeDKinkBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewUndershootTracksTool();

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

        const pandora::Cluster      *m_pClusterA;            ///< Address of non-shared cluster in element A
        const pandora::Cluster      *m_pClusterB;            ///< Address of non-shared cluster in element B
        const pandora::Cluster      *m_pCommonCluster;      ///< Address of the common cluster in view 1
        //const pandora::Cluster      *m_pCommonCluster2;      ///< Address of the common cluster in view 2
    };

    void GetIteratorListModifications(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Whether the provided particle is consistent with being a kink, when examined in three dimensions at the provided split position
     *
     *  @param  pAlgorithm the calling algorithm
     *  @param  particle the particle
     *  @param  splitPosition the candidate split position
     *  @param  isALowestInX whether cluster associated with matrix element a extends to lowest x positions
     *
     *  @return boolean
     */
    bool IsThreeDKink(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const Particle &particle, const pandora::CartesianVector &splitPosition,
        const bool isALowestInX) const;

    bool            m_splitMode;                        ///< Whether to run in cluster splitting mode, as opposed to cluster merging mode
    float           m_maxTransverseImpactParameter;     ///< The maximum transverse impact parameter for connecting broken clusters
    float           m_minImpactParameterCosTheta;       ///< The minimum cos theta (angle between vertex directions) for connecting broken clusters
    float           m_cosThetaCutForKinkSearch;         ///< The cos theta cut used for the kink search in three dimensions
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_UNDERSHOOT_TRACKS_TOOL_H
