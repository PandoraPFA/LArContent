/**
 *  @file   larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the cut based cluster characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CUT_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_CUT_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CutClusterCharacterisationAlgorithm class
 */
class CutClusterCharacterisationAlgorithm : public ClusterCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CutClusterCharacterisationAlgorithm();

    /**
     *  @brief  Get the distance between the interaction vertex (if present in the current vertex list) and a provided cluster
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pCluster address of the cluster
     *
     *  @return the vertex distance
     */
    static float GetVertexDistance(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get a measure of the width of a cluster, using a sliding shower fit result
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pCluster address of the cluster
     *  @param  showerFitWindow the layer window used for the sliding shower fit
     *
     *  @return the shower fit width
     */
    static float GetShowerFitWidth(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, const unsigned int showerFitWindow);

private:
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_slidingFitWindow;       ///< The layer window for the sliding linear fits
    unsigned int m_slidingShowerFitWindow; ///< The layer window for the sliding shower fits
    unsigned int m_minCaloHitsCut;         ///< The minimum number of calo hits to qualify as a track
    float m_maxShowerLengthCut;            ///< The maximum cluster length to qualify as a shower
    float m_pathLengthRatioCut;            ///< The maximum ratio of path length to straight line length to qualify as a track
    float m_rTWidthRatioCut;        ///< The maximum ratio of transverse fit position width to straight line length to qualify as a track
    float m_vertexDistanceRatioCut; ///< The maximum ratio of vertex separation to straight line length to qualify as a track
    float m_showerWidthRatioCut;    ///< The maximum ratio of shower fit width to straight line length to qualify as a track
};

} // namespace lar_content

#endif // #ifndef LAR_CUT_CLUSTER_CHARACTERISATION_ALGORITHM_H
