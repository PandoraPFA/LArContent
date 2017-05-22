/**
 *  @file   larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the cluster characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ClusterCharacterisationAlgorithm class
 */
class ClusterCharacterisationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterCharacterisationAlgorithm();

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
    pandora::StatusCode Run();

    /**
     *  @brief  Whether cluster is identified as a clear track
     *
     *  @param  pCluster address of the relevant cluster
     *
     *  @return boolean
     */
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_inputClusterListNames;        ///< The names of the input cluster lists

    bool                    m_zeroMode;                     ///< Whether to zero all existing cluster particle id, overrides all other parameters

    bool                    m_overwriteExistingId;          ///< Whether to consider any clusters that already have an assigned particle id
    bool                    m_useUnavailableClusters;       ///< Whether to consider clusters that are already constituents of a pfo

    unsigned int            m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    unsigned int            m_slidingShowerFitWindow;       ///< The layer window for the sliding shower fits
    unsigned int            m_minCaloHitsCut;               ///< The minimum number of calo hits to qualify as a track
    float                   m_maxShowerLengthCut;           ///< The maximum cluster length to qualify as a shower
    float                   m_pathLengthRatioCut;           ///< The maximum ratio of path length to straight line length to qualify as a track
    float                   m_rTWidthRatioCut;              ///< The maximum ratio of transverse fit position width to straight line length to qualify as a track
    float                   m_vertexDistanceRatioCut;       ///< The maximum ratio of vertex separation to straight line length to qualify as a track
    float                   m_showerWidthRatioCut;          ///< The maximum ratio of shower fit width to straight line length to qualify as a track
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H

