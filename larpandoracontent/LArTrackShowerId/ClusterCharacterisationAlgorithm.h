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
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    ClusterCharacterisationAlgorithm();

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

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterCharacterisationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
