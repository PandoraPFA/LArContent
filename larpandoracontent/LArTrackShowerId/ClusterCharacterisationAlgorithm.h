/**
 *  @file   larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the cluster characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClusterCharacterisationAlgorithm class
 */
class ClusterCharacterisationAlgorithm : public ShowerGrowingAlgorithm
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

    /**
     *  @brief  Destructor
     */
    ~ClusterCharacterisationAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Whether cluster is identified as a clear track
     *
     *  @param  pCluster address of the relevant cluster
     *  @param  pClusterList the address of the list of clusters in the same view
     *
     *  @return boolean
     */
    bool IsClearTrack(const pandora::Cluster *const pCluster, const pandora::ClusterList *const pClusterList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_inputClusterListNames;        ///< The names of the input cluster lists

    unsigned int            m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    unsigned int            m_minHitsInCluster;             ///< The minimum number of hits in a clear track candidate
    float                   m_maxLayerGapFraction;          ///< The maximum (sliding fit) layer gap fraction for a clear track candidate
    float                   m_maxWidthPerUnitLength;        ///< The maximum width per unit length for a clear track candidate
    float                   m_maxShowerLength;              ///< The maximum length for a cluster to be considered a shower candidate
    bool                    m_useDetectorGaps;              ///< Whether to account for registered detector gaps in characterisation

    bool                    m_overwriteExistingId;          ///< Whether to consider any clusters that already have an assigned particle id
    bool                    m_useUnavailableClusters;       ///< Whether to consider clusters that are already constituents of a pfo

    bool                    m_writeToTree;                  ///< Whether to write monitoring details to tree
    std::string             m_treeName;                     ///< Name of output tree
    std::string             m_fileName;                     ///< Name of output file
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterCharacterisationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
