/**
 *  @file   larpandoracontent/LArTrackShowerId/ClusterCharacterisationBaseAlgorithm.h
 *
 *  @brief  Header file for the cluster characterisation base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_CHARACTERISATION_BASE_ALGORITHM_H
#define LAR_CLUSTER_CHARACTERISATION_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ClusterCharacterisationBaseAlgorithm class
 */
class ClusterCharacterisationBaseAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterCharacterisationBaseAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~ClusterCharacterisationBaseAlgorithm();

protected:
    /**
     *  @brief  Whether cluster is identified as a clear track
     *
     *  @param  pCluster address of the relevant cluster
     *
     *  @return boolean
     */
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const = 0;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputClusterListNames; ///< The names of the input cluster lists

    bool m_zeroMode; ///< Whether to zero all existing cluster particle id, overrides all other parameters

    bool m_overwriteExistingId;    ///< Whether to consider any clusters that already have an assigned particle id
    bool m_useUnavailableClusters; ///< Whether to consider clusters that are already constituents of a pfo
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_CHARACTERISATION_BASE_ALGORITHM_H
