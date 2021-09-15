/**
 *  @file   larpandoracontent/LArTrackShowerId/DlClusterCharacterisationBaseAlgorithm.h
 *
 *  @brief  Header file for the deep learning cluster characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_DL_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationBaseAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DlClusterCharacterisationBaseAlgorithm class
 */
class DlClusterCharacterisationAlgorithm : public lar_content::ClusterCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlClusterCharacterisationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~DlClusterCharacterisationAlgorithm();

protected:
    /**
     *  @brief  Whether cluster is identified as a clear track
     *
     *  @param  pCluster address of the relevant cluster
     *
     *  @return boolean
     */
    bool IsClearTrack(const pandora::Cluster *const pCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_CLUSTER_CHARACTERISATION_ALGORITHM_H
