/**
 *  @file   LArContent/include/LArHelpers/LArParticleIdHelper.h
 * 
 *  @brief  Header file for the lar particle id helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_PARTICLE_ID_HELPER_H
#define LAR_PARTICLE_ID_HELPER_H 1

#include "Pandora/PandoraInternal.h"
#include "Pandora/PandoraInputTypes.h"
#include "Pandora/StatusCodes.h"

#include "Xml/tinyxml.h"

namespace lar
{

/**
 *  @brief  LArParticleIdHelper class
 */
class LArParticleIdHelper
{
public:
    /**
     *  @brief  Whether a cluster is a candidate electromagnetic shower
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArEmShowerId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Photon identification for use with fine granularity particle flow detectors
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArPhotonId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Electron identification for use with fine granularity particle flow detectors
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArElectronId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Muon identification for use with fine granularity particle flow detectors
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArMuonId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Read the lar particle id settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar

#endif // #ifndef LAR_PARTICLE_ID_HELPER_H
