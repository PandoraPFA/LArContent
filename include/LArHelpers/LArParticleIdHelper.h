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

#include "LArObjects/LArTwoDSlidingFitResult.h"

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
     *  @brief  Photon identification for use with lar tpcs
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArPhotonId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Electron identification for use with lar tpcs
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArElectronId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Muon identification for use with lar tpcs
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

private:
    /**
     *  @brief  Get the muon track width estimator for a provided sliding fit result
     * 
     *  @param  twoDSlidingFitResult the sliding fit result
     * 
     *  @return the muon track width estimator
     */
    static float GetMuonTrackWidth(const TwoDSlidingFitResult &twoDSlidingFitResult);

    static unsigned int     m_muonIdLayerFitHalfWindow;         ///< Layer fit half window, used for calculating sliding muon track width
    static float            m_muonIdMinLayerOccupancy;          ///< Min layer occupancy for for muon identification
    static float            m_muonIdMaxTrackWidth;              ///< Max muon track width estimator for muon identification
    static float            m_muonIdTrackResidualQuantile;      ///< Track residual quantile, used for calculating muon track width
};

} // namespace lar

#endif // #ifndef LAR_PARTICLE_ID_HELPER_H
