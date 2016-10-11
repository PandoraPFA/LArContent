/**
 *  @file   larpandoracontent/LArPlugins/LArParticleIdPlugins.h
 * 
 *  @brief  Header file for the lar particle id plugins class.
 * 
 *  $Log: $
 */
#ifndef LAR_PARTICLE_ID_PLUGINS_H
#define LAR_PARTICLE_ID_PLUGINS_H 1

#include "Plugins/ParticleIdPlugin.h"

namespace lar_content
{

class TwoDSlidingFitResult;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArParticleIdPlugins class
 */
class LArParticleIdPlugins
{
public:
    /**
     *   @brief  LArMuonId class
     */
    class LArMuonId : public pandora::ParticleIdPlugin
    {
    public:
        /**
         *  @brief  Default constructor
         */
        LArMuonId();

        bool IsMatch(const pandora::Cluster *const pCluster) const;
        bool IsMatch(const pandora::ParticleFlowObject *const pPfo) const;

    private:
        /**
         *  @brief  Get the muon track width estimator for a provided sliding fit result
         * 
         *  @param  twoDSlidingFitResult the sliding fit result
         * 
         *  @return the muon track width estimator
         */
        float GetMuonTrackWidth(const TwoDSlidingFitResult &twoDSlidingFitResult) const;

        pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

        unsigned int    m_layerFitHalfWindow;       ///< Layer fit half window, used for calculating sliding muon track width
        float           m_minLayerOccupancy;        ///< Min layer occupancy for for muon identification
        float           m_maxTrackWidth;            ///< Max muon track width estimator for muon identification
        float           m_trackResidualQuantile;    ///< Track residual quantile, used for calculating muon track width
        unsigned int    m_minClustersPassingId;     ///< Match pfo if at sufficient clusters in pfo pass the cluster particle id logic
    };
};

} // namespace lar_content

#endif // #ifndef LAR_PARTICLE_ID_PLUGINS_H
