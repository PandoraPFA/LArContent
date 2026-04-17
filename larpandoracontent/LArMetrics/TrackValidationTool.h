/**
 *  @file   larpandoracontent/LArMetrics/TrackValidationTool.h
 *
 *  @brief  Header file for the track validation tool class.
 *
 *  $Log: $
 */
#ifndef TRACK_VALIDATION_TOOL_H
#define TRACK_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  TrackValidationTool class
 */
class TrackValidationTool : public BaseValidationTool
{

struct TrackTreeVars
{
    int m_run;    ///< run number
    int m_subrun; ///< subrun number
    int m_event;  ///< event number
    pandora::FloatVector m_recoEndpointX;   ///< x-coordinate of the track endpoint obtained from a sliding fit [cm]
    pandora::FloatVector m_recoEndpointY;   ///< y-coordinate of the track endpoint obtained from a sliding fit [cm] 
    pandora::FloatVector m_recoEndpointZ;   ///< z-coordinate of the track endpoint obtained from a sliding fit [cm]
    pandora::FloatVector m_recoEndpointAcc; ///< signed true-reco endpoint separation (+ve = downstream, -ve = upstream) [cm]
    pandora::FloatVector m_recoStartDirX;   ///< x-component of the initial direction obtained from a sliding fit
    pandora::FloatVector m_recoStartDirY;   ///< y-component of the initial direction obtained from a sliding fit
    pandora::FloatVector m_recoStartDirZ;   ///< z-component of the initial direction obtained from a sliding fit
    pandora::FloatVector m_startDirAcc;     ///< opening angle between true and reco initial directions [radians]
    pandora::FloatVector m_recoEndDirX;     ///< x-component of the final direction obtained from a sliding fit
    pandora::FloatVector m_recoEndDirY;     ///< y-component of the final direction obtained from a sliding fit
    pandora::FloatVector m_recoEndDirZ;     ///< z-component of the final direction obtained from a sliding fit
    pandora::FloatVector m_endDirAcc;       ///< opening angle between true and reco final directions [radians]
    pandora::IntVector m_isCorrectOrientation; ///< whether the reco track has the correct orientation (i.e. endpoint is downstream of the start point)
    pandora::IntVector m_isLeaving;         ///< whether the track is leaving the TPC
    pandora::IntVector m_nEndpointMCHits;   ///< number of 2D MCParticle hits in the Xcm preceding the MCParticle endpoint
    pandora::IntVector m_nEndpointMCHitsU;  ///< number of U-view MCParticle hits in the Xcm preceding the MCParticle endpoint
    pandora::IntVector m_nEndpointMCHitsV;  ///< number of V-view MCParticle hits in the Xcm preceding the MCParticle endpoint
    pandora::IntVector m_nEndpointMCHitsW;  ///< number of W-view MCParticle hits in the Xcm preceding the MCParticle endpoint
    pandora::IntVector m_nEndpointPfoHits;  ///< number of 2D Pfo hits in the Xcm preceding the MCParticle endpoint
    pandora::IntVector m_nEndpointPfoHitsU; ///< number of U-view Pfo hits in the Xcm preceding the MCParticle endpoint
    pandora::IntVector m_nEndpointPfoHitsV; ///< number of V-view Pfo hits in the Xcm preceding the MCParticle endpoint
    pandora::IntVector m_nEndpointPfoHitsW; ///< number of W-view Pfo hits in the Xcm preceding the MCParticle endpoint
    pandora::FloatVector m_endpointCompleteness;  ///< completeness of the pfo in the region preceding the MCParticle endpoint
    pandora::FloatVector m_endpointCompletenessU; ///< completeness of the pfo in the region preceding the MCParticle endpoint in the U-view
    pandora::FloatVector m_endpointCompletenessV; ///< completeness of the pfo in the region preceding the MCParticle endpoint in the V-view
    pandora::FloatVector m_endpointCompletenessW; ///< completeness of the pfo in the region preceding the MCParticle endpoint in the W-view
    pandora::FloatVector m_endpointPurity;        ///< purity of the pfo in the region preceding the MCParticle endpoint
    pandora::FloatVector m_endpointPurityU;       ///< purity of the pfo in the region preceding the MCParticle endpoint in the U-view
    pandora::FloatVector m_endpointPurityV;       ///< purity of the pfo in the region preceding the MCParticle endpoint in the V-view
    pandora::FloatVector m_endpointPurityW;       ///< purity of the pfo in the region preceding the MCParticle endpoint in the W-view
    pandora::IntVector m_hasMichel;       ///< whether the MCParticle has a true michel (grand)child at its endpoint
    pandora::IntVector m_michelFromMuon;  ///< whether the parent is a muon or pion
    pandora::IntVector m_hasTargetMichel; ///< whether the MCParticle has a visible michel (grand)child at its endpoint
    pandora::IntVector m_hasRecoMichel;   ///< whether the michel has at least one matched pfo
    pandora::IntVector m_michelIndex;     ///< the index of the michel in the particle vector
    pandora::IntVector m_michelIsChild;   ///< whether the parent-child link was reconstructed correctly
    pandora::IntVector m_michelIsShower;  ///< whether the michel pfo is shower-like
};

public:
    /**
     *  @brief  Default constructor
     */
    TrackValidationTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

private:
    /**
     *  @brief  Fill the variables that describe the true vertex and endpoint
     *
     *  @param  pMCParticle the target MCParticle
     *  @param  trackTreeVars the track tree variables
    */
    void GetTrueVertexAndEndpointVars(const pandora::MCParticle *const pMCParticle, TrackTreeVars &trackTreeVars);

    /**
     *  @brief  Fit the pfo following the PandoraTrack implementation
     *
     *  @param  pPfo the pfo to fit
     *  @param  vertex the fitted vertex position
     *  @param  vertexDir the fitted vertex direction
     *  @param  endpoint the fitted endpoint position
     *  @param  endpointDir the fitted endpoint direction
     *
     *  @return whether the fit was successful
    */
    bool FitTrack(const pandora::Pfo *const pPfo, pandora::CartesianVector &vertex, pandora::CartesianVector &vertexDir, 
        pandora::CartesianVector &endpoint, pandora::CartesianVector &endpointDir);

    /**
     *  @brief  Fill the variables that describe the reconstructed vertex and endpoint
     *
     *  @param  pMCParticle the target MCParticle
     *  @param  recoVertex the fitted vertex position
     *  @param  recoVertexDir the fitted vertex direction
     *  @param  recoEndpoint the fitted endpoint position
     *  @param  recoEndpointDir the fitted endpoint direction
     *  @param  trackTreeVars the track tree variables
     */
    void GetRecoVertexAndEndpointVars(const pandora::MCParticle *const pMCParticle, const pandora::CartesianVector &recoVertex, 
        const pandora::CartesianVector &recoVertexDir, const pandora::CartesianVector &recoEndpoint, 
        const pandora::CartesianVector &recoEndpointDir, TrackTreeVars &trackTreeVars);

    /**
     *  @brief  Fill the variables that describe the completeness/purity of the track end region
     *
     *  @param  mcMatchesVec the mapping of true nodes (MCParticles) to reco node (pfo) matches
     *  @param  pMCParticle the target MCParticle
     *  @param  pPfo the best-matched pfo
     *  @param  trackTreeVars the track tree variables
     */
    void GetTrueEndRegionVars(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticle *const pMCParticle,
        const pandora::Pfo *const pPfo, TrackTreeVars &trackTreeVars);

    /**
     *  @brief  Fill the variables that describe (grand)child michel MCParticles and their reconstruction
     *
     *  @param  targetMC the vector of target MCParticles
     *  @param  bestRecoMatch the vector of best-matched pfo's
     *  @param  trackTreeVars the track tree variables
     */   
    void MichelValidation(const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch, TrackTreeVars &trackTreeVars);

    /**
     *  @brief  Fill relevant tree variables with default values if no best-matched pfo, or a failed track fit
     *
     *  @param  trackTreeVars the track tree variables
     */
    void FillForFailedPfo(TrackTreeVars &trackTreeVars);

    /**
     *  @brief  Fill the track tree
     *
     *  @param  trackTreeVars the track tree variables
     */
    void FillTree(TrackTreeVars &trackTreeVars);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_edgeBuffer;     ///< the distance from the TPC edge within which a track endpoint is considered to be leaving the TPC [cm]
    int m_slidingFitWindow; ///< the sliding fit window used in the track fit
    float m_endRegion;      ///< the distance from the MCParticle endpoint within which hits are considered for the end region metrics [cm]
};

} // namespace lar_content

#endif // #ifndef TRACK_VALIDATION_TOOL_H
