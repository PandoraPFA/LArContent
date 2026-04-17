/**
 *  @file   larpandoracontent/LArMetrics/ShowerValidationTool.h
 *
 *  @brief  Header file for the shower validation tool class.
 *
 *  $Log: $
 */
#ifndef SHOWER_VALIDATION_TOOL_H
#define SHOWER_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  ShowerValidationTool class
 */
class ShowerValidationTool : public BaseValidationTool
{
public:
/**
 *  @brief  ShowerTreeVars struct
 */
struct ShowerTreeVars
{
    int m_run;    ///< run number
    int m_subrun; ///< subrun number
    int m_event;  ///< event number
    pandora::FloatVector m_coreTrueLengthFromU; ///< 3D length that captures X% of the charge in the U-view
    pandora::FloatVector m_coreTrueLengthFromV; ///< 3D length that captures X% of the charge in the V-view
    pandora::FloatVector m_coreTrueLengthFromW; ///< 3D length that captures X% of the charge in the W-view
    pandora::FloatVector m_recoShrVtxX;    ///< x-coordinate of the pfo vertex given under a shower hypothesis
    pandora::FloatVector m_recoShrVtxY;    ///< y-coordinate of the pfo vertex given under a shower hypothesis
    pandora::FloatVector m_recoShrVtxZ;    ///< z-coordinate of the pfo vertex given under a shower hypothesis
    pandora::FloatVector m_recoShrDirX;    ///< x-component of the initial direction obtained by a PCA
    pandora::FloatVector m_recoShrDirY;    ///< y-component of the initial direction obtained by a PCA
    pandora::FloatVector m_recoShrDirZ;    ///< z-component of the initial direction obtained by a PCA
    pandora::FloatVector m_coreRecoLength; ///< 3D length that captures X% of the charge in the W-view
    pandora::FloatVector m_recoShrLength;  ///< 3D length obtained by a PCA fit to the shower
    pandora::FloatVector m_recoShrDirAcc;  //< opening angle between true and reco initial directions [radians]
    pandora::FloatVector m_moliereRadius;  ///< the transverse 3D distance (from the core shower axis) that captures 90% of the charge in the W-view
    pandora::IntVector m_nInitialMCHits;   ///< number of 2D MCParticle hits in the Xcm following the MCParticle vertex
    pandora::IntVector m_nInitialMCHitsU;  ///< number of U-view MCParticle hits in the Xcm following the MCParticle vertex
    pandora::IntVector m_nInitialMCHitsV;  ///< number of V-view MCParticle hits in the Xcm following the MCParticle vertex
    pandora::IntVector m_nInitialMCHitsW;  ///< number of W-view MCParticle hits in the Xcm following the MCParticle vertex
    pandora::IntVector m_nInitialPfoHits;  ///< number of 2D Pfo hits in the Xcm following the MCParticle vertex
    pandora::IntVector m_nInitialPfoHitsU; ///< number of U-view Pfo hits in the Xcm following the MCParticle vertex
    pandora::IntVector m_nInitialPfoHitsV; ///< number of V-view Pfo hits in the Xcm following the MCParticle vertex
    pandora::IntVector m_nInitialPfoHitsW; ///< number of W-view Pfo hits in the Xcm following the MCParticle vertex
    pandora::FloatVector m_initialCompleteness;  ///< completeness of the initial region
    pandora::FloatVector m_initialCompletenessU; ///< completeness of the initial region in the U-view
    pandora::FloatVector m_initialCompletenessV; ///< completeness of the initial region in the V-view
    pandora::FloatVector m_initialCompletenessW; ///< completeness of the initial region in the W-view
    pandora::FloatVector m_initialPurity;        ///< purity of the initial region
    pandora::FloatVector m_initialPurityU;       ///< purity of the initial region in the U-view
    pandora::FloatVector m_initialPurityV;       ///< purity of the initial region in the V-view
    pandora::FloatVector m_initialPurityW;       ///< purity of the initial region in the W-view
};

    /**
     *  @brief  Default constructor
     */
    ShowerValidationTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

private:
    /**
     *  @brief  Fill the variables that describe the true shower length in each view
     *
     *  @param  mcMatchesVec the mapping of true nodes (MCParticles) to reco node (pfo) matches
     *  @param  pMCParticle the target MCParticle
     *  @param  showerTreeVars the shower tree variables
     */
    void GetTrueLength(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticle *const pMCParticle,
        ShowerTreeVars &showerTreeVars);

    /**
     *  @brief  Fill the variables that describe the initial region of the shower
     *
     *  @param  mcMatchesVec the mapping of true nodes (MCParticles) to reco node (pfo) matches
     *  @param  pMCParticle the target MCParticle
     *  @param  pPfo the best-matched pfo
     *  @param  showerTreeVars the shower tree variables
     */
    void GetInitialRegionVars(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticle *const pMCParticle, 
        const pandora::Pfo *const pPfo, ShowerTreeVars &showerTreeVars);

    /**
     *  @brief  Fit the pfo following the PandoraShower implementation
     *
     *  @param  pPfo the best-matched pfo
     *  @param  showerVertex the fitted shower vertex
     *  @param  showerDirection the fitted shower direction
     *  @param  showerLength the fitted shower length
     *
     *  @return whether the fit was successful
     */
    bool FitShower(const pandora::Pfo *const pPfo, pandora::CartesianVector &showerVertex, 
        pandora::CartesianVector &showerDirection, float &showerLength);

    /**
     *  @brief  Fill the variables that describe the reconstructed shower vertex and direction
     *
     *  @param  recoShrVtx the fitted shower vertex
     *  @param  recoShrDir the fitted shower direction
     *  @param  recoShrLength the fitted shower length
     *  @param  pMC the target MCParticle (nullptr if the MCParticle has no direction)
     *  @param  showerTreeVars the shower tree variables
     */        
    void GetRecoVertexInfo(const pandora::CartesianVector &recoShrVtx, const pandora::CartesianVector &recoShrDir, const float recoShrLength,
        const pandora::MCParticle *const pMC, ShowerTreeVars &showerTreeVars);

    /**
     *  @brief  Fill the variables that describe the Moliere radius and core length of the shower
     *
     *  @param  pPfo the best-matched pfo
     *  @param  showerVertex the fitted shower vertex
     *  @param  showerDirection the fitted shower direction
     *  @param  showerTreeVars the shower tree variables
     */
    void GetMoliere(const pandora::Pfo *const pPfo, const pandora::CartesianVector &showerVertex, const pandora::CartesianVector &showerDirection,
        ShowerTreeVars &showerTreeVars);

    /**
     *  @brief  Fill relevant tree variables with default values for the case of a null MCParticle direction
     *
     *  @param  showerTreeVars the shower tree variables
     */
    void FillForNullMCDir(ShowerTreeVars &showerTreeVars);

    /**
     *  @brief  Fill relevant tree variables with default values if no best-matched pfo, or a failed shower fit
     *
     *  @param  showerTreeVars the shower tree variables
     */
    void FillForFailedPfo(ShowerTreeVars &showerTreeVars);

    /**
     *  @brief  Fill the shower tree
     *
     *  @param  showerTreeVars the shower tree variables
     */
    void FillTree(ShowerTreeVars &showerTreeVars);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_trueLengthEnergyFrac; ///< the fraction of the total energy in a given view that defines the true 3D shower length
    float m_initialRegion3D;      ///< the distance from the MCParticle vertex that defines the initial region
};

} // namespace lar_content

#endif // #ifndef SHOWER_VALIDATION_TOOL_H
