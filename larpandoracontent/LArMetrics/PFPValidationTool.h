/**
 *  @file   larpandoracontent/LArMetrics/PFPValidationTool.h
 *
 *  @brief  Header file for the pfp validation tool class.
 *
 *  $Log: $
 */
#ifndef PFP_VALIDATION_TOOL_H
#define PFP_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  PFPValidationTool class
 */
class PFPValidationTool : public BaseValidationTool
{
public:
/**
 *  @brief  PFPTreeVars struct
 */
struct PFPTreeVars
{
    int m_run;    ///< run number
    int m_subrun; ///< subrun number
    int m_event;  ///< event number
    pandora::IntVector m_truePDG;         ///< true PDG of the particle
    pandora::FloatVector m_trueEnergy;    ///< true particle energy 
    pandora::FloatVector m_trueVisEnergy; ///< true particle visible energy
    pandora::FloatVector m_trueThetaXZ;   ///< angle of the X-Z projected true direction to the Z-axis [radians]
    pandora::FloatVector m_trueThetaYZ;   ///< angle of the true direction to the Y-axis [radians]
    pandora::IntVector m_isTrack;    ///< whether the pfo is track-like
    pandora::IntVector m_isShower;   ///< whether the pfo is shower-like
    pandora::IntVector m_hasMatch;   ///< whether the MCParticle has at least one matched pfo
    pandora::IntVector m_nMCHitsU;   ///< number of U-view hits belonging to the MCParticle
    pandora::IntVector m_nMCHitsV;   ///< number of V-view hits belonging to the MCParticle
    pandora::IntVector m_nMCHitsW;   ///< number of W-view hits belonging to the MCParticle
    pandora::IntVector m_nMCHits2D;  ///< number of 2D hits belonging to the MCParticle
    pandora::IntVector m_nPfoHitsU;  ///< number of U-view hits belonging to the pfo
    pandora::IntVector m_nPfoHitsV;  ///< number of V-view hits belonging to the pfo
    pandora::IntVector m_nPfoHitsW;  ///< number of W-view hits belonging to the pfo
    pandora::IntVector m_nPfoHits2D; ///< number of 2D hits belonging to the pfo
    pandora::IntVector m_nPfoHits3D; ///< number of 3D hits belonging to the pfo
    pandora::FloatVector m_completeness;     ///< completeness of the pfo considering all 2D hits
    pandora::FloatVector m_completenessU;    ///< completeness of the pfo considering U-view hits
    pandora::FloatVector m_completenessV;    ///< completeness of the pfo considering V-view hits
    pandora::FloatVector m_completenessW;    ///< completeness of the pfo considering W-view hits
    pandora::FloatVector m_completenessADC;  ///< adc-weighted completeness of the pfo considering all 2D hits
    pandora::FloatVector m_completenessADCU; ///< adc-weighted completeness of the PFO considering U-view hits
    pandora::FloatVector m_completenessADCV; ///< adc-weighted completeness of the PFO considering V-view hits
    pandora::FloatVector m_completenessADCW; ///< adc-weighted completeness of the PFO considering W-view hits
    pandora::FloatVector m_purity;     ///< purity of the pfo considering all 2D hits
    pandora::FloatVector m_purityU;    ///< purity of the pfo considering U-view hits
    pandora::FloatVector m_purityV;    ///< purity of the pfo considering V-view hits
    pandora::FloatVector m_purityW;    ///< purity of the pfo considering W-view hits
    pandora::FloatVector m_purityADC;  ///< adc-weighted purity of the pfo considering all 2D hits
    pandora::FloatVector m_purityADCU; ///< adc-weighted purity of the pfo considering U-view hits
    pandora::FloatVector m_purityADCV; ///< adc-weighted purity of the pfo considering V-view hits
    pandora::FloatVector m_purityADCW; ///< adc-weighted purity of the pfo considering W-view hits
    pandora::FloatVector m_altCompleteness;      ///< completeness of the thief particle wrt the target particle
    pandora::FloatVector m_altPurity;            ///< purity of the thief particle wrt the target particle
    pandora::IntVector m_altPDG;                 ///< PDG of the MCParticle the thief best matches
    pandora::IntVector m_altIsUpstreamHierarchy; ///< is the MCParticle that the thief best matches, an ancestor of this MCParticle?
    pandora::IntVector m_altIsSameMC;            ///< is the MCParticle that the thief best matches, an ancestor of this MCParticle?
    pandora::FloatVector m_trueVertexX; ///< x-coordinate of the first trajectory point that lies within the detector [cm]
    pandora::FloatVector m_trueVertexY; ///< y-coordinate of the first trajectory point that lies within the detector [cm] 
    pandora::FloatVector m_trueVertexZ; ///< z-coordinate of the first trajectory point that lies within the detector [cm]
    pandora::FloatVector m_trueEndX;    ///< x-coordinate of the last trajectory point that lies within the detector [cm]
    pandora::FloatVector m_trueEndY;    ///< y-coordinate of the last trajectory point that lies within the detector [cm]
    pandora::FloatVector m_trueEndZ;    ///< z-coordinate of the last trajectory point that lies within the detector [cm]
    pandora::FloatVector m_trueDirX;    ///< x-component of the true direction at the MCParticle vertex
    pandora::FloatVector m_trueDirY;    ///< y-component of the true direction at the MCParticle vertex
    pandora::FloatVector m_trueDirZ;    ///< z-component of the true direction at the MCParticle vertex
    pandora::FloatVector m_trueEndDirX; ///< x-component of the true direction at the penultimate trajectory point
    pandora::FloatVector m_trueEndDirY; ///< y-component of the true direction at the penultimate trajectory point
    pandora::FloatVector m_trueEndDirZ; ///< z-component of the true direction at the penultimate trajectory point
    pandora::FloatVector m_trueLength;  ///< separation between the true vertex and true endpoint [cm]
    pandora::FloatVector m_trueDisplacement; ///< separation between the MCParticle and MCNu vertices [cm]
    pandora::FloatVector m_recoVertexX; ///< x-coordinate of the pfo's vertex [cm]
    pandora::FloatVector m_recoVertexY; ///< y-coordinate of the pfo's vertex [cm]
    pandora::FloatVector m_recoVertexZ; ///< z-coordinate of the pfo's vertex [cm]
    pandora::FloatVector m_vertexAcc;   ///< signed true-reco vertex separation (+ve = downstream, -ve = upstream) [cm]
    pandora::FloatVector m_recoLength;  ///< length from sqrt(LArPfoHelper::GetThreeDLengthSquared)
    pandora::FloatVector m_recoDisplacement; ///< separation between the reco particle-neutrino vertices [cm]
};

    /**
     *  @brief  Default constructor
     */
    PFPValidationTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Fill MCParticle variables
     *
     *  @param  pMCNu the neutrino MCParticle
     *  @param  pMCTarget the MCParticle
     *  @param  pfpTreeVars the pfo tree variables
     */
    void GetMCParticleInfo(const pandora::MCParticle *const pMCNu, const pandora::MCParticle *const pMCTarget, PFPTreeVars &pfpTreeVars);

    /**
     *  @brief  Fill the reco variables
     *
     *  @param  pMCTarget the MCParticle
     *  @param  pBestMatch the best-matched pfo
     *  @param  pfpTreeVars the pfo tree variables
     */
    void GetRecoParticleInfo(const pandora::MCParticle *const pMCTarget, const pandora::Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars);

    /**
     *  @brief  Fill the variables that describe MCParticle-pfo match
     *
     *  @param  mcMatchesVec the mapping of true nodes (MCParticles) to reco node (pfo) matches
     *  @param  pMCTarget the MCParticle
     *  @param  pBestMatch the best-matched pfo
     *  @param  pfpTreeVars the pfo tree variables
     */
    void GetMatchingInfo(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec,
        const pandora::MCParticle *const pMCTarget, const pandora::Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars);

    /**
     *  @brief  Fill the variables that describe the second-best MCParticle-pfo matches
     *
     *  @param  mcMatchesVec the mapping of true nodes (MCParticles) to reco node (pfo) matches
     *  @param  pMCTarget the MCParticle
     *  @param  pBestMatch the best-matched pfo
     *  @param  pfpTreeVars the pfo tree variables
     */
    void GetAltMatchInfo(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC,
        const pandora::MCParticle *const pMCTarget, const pandora::Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars);

    /**
     *  @brief  Calculate the matching metrics (completeness/purity) for a given match
     *
     *  @param  pMCNode the true node associated to the MCParticle
     *  @param  pRecoNode the reco node associated to the pfo
     *  @param  completeness the output completeness
     *  @param  purity the output purity
     */
    void GetAltMetrics(const LArHierarchyHelper::MCHierarchy::Node *const pMCNode, const LArHierarchyHelper::RecoHierarchy::Node *const pRecoNode,
        float &completeness, float &purity);

    /**
     *  @brief  Fill the pfp tree
     *
     *  @param  pfpTreeVars the pfo tree variables
     */
    void FillTree(PFPTreeVars &pfpTreeVars);

    const pandora::VertexList *m_pNuVertexList; ///< a pointer to the neutrino vertex
    std::string m_nuVertexListName; ///< the name of the neutrino vertex list
    float m_maxMichelSep;           ///< the maximum separation between a michel vertex and parent muon endpoint
};

} // namespace lar_content

#endif // #ifndef PFP_VALIDATION_TOOL_H
