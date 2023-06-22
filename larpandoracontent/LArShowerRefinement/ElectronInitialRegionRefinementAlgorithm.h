/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronInitialRegionRefinementAlgorithm.h
 *
 *  @brief  Header file for the electron initial region refinement algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ELECTRON_INITIAL_REGION_REFINEMENT_ALGORITHM_H
#define LAR_ELECTRON_INITIAL_REGION_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArShowerRefinement/ConnectionPathwayFeatureTool.h"
#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h"

namespace lar_content
{

/**
 *  @brief  ElectronInitialRegionRefinementAlgorithm class
 */
class ElectronInitialRegionRefinementAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ElectronInitialRegionRefinementAlgorithm();

private:
    typedef std::map<const pandora::MCParticle*, pandora::CaloHitList> HitOwnershipMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Obtain a sorted vector of the reconstructed shower pfos
     *
     *  @param  showerPfoVector the sorted shower vector to be filled
     */
    void FillShowerPfoVector(pandora::PfoVector &showerPfoVector) const;

    /**
     *  @brief  Find and evaluate shower connection pathway, removing if necessary  
     *
     *  @param  pShowerPfo the input shower pfo
     */
    void RefineShower(const pandora::ParticleFlowObject *const pShowerPfo) const;

    /**
     *  @brief  Obtain the reconstructed neutrino vertex
     *
     *  @param  nuVertex3D the output neutrino vertex
     *
     *  @return whether a reconstructed neutrino vertex could be found
     */
    pandora::StatusCode GetNeutrinoVertex(pandora::CartesianVector &nuVertex3D) const;

    /**
     *  @brief  Build the 2D ProtoShower objects for a given view
     *
     *  @param  pShowerPfo the input shower pfo
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  hitType the 2D view
     *  @param  protoShowerVector the output vector of ProtoShower objects
     */
    void BuildViewProtoShowers(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, 
        pandora::HitType hitType, ProtoShowerVector &protoShowerVector) const;

    /**
     *  @brief  Obtain the event hit list of a given view
     *
     *  @param  hitType the 2D view 
     *  @param  pCaloHitList the output 2D hit list
     *
     *  @return whether a valid 2D hit list could be found
     */
    pandora::StatusCode GetHitListOfType(const pandora::HitType hitType, const pandora::CaloHitList *&pCaloHitList) const;

    /**
     *  @brief  Fit the shower to obtain a 2D shower vertex
     *
     *  @param  pShowerPfo the input shower pfo
     *  @param  hitType the 2D view 
     *  @param  nuVertex3D the 3D neutrino vertex 
     *
     *  @return the 2D shower vertex position
     */
    pandora::CartesianVector GetShowerVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, 
        const pandora::CartesianVector &nuVertex3D) const;

    /**
     *  @brief  Move the shower vertex closer to the connection pathway
     *
     *  @param  pShowerPfo the input shower pfo
     *  @param  hitType the 2D view 
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  peakDirection the initial direction of the connection pathway
     *  @param  showerVertexPosition the (perhaps) refined 2D shower vertex position 
     */
    void RefineShowerVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, 
        const pandora::CartesianVector &nuVertex3D, const pandora::CartesianVector &peakDirection, pandora::CartesianVector &showerVertexPosition) const;

    /**
     *  @brief  To determine whether the shower vertex lies on the connection pathway
     *
     *  @param  showerVertexPosition the input 2D shower vertex position
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  peakDirection the initial direction of the connection pathway
     *
     *  @return whether the shower vertex lies on the connection pathway
     */
    bool IsShowerConnected(const pandora::CartesianVector &showerVertexPosition, const pandora::CartesianVector &nuVertex2D, 
        const pandora::CartesianVector &peakDirection) const;

    /**
     *  @brief  To determine if the hits downstream of the shower vertex lie within the shower
     *
     *  @param  pShowerPfo the input shower pfo 
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  hitType the 2D view
     *  @param  showerVertex the 2D shower vertex
     *  @param  showerSpineHitList the hits of the found shower spine
     *
     *  @return whether the hits downstream of the shower vertex lie within the shower
     */
    bool IsSpineCoincident(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, 
        const pandora::HitType hitType, const pandora::CartesianVector &showerVertex, const pandora::CaloHitList &showerSpineHitList) const;

    /**
     *  @brief  Build the connection pathways of all other particles in the event
     *
     *  @param  pShowerPfo the input shower pfo 
     *  @param  protectedHits the list of protected hits which will not be considered
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  hitType the 2D view
     *  @param  viewPathways the output vector of found connection pathways
     */
    void BuildViewPathways(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CaloHitList &protectedHits, const pandora::CartesianVector &nuVertex3D,
        pandora::HitType hitType, ConnectionPathwayVector &viewPathways) const;

    /**
     *  @brief  Determine the continuous and unambiguous hits to add to an electron-like shower pfo
     *
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  hitType the 2D view
     *  @param  viewPathways the vector of other event connection pathways
     *  @param  protoShower the ProtoShower object associated with the input shower
     */
    void RefineHitsToAdd(const pandora::CartesianVector &nuVertex3D, const pandora::HitType hitType, const ConnectionPathwayVector &viewPathways, 
        ProtoShower &protoShower) const;

    /**
     *  @brief  Find the continuous path of hits that lies closest to a given point
     *
     *  @param  refinedHitList the input hit list
     *  @param  nuVertex2D the given point (the 2D neutrino vertex)
     *
     *  @return a continuous pathway of hits
     */
    pandora::CaloHitList FindContinuousPath(const pandora::CaloHitList &refinedHitList, const pandora::CartesianVector &nuVertex2D) const;

    /**
     *  @brief  Add the shower characterisation information to the pfo metadata
     *
     *  @param  pShowerPfo the input shower pfo 
     *  @param  featureMap the map of [characterisation variable -> value]
     */
    void SetMetadata(const pandora::ParticleFlowObject *const pShowerPfo, const LArMvaHelper::MvaFeatureMap &featureMap) const;

    /**
     *  @brief  Determine the one-to-one mapping of leading MCParticle electrons and the hits which contain their energy depositions
     *
     *  @param  electronHitMap the output [MCParticle -> CaloHitList] map
     */
    void FillElectronHitMap(HitOwnershipMap &electronHitMap) const;

    /**
     *  @brief  To determine whether a pfo is a true leading electron via its completeness and purity
     *
     *  @param  pShowerPfo the input shower pfo 
     *  @param  electronHitMap the mapping of [MCParticle leading electrons -> their associated hits]
     */
    bool IsElectron(const pandora::ParticleFlowObject *const pShowerPfo, const HitOwnershipMap &electronHitMap) const;

    PeakDirectionFinderTool* m_pShowerPeakDirectionFinderTool;      ///< The shower initial pathway direction finder tool
    PeakDirectionFinderTool* m_pEventPeakDirectionFinderTool;       ///< The other (not incl. shower) initial pathway direction finder tool
    ShowerSpineFinderTool* m_pShowerSpineFinderTool;                ///< The shower spine finder tool for the shower
    ShowerSpineFinderTool* m_pEventPathwayFinderTool;               ///< The shower spine finder tool for all other event particles
    ShowerStartFinderTool* m_pShowerStartFinderTool;                ///< The shower start finder tool
    ProtoShowerMatchingTool* m_pProtoShowerMatchingTool;            ///< The 2D -> 3D ProtoShower matching tool
    std::string m_showerPfoListName;                                ///< The shower pfo list name
    std::string m_neutrinoVertexListName;                           ///< The neutrino vertex list name
    std::string m_caloHitListNameU;                                 ///< The U calo hit list name
    std::string m_caloHitListNameV;                                 ///< The V calo hit list name
    std::string m_caloHitListNameW;                                 ///< The W calo hit list name
    unsigned int m_minShowerHits3D;                                 ///< The min. number of hits of a significant shower
    unsigned int m_showerSlidingFitWindow;                          ///< The sliding fit window for shower fits
    float m_maxCoincidenceTransverseSeparation;                     ///< The max. transverse distance from the pathway direction of a coincident shower vertex 
    float m_minSpinePurity;                                         ///< The min. purity of a coincident shower spine downstream of the shower vertex
    bool m_trainingMode;                                            ///< Whether to run the algorithm to train the BDT
    float m_unambiguousThreshold;                                   ///< The min. transverse distance of an unambiguous shower hit from another pathway direction
    float m_maxConnectionDistance;                                  ///< The max. distance between connected hits
    unsigned int m_minNConnectedHits;                               ///< The number of connected hits needed for a conntected pathway
    float m_minElectronCompleteness;                                ///< The min. completeness of an electron-like pfo
    float m_minElectronPurity;                                      ///< The min. purity of an electron-like pfo
    float m_maxSeparationFromHit;                                   ///< The max. separation between the projected 3D shower start and the closest 2D shower hit
    float m_maxProjectionSeparation;                                ///< The max. separation between the projected 3D shower start and the shower start of that view
    float m_maxXSeparation;                                         ///< The max. drift-coordinate separation between a 3D shower start and a matched 2D shower hit 
    ConnectionPathwayFeatureTool::FeatureToolMap m_featureToolMap;  ///< The feature tool map
    pandora::StringVector m_algorithmToolNames;                     ///< The algorithm tool names
};

} // namespace lar_content

#endif // #ifndef LAR_ELECTRON_INITIAL_REGION_REFINEMENT
