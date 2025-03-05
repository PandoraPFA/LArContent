/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLNeutrinoHierarchyAlgorithm.h
 *
 *  @brief  Header file for the DL neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_NEUTRINO_HIERARCHY_ALGORITHM_H
#define LAR_DL_NEUTRINO_HIERARCHY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoradlcontent/LArCheating/DLCheatHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLLaterTierHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLPrimaryHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DLNeutrinoHierarchyAlgorithm class
 */
class DLNeutrinoHierarchyAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DLNeutrinoHierarchyAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~DLNeutrinoHierarchyAlgorithm();

private:
    typedef std::map<const pandora::ParticleFlowObject *, const pandora::MCParticle *> PfoToMCParticleMap;
    typedef std::map<const pandora::ParticleFlowObject *, std::pair<const pandora::ParticleFlowObject *, int>> ChildToParentPfoMap;
    typedef std::vector<pandora::PfoVector> Hierarchy;

    pandora::StatusCode Run();

    /**
     *  @brief  Return the neutrino pfo
     *
     *  @param  the pointer to the neutrino pfo to fill
     *
     *  @return whether the neutrino pfo can be found
     */
    bool GetNeutrinoPfo(const pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
     *  @brief  Fill the track and shower-like HierarchyPfoVector
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     */
    void FillTrackShowerVectors(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoVector &trackPfos, HierarchyPfoVector &showerPfos) const;

    /**
     *  @brief  Identify the upstream (closest to the nu vertex) and downstream (furthest 
     *          from the nu vertex) endpoints and directions of an input pfo  
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  pPfo a pointer to the pfo
     *  @param  slidingFitResult the sliding fit result of the input pfo 
     *  @param  upstreamPoint the upstream endpoint
     *  @param  downstreamPoint the downstream endpoint
     *
     *  @return whether this was successful
     */
    bool GetExtremalVerticesAndDirections(const pandora::ParticleFlowObject *const pNeutrinoPfo, const pandora::ParticleFlowObject *const pPfo,
        const ThreeDSlidingFitResult &slidingFitResult, ExtremalPoint &upstreamPoint, ExtremalPoint &downstreamPoint) const;

    /**
     *  @brief  Obtain the direction of a shower at a given endpoint from
     *          the angular decomposition of the shower's 'close' 3D hits
     *
     *  @param  pPfo a pointer to the pfo
     *  @param  vertex the position at which to find the direction
     *  @param  direction the found direction
     *
     *  @return whether a shower direction could be found
     */
    bool GetShowerDirection(
        const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &vertex, pandora::CartesianVector &direction) const;

    /**
     *  @brief  Set the primary network scores of the HierarchyPfos
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     */
    void SetPrimaryScores(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoVector &trackPfos, HierarchyPfoVector &showerPfos) const;

    /**
     *  @brief  Call the primary network tool to obtain the primary score of a HierarchyPfo
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  hierarchyPfo the HierarchyPfo object of the input pfo 
     *
     *  @return the primary network score
     */
    float GetPrimaryScore(
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoVector &trackPfos, const HierarchyPfo &hierarchyPfo) const;

    /**
     *  @brief  Add particles to the hierarchy if they satisfy a criteria
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  buildPrimaryTier whether to build the primary tier (or the later tier)
     *  @param  usePrimaryScore whether to assess the primary score (or the later tier score)
     *  @param  trackThreshold the threshold network score for tracks
     *  @param  showerThreshold the threshold network score for showers
     *  @param  isLowerThreshold whether the threshold is a minimum or maximum requirement
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     *  @param  hierarchy the hierarchy object to modify
     */
    void UpdateHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo, const bool buildPrimaryTier, const bool usePrimaryScore,
        const float trackThreshold, const float showerThreshold, const bool isLowerThreshold, HierarchyPfoVector &trackPfos,
        HierarchyPfoVector &showerPfos, Hierarchy &hierarchy) const;

    /**
     *  @brief  Set the predicted parent and later tier score of the
     *          HierarchyPfos to be considered for the later tiers
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     */
    void SetLaterTierScores(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoVector &trackPfos, HierarchyPfoVector &showerPfos) const;

    /**
     *  @brief  Call the later tier network tool to obtain the later tier score
     *          of a parent-child link
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  parentPfo the HierarchyPfo object of the parent pfo
     *  @param  childPfo the HierarchyPfo object of the child pfo     
     *
     *  @return the parent-child network score
     */
    float GetLaterTierScore(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo) const;

    /**
     *  @brief  Register the parent-child relationships in Pandora
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     */
    void BuildPandoraHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoVector &trackPfos,
        const HierarchyPfoVector &showerPfos) const;

    //------------------------------------------------------------------------------------------------------------------------------------------
    // Training functions
    //------------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Whether the event should be trained on
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *
     *  @return whether to train on the links in the event
     */
    bool ShouldTrainOnEvent(const pandora::ParticleFlowObject *const pNeutrinoPfo) const;

    /**
     *  @brief  Obtain the [pfo -> unique ID] map
     *
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector     
     *  @param  particleIDMap the [pfo -> unique ID] map to populate
     */
    void GetParticleIDMap(const HierarchyPfoVector &trackPfos, const HierarchyPfoVector &showerPfos,
        std::map<const pandora::ParticleFlowObject *, int> &particleIDMap) const;

    /**
     *  @brief  Get the training cuts for the later tier parent-child link
     *          these describe how well the child points back to the parent
     *
     *  @param  parentHierarchyPfo the HierarchyPfo object of the parent
     *  @param  childHierarchyPfo the HierarchyPfo object of the child
     *  @param  trueParentOrientation whether the parent POI is in the upstream position
     *  @param  trueChildOrientation whether the child POI is in the upstream position
     *
     *  @return the training cut pair
     */
    std::pair<float, float> GetTrainingCuts(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo,
        const bool trueParentOrientation, const bool trueChildOrientation) const;

    /**
     *  @brief  Fill the output event tree (containing event-level info)
     *
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector          
     *  @param  nPrimaryTrackLinks the number of recorded primary track links
     *  @param  nPrimaryShowerLinks the number of recorded primary shower links
     *  @param  nLaterTierTrackTrackLinks the number of recorded track-track later tier links
     *  @param  nLaterTierTrackShowerLinks the number of recorded track-shower later tier links     
     */
    void FillEventTree(const HierarchyPfoVector &trackPfos, const HierarchyPfoVector &showerPfos, const int nPrimaryTrackLinks,
        const int nPrimaryShowerLinks, const int nLaterTierTrackTrackLinks, const int nLaterTierTrackShowerLinks) const;

    /**
     *  @brief  Fill the primary trees (containing primary link info)
     *
     *  @param  matchingMap the [pfo -> mcparticle] matching map
     *  @param  childToParentPfoMap the [child pfo -> true parent pfo] map
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     *  @param  particleIDMap the [pfo -> unique ID] map
     *  @param  nPrimaryTrackLinks the number of recorded primary track links
     *  @param  nPrimaryShowerLinks the number of recorded primary shower links
     */
    void FillPrimaryTrees(const PfoToMCParticleMap &matchingMap, const ChildToParentPfoMap &childToParentPfoMap,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoVector &trackPfos, const HierarchyPfoVector &showerPfos,
        const std::map<const pandora::ParticleFlowObject *, int> &particleIDMap, int &nPrimaryTrackLinks, int &nPrimaryShowerLinks) const;

    /**
     *  @brief  Enter an entry into the specified primary tree
     *
     *  @param  treeName the name of the tree to fill
     *  @param  isTrainingLink whether the input link is one to train on
     *  @param  isTrueLink whether the input link is a true parent-child link
     *  @param  isOrientationCorrect whether the assumed orientation of the particle is correct
     *  @param  trueVisibleGen the true visible generation of the particle
     *  @param  truePDG the PDG of the matched MCParticle
     *  @param  trueParentID the unique ID of the true parent
     *  @param  particleID the unique ID of the particle     
     *  @param  primaryNetworkParams the network parameter variable values of the input link
     */
    void FillPrimaryTree(const std::string &treeName, const bool isTrainingLink, const bool isTrueLink, const bool isOrientationCorrect,
        const int trueVisibleGen, const int truePDG, const int trueParentID, const int particleID,
        const DLPrimaryHierarchyTool::DLPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Fill the later tier trees (containing later tier link info)
     *
     *  @param  matchingMap the [pfo -> mcparticle] matching map
     *  @param  childToParentPfoMap the [child pfo -> true parent pfo] map
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     *  @param  particleIDMap the [pfo -> unique ID] map
     *  @param  nTrackLinks the number of recorded track-track later tier links
     *  @param  nShowerLinks the number of recorded track-shower later tier links
     */
    void FillLaterTierTrees(const PfoToMCParticleMap &matchingMap, const ChildToParentPfoMap &childToParentPfoMap,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoVector &trackPfos, const HierarchyPfoVector &showerPfos,
        const std::map<const pandora::ParticleFlowObject *, int> &particleIDMap, int &nTrackLinks, int &nShowerLinks) const;

    /**
     *  @brief  Enter an entry into the specified later tier tree
     *
     *  @param  treeName the name of the tree to fill
     *  @param  isTrainingLink whether the input link is one to train on
     *  @param  isTrueLink whether the input link is a true parent-child link
     *  @param  isOrientationCorrect whether the assumed orientation of the particles is correct
     *  @param  trueChildGen the true visible generation of the child particle
     *  @param  the training cut pair of the parent-child link
     *  @param  parentID the unique ID of the parent particle
     *  @param  childID the unique ID of the child particle     
     *  @param  networkParams the network parameter variable values of the input link
     */
    void FillLaterTierTree(const std::string &treeName, const bool isTrainingLink, const bool isTrueLink, const bool isOrientationCorrect,
        const int childTrueGen, const std::pair<float, float> &trainingCuts, const int parentID, const int childID,
        const DLLaterTierHierarchyTool::DLLaterTierNetworkParams &networkParams) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    int m_eventID;                                      ///< the event counter
    bool m_trainingMode;                                ///< whether to run in training mode
    std::string m_trainingFileName;                     ///< the name of the output training file
    std::string m_eventTreeName;                        ///< the name of the output event tree
    std::string m_primaryTrackTreeName;                 ///< the name of the output primary track link tree
    std::string m_primaryShowerTreeName;                ///< the name of the output primary shower link tree
    std::string m_laterTierTrackTrackTreeName;          ///< the name of the output track-track later tier link tree
    std::string m_laterTierTrackShowerTreeName;         ///< the name of the output track-shower later tier link tree
    std::string m_mcParticleListName;                   ///< the name of the MCParticle list
    float m_trainingVertexAccuracy;                     ///< the maximum true-reco nu vertex displacement allowed for training
    std::string m_neutrinoPfoListName;                  ///< the name of the neutrino pfo list
    pandora::StringVector m_pfoListNames;               ///< the name of the pfo lists
    unsigned int m_minClusterSize;                      ///< the minimum threshold of 3D hits of a considered pfo
    int m_slidingFitWindow;                             ///< the sliding fit window to use in pfo sliding linear fits
    float m_regionForDirFit;                            ///< the radius of the region used in shower direction fits
    int m_nAngularBins;                                 ///< the number of angle bins used by the shower direction fitter
    float m_primaryRegion;                              ///< the radius of the nu vertex region where particles are assumed to be primaries
    float m_primaryThresholdTrackPass1;                 ///< The threshold applied to tracks in pass 1 primary tier building
    float m_primaryThresholdShowerPass1;                ///< The threshold applied to showers in pass 1 primary tier building
    float m_laterTierThresholdTrackPass1;               ///< The threshold applied to track-track links in pass 1 later tier building
    float m_laterTierThresholdShowerPass1;              ///< The threshold applied to track-shower links in pass 1 later tier building
    float m_primaryThresholdTrackPass2;                 ///< The threshold applied to tracks in pass 2 primary tier building
    float m_primaryThresholdShowerPass2;                ///< The threshold applied to showers in pass 2 primary tier building
    float m_laterTierThresholdTrackPass2;               ///< The threshold applied to track-track links in pass 2 later tier building
    float m_laterTierThresholdShowerPass2;              ///< The threshold applied to track-shower links in pass 2 later tier building
    DLPrimaryHierarchyTool *m_primaryHierarchyTool;     ///< The tool used to build the primary tiers
    DLLaterTierHierarchyTool *m_laterTierHierarchyTool; ///< The tool used to build the later tiers
    DLCheatHierarchyTool *m_cheatHierarchyTool;         ///< The tool used to obtain the true hierarchy
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_NEUTRINO_HIERARCHY_ALGORITHM_H
