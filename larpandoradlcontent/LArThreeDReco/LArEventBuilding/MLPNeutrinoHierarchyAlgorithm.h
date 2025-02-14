/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPNeutrinoHierarchyAlgorithm.h
 *
 *  @brief  Header file for the MLP neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H
#define LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoradlcontent/LArCheating/MLPCheatHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPLaterTierHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.h"

using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MLPNeutrinoHierarchyAlgorithm class
 */
class MLPNeutrinoHierarchyAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MLPNeutrinoHierarchyAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~MLPNeutrinoHierarchyAlgorithm();

private:
    typedef std::map<const pandora::ParticleFlowObject*, const pandora::MCParticle*> PfoToMCParticleMap;
    typedef std::map<const pandora::ParticleFlowObject*, std::pair<const pandora::ParticleFlowObject*, int>> ChildToParentPfoMap;
    typedef std::vector<pandora::PfoVector> Hierarchy;

    pandora::StatusCode Run();

    float GetSeparation(const pandora::ParticleFlowObject *const pParentPfo, const pandora::ParticleFlowObject *const pChildPfo) const;
    
    /**
     *  @brief  Return the neutrino pfo
     *
     *  @param  the pointer to the neutrino pfo to fill
     *
     *  @return whether the neutrino pfo can be found
     */
    bool GetNeutrinoPfo(const pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
     *  @brief  Fill the track and shower-like HierarchyPfoVectors
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     */
    void FillTrackShowerVectors(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, 
        HierarchyPfoMap &showerPfos) const;

    /**
     *  @brief  Get the number of 3D hits owned by a pfo
     *
     *  @param  pPfo a pointer to the pfo
     *
     *  @return The number of 3D hits
     */
    float GetNSpacepoints(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Identify the upstream (closest to the nu vertex) and downstream (furthest 
     *          from the nu vertex) endpoints and directions of an input pfo  
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  pPfo a pointer to the pfo
     *  @param  slidingFitResult the sliding fit result of the input pfo 
     *  @param  upstreamVertex the upstream endpoint
     *  @param  upstreamDirection the direction at the upstream point
     *  @param  downstreamVertex the downstream endpoint
     *  @param  downstreamDirection the direction at the downstream point
     *
     *  @return whether this was successful
     */
    bool GetExtremalVerticesAndDirections(const pandora::ParticleFlowObject *const pNeutrinoPfo, const pandora::ParticleFlowObject *const pPfo, 
        const ThreeDSlidingFitResult &slidingFitResult, pandora::CartesianVector &upstreamVertex, pandora::CartesianVector &upstreamDirection, 
        pandora::CartesianVector &downstreamVertex, pandora::CartesianVector &downstreamDirection) const;

    /**
     *  @brief  Obtain the direction of a shower at a given endpoint from
     *          the angular decomposition of the shower's 'close' 3D hits
     *
     *  @param  pPfo a pointer to the pfo
     *  @param  vertex the endpoint at which to find the direction
     *  @param  direction the found direction
     *
     *  @return whether a shower direction could be found
     */
    bool GetShowerDirection(const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &vertex, pandora::CartesianVector &direction) const;

    /**
     *  @brief  Set the primary network score of the HierarchyPfo
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector 
     */
    void SetPrimaryScores(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, HierarchyPfoMap &showerPfos) const;

    float GetPrimaryScore(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const HierarchyPfo &hierarchyPfo) const;

    void UpdateHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo, const bool buildPrimaryTier, const bool usePrimaryScore, 
        const float trackThreshold, const float showerThreshold, const bool isLowerThreshold, HierarchyPfoMap &trackPfos, 
        HierarchyPfoMap &showerPfos, Hierarchy &hierarchy) const;

    void SetLaterTierScores(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, HierarchyPfoMap &showerPfos) const;

    float GetLaterTierScore(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo) const;

    void BuildPandoraHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos) const;

    bool ShouldTrainOnEvent(const pandora::ParticleFlowObject *const pNeutrinoPfo) const;

    void GetParticleIDMap(const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos, 
        std::map<const pandora::ParticleFlowObject *, int> &particleIDMap);

    std::pair<float, float> GetTrainingCuts(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo,
        const bool trueParentOrientation, const bool trueChildOrientation) const;

    void FillEventTree(const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos, const int nPrimaryTrackLinks,
        const int nPrimaryShowerLinks, const int nLaterTierTrackTrackLinks, const int nLaterTierTrackShowerLinks) const;

    void FillPrimaryTrees(const PfoToMCParticleMap &matchingMap, const ChildToParentPfoMap &childToParentPfoMap,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos, 
        const std::map<const pandora::ParticleFlowObject *, int> &particleIDMap, int &nPrimaryTrackLinks, int &nPrimaryShowerLinks) const;

    void FillPrimaryTree(const std::string &treeName, const bool isTrainingLink, const bool isTrueLink, const bool isOrientationCorrect, 
        const int trueVisibleGen, const int trueParentID, const int particleID, const pandora::CartesianVector &upstreamVertex,
        const pandora::CartesianVector &downstreamVertex, const MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams &primaryNetworkParams) const;

    void FillLaterTierTrees(const PfoToMCParticleMap &matchingMap, const ChildToParentPfoMap &childToParentPfoMap,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos, 
        const std::map<const pandora::ParticleFlowObject *, int> &particleIDMap,  int &nTrackLinks, int &nShowerLinks) const;

    void FillLaterTierTree(const std::string &treeName, const bool isTrainingLink, const bool isTrueLink, const bool isOrientationCorrect, 
        const int childTrueGen, const std::pair<float, float> &trainingCuts, const int parentID, const int childID,
        const MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams &networkParams) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    int m_eventID;

    bool m_trainingMode;
    std::string m_trainingFileName;
    std::string m_eventTreeName;
    std::string m_primaryTrackTreeName;
    std::string m_primaryShowerTreeName;
    std::string m_laterTierTrackTrackTreeName;
    std::string m_laterTierTrackShowerTreeName;
    std::string m_mcParticleListName;
    float m_trainingVertexAccuracy;
    
    std::string m_neutrinoPfoListName;
    pandora::StringVector m_pfoListNames;
    float m_bogusFloat;
    int m_minClusterSize;
    int m_slidingFitWindow;
    float m_regionForDirFit;
    int m_nAngularBins;
    float m_primaryRegion;
    float m_primaryThresholdTrackPass1;
    float m_primaryThresholdShowerPass1;
    float m_laterTierThresholdTrackPass1;
    float m_laterTierThresholdShowerPass1;
    float m_primaryThresholdTrackPass2;
    float m_primaryThresholdShowerPass2;
    float m_laterTierThresholdTrackPass2;
    float m_laterTierThresholdShowerPass2;
    MLPPrimaryHierarchyTool *m_primaryHierarchyTool;
    MLPLaterTierHierarchyTool *m_laterTierHierarchyTool;
    MLPCheatHierarchyTool *m_cheatHierarchyTool;
};

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H
