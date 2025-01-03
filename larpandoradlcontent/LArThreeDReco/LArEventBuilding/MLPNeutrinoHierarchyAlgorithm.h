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

private:
    typedef std::vector<pandora::PfoVector> Hierarchy;

    pandora::StatusCode Run();

    void DetermineIsobelID();

    void PrintHierarchy(const Hierarchy &hierarchy) const;

    bool GetNeutrinoPfo(const pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    void FillTrackShowerVectors(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, HierarchyPfoMap &showerPfos) const;

    float GetNSpacepoints(const pandora::ParticleFlowObject *const pPfo) const;

    bool GetExtremalVerticesAndDirections(const pandora::ParticleFlowObject *const pNeutrinoPfo, const pandora::ParticleFlowObject *const pPfo, 
        const ThreeDSlidingFitResult &slidingFitResult, pandora::CartesianVector &upstreamVertex, pandora::CartesianVector &upstreamDirection, 
        pandora::CartesianVector &downstreamVertex, pandora::CartesianVector &downstreamDirection) const;

    bool GetShowerDirection(const pandora::ParticleFlowObject *const pPfp, const pandora::CartesianVector &vertex, pandora::CartesianVector &direction) const;

    void SetPrimaryScores(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, HierarchyPfoMap &showerPfos) const;

    float GetPrimaryScore(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const HierarchyPfo &hierarchyPfo) const;

    void UpdateHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo, const bool buildPrimaryTier, const bool usePrimaryScore, 
        const float trackThreshold, const float showerThreshold, const bool isLowerThreshold, HierarchyPfoMap &trackPfos, 
        HierarchyPfoMap &showerPfos, Hierarchy &hierarchy) const;

    void SetLaterTierScores(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, HierarchyPfoMap &showerPfos) const;

    float GetLaterTierScore(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo) const;

    void BuildPandoraHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos) const;

    void PrintPandoraHierarchy(const pandora::ParticleFlowObject *const pNeutrinoPfo) const;

    void CheckForOrphans() const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_neutrinoPfoListName;
    pandora::StringVector m_pfoListNames;
    bool m_trainingMode;
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

    ///////////////////////////
    std::map<const pandora::ParticleFlowObject*, int> m_isobelID;
    ///////////////////////////
};

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H
