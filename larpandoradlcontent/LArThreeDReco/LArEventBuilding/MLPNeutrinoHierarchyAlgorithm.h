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
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
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


    pandora::StatusCode Run();

    void DetermineIsobelID();

    void PrintHierarchy();

    bool GetNeutrinoPfo();

    void FillTrackShowerVectors();

    bool GetExtremalVerticesAndDirections(const pandora::ParticleFlowObject *const pPfo, pandora::CartesianVector &upstreamVertex, 
        pandora::CartesianVector &upstreamDirection, pandora::CartesianVector &downstreamVertex, pandora::CartesianVector &downstreamDirection);

    bool GetShowerDirection(const pandora::ParticleFlowObject *const pPfp, const pandora::CartesianVector &vertex, const float searchRegion, 
        pandora::CartesianVector &direction);

    void SetPrimaryScores();

    void BuildPrimaryTierPass1();

    void SetLaterTierScores();

    void BuildLaterTierPass1();

    float GetLaterTierScoreTrackToTrack(const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo, 
        int &parentOrientation, int &childOrientation) const;
    float GetLaterTierScoreTrackToShower(const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo,
        int &parentOrientation, int &childOrientation) const;

    void BuildPandoraHierarchy();

    void PrintPandoraHierarchy();

    void CheckForOrphans();

    float GetRandomNumber() const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //typedef std::vector<PfoRelationTool *> PfoRelationToolVector;

    std::string m_neutrinoPfoListName;
    pandora::StringVector m_pfoListNames;

    MLPPrimaryHierarchyTool *m_primaryHierarchyTool;

    float m_primaryThresholdTrackPass1;
    float m_primaryThresholdShowerPass1;
    float m_laterTierThresholdTrackPass1;
    float m_laterTierThresholdShowerPass1;

    ///////////////////////////
    std::map<const pandora::ParticleFlowObject*, int> m_isobelID;
    ///////////////////////////

    const pandora::ParticleFlowObject* m_pNeutrinoPfo;
    HierarchyPfoMap m_trackPfos;
    HierarchyPfoMap m_showerPfos;
    std::vector<pandora::PfoVector> m_hierarchy;
};

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H
