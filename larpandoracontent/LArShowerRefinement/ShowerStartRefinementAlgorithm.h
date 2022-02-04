/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h
 *
 *  @brief  Header file for the shower start refinement algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_START_REFINEMENT_ALGORITHM_H
#define LAR_SHOWER_START_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

namespace lar_content
{

class ShowerStartRefinementBaseTool;

class ShowerStartRefinementAlgorithm : public PfoMopUpBaseAlgorithm
{
public:
    ShowerStartRefinementAlgorithm();
    ~ShowerStartRefinementAlgorithm();

    typedef std::map<int, int> DeviationAngleMap;
    DeviationAngleMap m_thetaMapU;
    DeviationAngleMap m_thetaMapV;
    DeviationAngleMap m_thetaMapW;

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;
    typedef std::unordered_map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    HitToClusterMap m_hitToClusterMapU;
    HitToClusterMap m_hitToClusterMapV;
    HitToClusterMap m_hitToClusterMapW;
    ClusterToPfoMap m_clusterToPfoMapU;
    ClusterToPfoMap m_clusterToPfoMapV;
    ClusterToPfoMap m_clusterToPfoMapW;

    pandora::PfoList m_deletedPfos;

    float m_binSize;

    void FillOwnershipMaps();

    bool IsElectronPathway(const pandora::CaloHitList &hitsToAdd);
    void AddElectronPathway(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CaloHitList &pathwayHitList);

private:
    pandora::StatusCode Run();
    void FillPfoVector(pandora::PfoVector &pfoVector);
    pandora::StatusCode GetNeutrinoVertex(pandora::CartesianVector &neutrinoVertex);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;
    std::string m_neutrinoVertexListName;

    typedef std::vector<ShowerStartRefinementBaseTool *> ShowerStartRefinementToolVector;
    ShowerStartRefinementToolVector m_algorithmToolVector;
    //private:


};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
