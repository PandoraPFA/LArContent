/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
#define LAR_GAMMA_START_REFINEMENT_TOOL_H 1

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

namespace lar_content
{

/**
 *  @brief  GammaStartRefinementTool class
 */
class GammaStartRefinementTool : public ShowerStartRefinementBaseTool
{
public:
    GammaStartRefinementTool();
    ~GammaStartRefinementTool();

    typedef std::map<int, ProtoShowerVector> MatchedConnectionPathwayMap;

    bool Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition,
        const pandora::CaloHitList *const pCaloHitListU, const pandora::CaloHitList *const pCaloHitListV, const pandora::CaloHitList *const pCaloHitListW);

    void BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType tpcView, ProtoShowerVector &protoShowerVector, pandora::CaloHitList &usedHitList);

    void BuildHelperProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const pandora::CartesianVector &nuVertexPosition, const pandora::HitType tpcView, ProtoShowerVector &protoShowerVector, pandora::CaloHitList &usedHitList);

    void FillOutPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CaloHitList &showerPfoHits, 
        pandora::CaloHitList &unavailableHits, ProtoShowerVector &protoShowerVector);

    bool IsShowerTruncatable(const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
        const ProtoShowerVector &protoShowerVectorW);

    void AssignShowerHits(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const pandora::HitType &hitType, ProtoShowerVector &protoShowerVector);

    void MatchConnectionPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &nuVertexPosition, 
        const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
        const ProtoShowerVector &protoShowerVectorW, MatchedConnectionPathwayMap &threeViewConnectionPathways, MatchedConnectionPathwayMap &twoViewConnectionPathways, 
        MatchedConnectionPathwayMap &oneViewConnectionPathways);

    void MatchShowerStart(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
        const ProtoShowerVector &protoShowerVectorW, pandora::IntVector &matchedProtoShowersU, pandora::IntVector &matchedProtoShowersV, 
        pandora::IntVector &matchedProtoShowersW, MatchedConnectionPathwayMap &matchedConnectionPathwayMap);

    void MatchShowerStart2(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
        const ProtoShowerVector &protoShowerVectorW, pandora::IntVector &matchedProtoShowersU, pandora::IntVector &matchedProtoShowersV, 
        pandora::IntVector &matchedProtoShowersW, MatchedConnectionPathwayMap &matchedConnectionPathwayMap);

    void MatchPeakDirection(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
        const ProtoShowerVector &protoShowerVectorW, pandora::IntVector &matchedProtoShowersU, pandora::IntVector &matchedProtoShowersV,  
        pandora::IntVector &matchedProtoShowersW, MatchedConnectionPathwayMap &matchedConnectionPathwayMap);

    //void MatchShowerStartTwoViews(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVectorU, 
    //const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW, pandora::IntVector &matchedProtoShowersU, 
    //pandora::IntVector &matchedProtoShowersV, pandora::IntVector &matchedProtoShowersW, MatchedConnectionPathwayMap &matchedConnectionPathwayMap);

    //bool MatchShowerStartTwoViews(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVector1, 
    //const ProtoShowerVector &protoShowerVector2, pandora::IntVector &matchedProtoShowers1, pandora::IntVector &matchedProtoShowers2, int &bestProtoShower1, 
    //int &bestProtoShower2, float &bestSeparation);

    //void MatchPeakDirectionTwoViews(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &nuVertexPosition, const pandora::ParticleFlowObject *const pShowerPfo, 
    //const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW, 
    //pandora::IntVector &matchedProtoShowersU, pandora::IntVector &matchedProtoShowersV, pandora::IntVector &matchedProtoShowersW, 
    //MatchedConnectionPathwayMap &matchedConnectionPathwayMap);

    //bool MatchPeakDirectionTwoViews(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::CartesianVector &nuVertexPosition, 
    //const pandora::ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVector1, const ProtoShowerVector &protoShowerVector2, 
    //pandora::IntVector &matchedProtoShowers1, pandora::IntVector &matchedProtoShowers2, int &bestProtoShower1, int &bestProtoShower2, unsigned int &bestHitCount, 
    //float &bestClosestDistance);

    void AddUnmatched(const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW, 
        pandora::IntVector &matchedProtoShowersU, pandora::IntVector &matchedProtoShowersV, pandora::IntVector &matchedProtoShowersW, 
        MatchedConnectionPathwayMap &matchedConnectionPathwayMap);

    void AddUnmatched(const ProtoShowerVector &protoShowerVector, pandora::IntVector &matchedProtoShowers, MatchedConnectionPathwayMap &matchedConnectionPathwayMap);

    void RemoveThreeViewConnectionPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const MatchedConnectionPathwayMap &threeViewConnectionPathways, 
        const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition);

    bool FindShowerVertex(ShowerStartRefinementAlgorithm *const pAlgorithm, const MatchedConnectionPathwayMap &threeViewConnectionPathways, const pandora::IntVector &toRefineIndices, 
        const pandora::CartesianVector &nuVertexPosition, const pandora::ParticleFlowObject *const pShowerPfo, pandora::CartesianVector &showerVertex);

    void RemoveOneViewConnectionPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const MatchedConnectionPathwayMap &oneViewConnectionPathways, 
        const pandora::ParticleFlowObject *const pShowerPfo);

    void FillTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition);

    void AngularDistributionTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const AngularDecompositionMap &angularDecompositionMapU, 
        const AngularDecompositionMap &angularDecompositionMapV, const AngularDecompositionMap &angularDecompositionMapW);

    bool FindShowerVertexFromPosition(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShowerVector &protoShowerVector, 
        const pandora::CartesianVector &nuVertexPosition, pandora::CartesianVector &showerStart3D);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void RemoveConnectionPathway(const ProtoShower &protoShower);

    int m_counter;
    float m_maxXSeparation;
    float m_maxSeparation;
    float m_maxAngularDeviation;
    float m_maxProjectedShowerStartSeparation;
};

} // namespace lar_content

#endif // #ifndef LAR_GAMMA_START_REFINEMENT_TOOL_H
