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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void FillShowerPfoVector(pandora::PfoVector &showerPfoVector);

    void RefineShower(const pandora::ParticleFlowObject *const pShowerPfo);

    pandora::StatusCode GetNeutrinoVertex(pandora::CartesianVector &nuVertex3D);

    void BuildViewProtoShowers(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, 
        pandora::HitType hitType, ElectronProtoShowerVector &protoShowerVector);

    pandora::StatusCode GetHitListOfType(const pandora::HitType hitType, const pandora::CaloHitList *&pCaloHitList) const;
    
    pandora::CartesianVector GetShowerVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, 
        const pandora::CartesianVector &nuVertex3D) const;

    void RefineShowerVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::HitType hitType, 
        const pandora::CartesianVector &nuVertex3D, const pandora::CartesianVector &peakDirection, pandora::CartesianVector &showerVertexPosition) const;

    bool IsSpineCoincident(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, 
        const pandora::HitType hitType, const pandora::CartesianVector &showerVertex, const pandora::CaloHitList &showerSpineHitList) const;

    bool IsShowerConnected(const pandora::CartesianVector &showerVertexPosition, const pandora::CartesianVector &nuVertex2D, 
        const pandora::CartesianVector &peakDirection) const;

    void BuildViewPathways(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CaloHitList &protectedHits, const pandora::CartesianVector &nuVertex3D,
        pandora::HitType hitType, ConnectionPathwayVector &viewPathways) const;

    void RefineHitsToAdd(const pandora::CartesianVector &nuVertex3D, const pandora::HitType hitType, const ConnectionPathwayVector &viewPathways, 
        ElectronProtoShower &protoShower) const;

    pandora::CaloHitList FindContinuousPath(const pandora::CaloHitList &refinedHitList, const pandora::CartesianVector &nuVertex2D) const;

    PeakDirectionFinderTool* m_pShowerPeakDirectionFinderTool;
    PeakDirectionFinderTool* m_pEventPeakDirectionFinderTool;
    ShowerSpineFinderTool* m_pShowerSpineFinderTool;
    ShowerSpineFinderTool* m_pEventPathwayFinderTool;
    ShowerStartFinderTool* m_pShowerStartFinderTool;
    ProtoShowerMatchingTool* m_pProtoShowerMatchingTool;

    std::string m_showerPfoListName;
    std::string m_neutrinoVertexListName;
    std::string m_caloHitListNameU;
    std::string m_caloHitListNameV;
    std::string m_caloHitListNameW;

    unsigned int m_showerSlidingFitWindow;
    float m_maxCoincideneTransverseSeparation;
    float m_minSpinePurity;
};

} // namespace lar_content

#endif // #ifndef LAR_ELECTRON_INITIAL_REGION_REFINEMENT
