/**
 *  @file   larpandoracontent/LArMetrics/EventValidationTool.h
 *
 *  @brief  Header file for the track validation tool class.
 *
 *  $Log: $
 */
#ifndef EVENT_VALIDATION_TOOL_H
#define EVENT_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  EventValidationTool class
 */
class EventValidationTool : public BaseValidationTool
{
public:
/**
 *  @brief  EventTreeVars struct
 */
struct EventTreeVars
{
    /**
     *  @brief  Default constructor
     */
    EventTreeVars();

    int m_run;           ///< run number
    int m_subrun;        ///< subrun number
    int m_event;         ///< event number
    int m_nTargets;      ///< number of 'reconstructable' MCParticles
    int m_nuPDG;         ///< nu PDG
    float m_nuEnergy;    ///< true nu energy [GeV]
    float m_nuVisEnergy; ///< true nu visible energy [GeV]
    int m_isCC;          ///< whether CC interaction 
    pandora::CartesianVector m_trueNuVertex;      ///< true neutrino vertex
    pandora::CartesianVector m_recoNuVertexPass1; ///< pass 1 reco neutrino vertex 
    pandora::CartesianVector m_recoNuVertexPass2; ///< pass 2 reco neutrino vertex 
    float m_recoNuVertexAccPass1; ///< accuracy of the pass 1 reco neutrino vertex
    float m_recoNuVertexAccPass2; ///< accuracy of the pass 1 reco neutrino vertex
};
    /**
     *  @brief  Default constructor
     */
    EventValidationTool();

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Fill neutrino vertex variables
     *
     *  @param  pAlgorithm a pointer to the parent algorithm 
     *  @param  pMCNu a pointer to the neutrino MCParticle
     *  @param  eventTreeVars the event tree variables to fill
     */
    void GetVertexVariables(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, EventTreeVars &eventTreeVars);

    /**
     *  @brief  Fill neutrino interaction variables
     *
     *  @param  pMCNu a pointer to the neutrino MCParticle
     *  @param  eventTreeVars the event tree variables to fill
     */
    void GetInteractionTypeVariables(const pandora::MCParticle *const pMCNu, EventTreeVars &eventTreeVars);

    /**
     *  @brief  Fill the event tree
     *
     *  @param  eventTreeVars the event tree variables to fill
     */
    void FillTree(EventTreeVars &eventTreeVars);

    std::string m_nuVertexPass1ListName; ///< name of the pass 1 neutrino vertex list
    std::string m_nuVertexPass2ListName; ///< name of the pass 2 neutrino vertex list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline EventValidationTool::EventTreeVars::EventTreeVars() :
    m_run(-1),
    m_subrun(-1),
    m_event(-1),
    m_nTargets(-1),
    m_nuPDG(-1),
    m_nuEnergy(-1.f),
    m_nuVisEnergy(-1.f),
    m_isCC(-1),
    m_trueNuVertex(pandora::CartesianVector(-9999.f, -9999.f, -9999.f)),
    m_recoNuVertexPass1(pandora::CartesianVector(-9999.f, -9999.f, -9999.f)),
    m_recoNuVertexPass2(pandora::CartesianVector(-9999.f, -9999.f, -9999.f)),
    m_recoNuVertexAccPass1(-1.f),
    m_recoNuVertexAccPass2(-1.f)
{
}

} // namespace lar_content

#endif // #ifndef EVENT_VALIDATION_TOOL_H
