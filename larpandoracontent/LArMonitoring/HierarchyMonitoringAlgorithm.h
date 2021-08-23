/**
 *  @file   larpandoracontent/LArMonitoring/HierarchyMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_HIERARCHY_MONITORING_ALGORITHM_H
#define LAR_HIERARCHY_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraEnumeratedTypes.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

namespace lar_content
{

/**
 *  @brief  HierarchyMonitoringAlgorithm class
 */
class HierarchyMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HierarchyMonitoringAlgorithm();

    virtual ~HierarchyMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Visualize MC nodes.
     *
     *  @param  hierarchy The MC hierarchy to render
     **/
    void VisualizeMC(const LArHierarchyHelper::MCHierarchy &hierarchy) const;

    /**
     *  @brief  Visualize MC nodes without grouping by particle id.
     *
     *  @param  hierarchy The MC hierarchy to render
     **/
    void VisualizeMCDistinct(const LArHierarchyHelper::MCHierarchy &hierarchy) const;

    /**
     *  @brief  Visualize the reco nodes.
     *
     *  @param  hierarchy The reco hierarchy to render
     **/
    void VisualizeReco(const LArHierarchyHelper::RecoHierarchy &hierarchy) const;

    /**
     *  @brief  Visualize reco to MC matches.
     *
     *  @param  matchInfo The match information between reco and MC hierarchies
     */
    void VisualizeMatches(const LArHierarchyHelper::MatchInfo &matchInfo) const;

    /**
     *  @brief  Visualize the reco nodes matched to a single MC node
     *
     *  @param  matches The MC to reco matches
     *  @param  mcIdx The unique identifier for the MC particle
     */
    void VisualizeMatchedMC(const LArHierarchyHelper::MCMatches &matches, const int mcIdx) const;

    /**
     *  @brief  Visualize the unmatched MC node
     *
     *  @param  pNode The MC without reco matches
     *  @param  mcIdx The unique identifier for the MC particle
     */
    void VisualizeUnmatchedMC(const LArHierarchyHelper::MCHierarchy::Node *pNode, const int mcIdx) const;

    /**
     *  @brief  Visualize the unmatched reco node
     *
     *  @param  pNode The unmatched reco node
     */
    void VisualizeUnmatchedReco(const LArHierarchyHelper::RecoHierarchy::Node *pNode) const;

    std::string m_caloHitListName;  ///< Name of input calo hit list
    std::string m_pfoListName;      ///< Name of input PFO list
    bool m_visualizeMC;             ///< Whether or not to visualize the MC nodes
    bool m_visualizeReco;           ///< Whether or not to visualize the reco nodes
    bool m_visualizeDistinct;       ///< If true, allocate colours without consideration of particle id
    bool m_match;                   ///< Whether or not to visualize the reco to MC matches
    float m_transparencyThresholdE; ///< Cell energy for which transparency is saturated (0%, fully opaque)
    float m_energyScaleThresholdE;  ///< Cell energy for which color is at top end of continous color palette
    float m_scalingFactor;          ///< TEve works with [cm], Pandora usually works with [mm] (but LArContent went with cm too)
};

} // namespace lar_content

#endif // LAR_HIERARCHY_MONITORING_ALGORITHM_H
