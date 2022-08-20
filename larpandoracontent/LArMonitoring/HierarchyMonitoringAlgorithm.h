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
     *  @brief  Visualize MC nodes based on the MC process that created them.
     *
     *  @param  hierarchy The MC hierarchy to render
     **/
    void VisualizeMCProcess(const LArHierarchyHelper::MCHierarchy &hierarchy) const;

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
     *  @brief  Visualize the unmatched reco node
     *
     *  @param  pNode The unmatched reco node
     */
    void VisualizeUnmatchedReco(const LArHierarchyHelper::RecoHierarchy::Node *pNode) const;

    /**
     *  @brief  Visualize a calo hit list
     *
     *  @param  hits The hits to visualize
     *  @param  label The label to apply to the hits
     *  @param  color The color to apply to the hits
     */
    void Visualize(const pandora::CaloHitList &hits, const std::string &label, const int color) const;

    /**
     *  @brief  Fill per view hit lists
     *
     *  @param  hits The input list of hits
     *  @param  uHits The output list of hits in U
     *  @param  vHits The output list of hits in V
     *  @param  wHits The output list of hits in W
     */
    void FillHitLists(const pandora::CaloHitList &hits, pandora::CaloHitList &uHits, pandora::CaloHitList &vHits, pandora::CaloHitList &wHits) const;

    std::string ToStringSF(const float val, const int sf = 3) const;

    std::string m_caloHitListName;  ///< Name of input calo hit list
    std::string m_pfoListName;      ///< Name of input PFO list
    std::string m_rootFileName;     ///< Name of the output ROOT file (optional)
    bool m_visualizeMC;             ///< Whether or not to visualize the MC nodes
    bool m_visualizeReco;           ///< Whether or not to visualize the reco nodes
    bool m_visualizeDistinct;       ///< If true, allocate colours without consideration of particle id
    bool m_visualizeProcess;        ///< If true, allocate colours based on the MC process
    bool m_match;                   ///< Whether or not to visualize the reco to MC matches
    bool m_collectionOnly;          ///< Limit display to the collection plane only
    bool m_foldToPrimaries;         ///< Whether or not to fold everything back to primaries
    bool m_foldDynamic;             ///< Whether or not to fold based on process information
    float m_transparencyThresholdE; ///< Cell energy for which transparency is saturated (0%, fully opaque)
    float m_energyScaleThresholdE;  ///< Cell energy for which color is at top end of continous color palette
    float m_scalingFactor;          ///< TEve works with [cm], Pandora usually works with [mm] (but LArContent went with cm too)
};

} // namespace lar_content

#endif // LAR_HIERARCHY_MONITORING_ALGORITHM_H
