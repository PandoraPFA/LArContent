/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringTool.h
 *
 *  @brief  Header file for the slice monitoring tool.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_MONITORING_TOOL_H
#define LAR_SLICE_MONITORING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/SlicingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  SliceMonitoringTool class
 */
class SliceMonitoringTool : public pandora::AlgorithmTool
/**
 * @class SliceMonitoringTool
 * @brief A class for monitoring and processing slices in the LAr monitoring tool.
 */
{
public:
    /**
     * @brief Default constructor
     */
    SliceMonitoringTool();

    /**
     * @brief Destructor
     */
    virtual ~SliceMonitoringTool();
    
    /**
     * @brief Process the slices to produce monitoring information, and possibly training examples.
     * 
     * @param pAlgorithm The algorithm used for processing.
     * @param inputSliceList The list of input slices to be processed.
     */
    void ProcessSlices(const pandora::Algorithm *const pAlgorithm, SlicingAlgorithm::SliceList &inputSliceList);

private:
    /**
     * @brief Write out the hits from the input slice list using the MC contribution map.
     * 
     * @param inputSliceList The list of input slices.
     * @param mcToTrueHitListMap The map of MC contributions to true hit lists.
     */
    void WriteOutHits(SlicingAlgorithm::SliceList &inputSliceList, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    bool m_trainingMode;              ///< Training mode toggle.
    std::string m_filename;           ///< The filename of the ROOT output file
    std::string m_treename;           ///< The name of the ROOT tree
    std::string m_trainingOutputFile; ///< Output name for training examples.
};

} // namespace lar_content

#endif // LARe_SLICE_MONITORING_TOOL_H
