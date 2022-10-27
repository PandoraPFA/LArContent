/**
 *  @file   larpandoracontent/LArCheating/CheatingSliceSelectionTool.h
 *
 *  @brief  Header file for the cheating slice selection tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_SLICE_SELECTION_TOOL_H
#define LAR_CHEATING_SLICE_SELECTION_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/SliceSelectionBaseTool.h"

namespace lar_content
{

/**
 *  @brief  CheatingSliceSelectionTool class
 */
class CheatingSliceSelectionTool : public SliceSelectionBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingSliceSelectionTool();

    /**
     *  @brief  Select which slice(s) to use
     *
     *  @param  pAlgorithm the address of the master instance, used to access MCParticles when in training mode
     *  @param  inputSliceVector the initial slice vector
     *  @param  outputSliceVector the output slice vector
     */
    void SelectSlices(const pandora::Algorithm *const pAlgorithm, const SliceVector &inputSliceVector, SliceVector &outputSliceVector);

    typedef std::map<float, int, std::greater<float>> MetricSliceIndexMap;

protected:
    /**
     *  @brief  Template method to determine if an MC particle matches the target criteria for slice selection. Return true if match.
     *
     *  @param  mcParticle the MC particle to check
     */
    virtual bool IsTarget(const pandora::MCParticle *const mcParticle) const = 0;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

protected:
    int m_maxSlices;           ///< The maximum number of slices to retain (0 to retain all) - default 0
    float m_threshold;         ///< The minimum cut threshold to retain a slice (< 0 for no threshold) - default -1
    std::string m_cutVariable; ///< The variable to cut on ("purity" or "completeness") - default "completeness"
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SLICE_SELECTION_TOOL_H
