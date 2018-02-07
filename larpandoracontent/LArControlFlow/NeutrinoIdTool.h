/**
 *  @file   larpandoracontent/LArControlFlow/NeutrinoIdTool.h
 *
 *  @brief  Header file for the neutrino id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_ID_TOOL_H
#define LAR_NEUTRINO_ID_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "Pandora/PdgTable.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoIdTool class
 *
 *  Compares the neutrino and cosmic hypotheses of all of the slices in the event. Uses an SVM to calculate the probability of each slice 
 *  containing a neutrino interaction. If the slice with the highest probability is sufficiently probable, then identify it as the neutrino 
 *  slice, otherwise assume all slices are cosmogenic.
 *
 *  If training mode is switched on, then the tool will write SVM training exmples to the specified output file. The events selected for 
 *  training must pass (user congigurable) slicing quality cuts. Users may also select events based on their interaction type (nuance code).
 */
class NeutrinoIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoIdTool();

    void SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    /**
     *  @brief  Get the slice with the most neutrino induced hits
     */
    unsigned int GetBestSliceIndex(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, float &purity, float &completeness) const;

    /**
     *  @brief  Collect all 2D hits in a supplied list of Pfos and push them on to an existing hit list, check so not to double count
     */
    void Collect2DHits(const pandora::PfoList &pfos, pandora::CaloHitList &hitList) const;

    /**
     *  @brief  Count the number of neutrino induced hits in a given list
     */
    unsigned int CountNeutrinoInducedHits(const pandora::CaloHitList &hitList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Training
    bool         m_useTrainingMode;     // Should use training mode. If true, training examples will be written to the output file
    std::string  m_trainingOutputFile;  // Output file name for training examples
    bool         m_selectNuanceCode;    // Should select training events by nuance code
    int          m_nuance;              // Nuance code to select for training
    float        m_minPurity;           // Minimum purity of the best slice to use event for training
    float        m_minCompleteness;     // Minimum completeness of the best slice to use event for training

    // Classification
    float        m_minProbability;      // Minimum probability required to classify a slice as the neutrino
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_ID_TOOL_H
