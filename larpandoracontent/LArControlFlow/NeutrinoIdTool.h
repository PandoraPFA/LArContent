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

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"
#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"

#include <functional>

namespace lar_content
{

/**
 *  @brief  NeutrinoIdTool class
 *
 *          Compares the neutrino and cosmic hypotheses of all of the slices in the event. Uses an MVA to calculate the probability of each slice
 *          containing a neutrino interaction. The N slices with the highest probabilities are identified as a neutrino (if sufficiently probable)
 *          all other slices are deemed cosmogenic.
 *
 *          If training mode is switched on, then the tool will write MVA training exmples to the specified output file. The events selected for
 *          training must pass (user configurable) slicing quality cuts. Users may also select events based on their interaction type (nuance code).
 */
template <typename T>
class NeutrinoIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoIdTool();

    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    /**
     *  @brief  Slice features class
     */
    class SliceFeatures
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  nuPfos input list of Pfos reconstructed under the neutrino hypothesis
         *  @param  crPfos input list of Pfos reconstructed under the cosmic ray hypothesis
         *  @param  pTool address of the tool using this class
         */
        SliceFeatures(const pandora::PfoList &nuPfos, const pandora::PfoList &crPfos, const NeutrinoIdTool *const pTool);

        /**
         *  @brief  Check if all features were calculable
         *
         *  @return true if the feature vector is available
         */
        bool IsFeatureVectorAvailable() const;

        /**
         *  @brief  Get the feature vector for the MVA
         *
         *  @param  featuresVector empty feature vector to populate
         */
        void GetFeatureVector(LArMvaHelper::MvaFeatureVector &featureVector) const;

        /**
         *  @brief  Get the feature map for the MVA
         *
         *  @param  featuresMap empty feature map to populate
         */
        void GetFeatureMap(std::map<std::string, double> &featureMap) const;

        /**
         *  @brief  Get the probability that this slice contains a neutrino interaction
         *
         *  @param  t the MVA used to calculate the probability
         *
         *  @return the probability that the slice contains a neutrino interaction
         */
        float GetNeutrinoProbability(const T &t) const;

    private:
        /**
         *  @brief  Get the recontructed neutrino the input list of neutrino Pfos
         *
         *  @param  nuPfos input list of neutrino pfos
         */
        const pandora::ParticleFlowObject *GetNeutrino(const pandora::PfoList &nuPfos) const;

        /**
         *  @brief  Get the 3D space points in a given pfo
         *
         *  @param  pPfo input pfo
         *  @param  spacePoints vector to hold the 3D space points associated with the input pfo
         */
        void GetSpacePoints(const pandora::ParticleFlowObject *const pPfo, pandora::CartesianPointVector &spacePoints) const;

        /**
         *  @brief  Use a sliding fit to get the direction of a collection of spacepoints
         *
         *  @param  spacePoints the input spacepoints to fit
         *  @param  fShouldChooseA a function that when given two fitted endpoints A and B, will return true if A is the endpoint at which to calculate the direction
         *
         *  @return the direction of the input spacepoints
         */
        pandora::CartesianVector GetDirection(const pandora::CartesianPointVector &spacePoints,
            std::function<bool(const pandora::CartesianVector &pointA, const pandora::CartesianVector &pointB)> fShouldChooseA) const;

        /**
         *  @brief  Use a sliding fit to get the direction of a collection of spacepoint near a vertex position
         *
         *  @param  spacePoints the input spacepoints to fit
         *  @param  vertex the position from which the fitted direction should be calculated
         *
         *  @return the direction of the input space points from the vertex supplied
         */
        pandora::CartesianVector GetDirectionFromVertex(const pandora::CartesianPointVector &spacePoints, const pandora::CartesianVector &vertex) const;

        /**
         *  @brief  Use a sliding fit to get the upper direction of a collection of spacepoints
         *
         *  @param  spacePoints the input spacepoints to fit
         *
         *  @return the direction of the upper input space points
         */
        pandora::CartesianVector GetUpperDirection(const pandora::CartesianPointVector &spacePoints) const;

        /**
         *  @brief  Use a sliding fit to get the lower direction of a collection of spacepoints
         *
         *  @param  spacePoints the input spacepoints to fit
         *
         *  @return the direction of the lower input space points
         */
        pandora::CartesianVector GetLowerDirection(const pandora::CartesianPointVector &spacePoints) const;

        /**
         *  @brief  Get a vector of spacepoints within a given radius of a vertex point
         *
         *  @param  spacePoints the input spacepoints
         *  @param  vertex the center of the sphere
         *  @param  radius the radius of the sphere
         *  @param  spacePointsInSphere the vector to hold the spacepoint in the sphere
         */
        void GetPointsInSphere(const pandora::CartesianPointVector &spacePoints, const pandora::CartesianVector &vertex, const float radius,
            pandora::CartesianPointVector &spacePointsInSphere) const;

        bool m_isAvailable;                             ///< Is the feature vector available
        LArMvaHelper::MvaFeatureVector m_featureVector; ///< The MVA feature vector
        std::map<std::string, double> m_featureMap;     ///< A map between MVA features and their names
        const NeutrinoIdTool *const m_pTool;            ///< The tool that owns this
    };

    typedef std::pair<unsigned int, float> UintFloatPair;
    typedef std::vector<SliceFeatures> SliceFeaturesVector;

    /**
     *  @brief  Get the features of each slice
     *
     *  @param  pTool the address of the this NeutrinoId tool
     *  @param  nuSliceHypotheses the input neutrino slice hypotheses
     *  @param  crSliceHypotheses the input cosmic slice hypotheses
     *  @param  sliceFeaturesVector vector to hold the slice features
     */
    void GetSliceFeatures(const NeutrinoIdTool *const pTool, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, SliceFeaturesVector &sliceFeaturesVector) const;

    /**
     *  @brief  Get the slice with the most neutrino induced hits using Monte-Carlo information
     *
     *  @param  pAlgorithm address of the master algorithm
     *  @param  nuSliceHypotheses the input neutrino slice hypotheses
     *  @param  crSliceHypotheses the input cosmic slice hypotheses
     *  @param  bestSliceIndex the index of the slice with the most neutrino hits
     *
     *  @return does the best slice pass the quality cuts for training?
     */
    bool GetBestMCSliceIndex(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, unsigned int &bestSliceIndex) const;

    /**
     *  @brief  Determine if the event passes the selection cuts for training and has the required NUANCE code
     *
     *  @param  pAlgorithm address of the master algorithm
     *  @param  purity purity of best slice
     *  @param  completeness completeness of best slice
     *
     *  @return does the evenr pass the quality cuts on purity and completeness and has the required NUANCE code
     */
    bool PassesQualityCuts(const pandora::Algorithm *const pAlgorithm, const float purity, const float completeness) const;

    /**
     *  @brief  Collect all 2D hits in a supplied list of Pfos and push them on to an existing hit list, check so not to double count
     *
     *  @param  pfos input list of pfos
     *  @param  reconstructedCaloHitList output list of all 2d hits in the input pfos
     *  @param  reconstructableCaloHitSet set of reconstructable calo hits
     */
    void Collect2DHits(const pandora::PfoList &pfos, pandora::CaloHitList &reconstructedCaloHitList,
        const pandora::CaloHitSet &reconstructableCaloHitSet) const;

    /**
     *  @brief  Count the number of neutrino induced hits in a given list using MC information
     *
     *  @param  caloHitSet input list of calo hits
     *
     *  @return the number of neutrino induced hits in the input list
     */
    unsigned int CountNeutrinoInducedHits(const pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Use the current MCParticle list to get the nuance code of the neutrino in the event
     *
     *  @param  pAlgorithm address of the master algorithm
     *
     *  @return the nuance code of the event
     */
    int GetNuanceCode(const pandora::Algorithm *const pAlgorithm) const;

    /**
     *  @brief  Select all pfos under the same hypothesis
     *
     *  @param  pAlgorithm address of the master algorithm
     *  @param  hypotheses the lists of slices under a certain hypothesis
     *  @param  selectedPfos the list of pfos to populate
     */
    void SelectAllPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &hypotheses, pandora::PfoList &selectedPfos) const;

    /**
     *  @brief  Select pfos based on the probability that their slice contains a neutrino interaction
     *
     *  @param  pAlgorithm address of the master algorithm
     *  @param  nuSliceHypotheses the input neutrino slice hypotheses
     *  @param  crSliceHypotheses the input cosmic slice hypotheses
     *  @param  sliceFeaturesVector vector holding the slice features
     *  @param  selectedPfos the list of pfos to populate
     */
    void SelectPfosByProbability(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, const SliceFeaturesVector &sliceFeaturesVector, pandora::PfoList &selectedPfos) const;

    /**
     *  @brief  Add the given pfos to the selected Pfo list
     *
     *  @param  pfos the pfos to select
     *  @param  selectedPfos the list of pfos to populate
     */
    void SelectPfos(const pandora::PfoList &pfos, pandora::PfoList &selectedPfos) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Training
    bool m_useTrainingMode;           ///< Should use training mode. If true, training examples will be written to the output file
    std::string m_trainingOutputFile; ///< Output file name for training examples
    bool m_selectNuanceCode;          ///< Should select training events by nuance code
    int m_nuance;                     ///< Nuance code to select for training
    float m_minPurity;                ///< Minimum purity of the best slice to use event for training
    float m_minCompleteness;          ///< Minimum completeness of the best slice to use event for training

    // Classification
    float m_minProbability;      ///< Minimum probability required to classify a slice as the neutrino
    unsigned int m_maxNeutrinos; ///< The maximum number of neutrinos to select in any one event

    // Persistence
    bool m_persistFeatures; ///< If true, the mva features will be persisted in the metadata

    T m_mva;                                   ///< The mva
    std::string m_filePathEnvironmentVariable; ///< The environment variable providing a list of paths to mva files
};

typedef NeutrinoIdTool<AdaBoostDecisionTree> BdtNeutrinoIdTool;
typedef NeutrinoIdTool<SupportVectorMachine> SvmNeutrinoIdTool;

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_ID_TOOL_H
