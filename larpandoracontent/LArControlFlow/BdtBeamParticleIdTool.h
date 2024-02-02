/**
 *  @file   larpandoracontent/LArControlFlow/BdtBeamParticleIdTool.h
 *
 *  @brief  Header file for the beam particle id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_BDT_BEAM_PARTICLE_ID_TOOL_H
#define LAR_BDT_BEAM_PARTICLE_ID_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/SliceIdBaseTool.h"

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"

namespace lar_content
{

/**
 *  @brief  BdtBeamParticleIdTool class
 */
class BdtBeamParticleIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Constructor
     */
    BdtBeamParticleIdTool();

    /**
     *  @brief  Copy constructor
     *
     *  @param  BdtBeamParticleIdTool to copy
     */
    BdtBeamParticleIdTool(const BdtBeamParticleIdTool &) = default;

    /**
     *  @brief  Assignment operator
     *
     *  @param  The BdtBeamParticleIdTool to assign
     */
    BdtBeamParticleIdTool &operator=(const BdtBeamParticleIdTool &) = default;

    /**
     *  @brief  Destructor
     */
    ~BdtBeamParticleIdTool() = default;

    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &beamSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    /**
     *  @brief  Plane class
     */
    class Plane
    {
    public:
        /**
         *  @brief  Constructor, using equation of plane: m_a*x + m_b*y + m_c*z + m_d = 0.
         *
         *  @param  normal a Cartesian vector that points in a direction that is normal to the plane
         *  @param  point a Cartesian vector that corresponds to any point on the plane
         */
        Plane(const pandora::CartesianVector &normal, const pandora::CartesianVector &point);

        /**
         *  @brief  Return the intersection between the plane and a line
         *
         *  @param  point a point on the line
         *  @param  direction vector pointing along the line
         */
        pandora::CartesianVector GetLineIntersection(const pandora::CartesianVector &point, const pandora::CartesianVector &direction) const;

    private:
        pandora::CartesianVector m_unitNormal; ///< Unit normal to plane
        pandora::CartesianVector m_point;      ///< A point on the plane
        double m_d;                            ///< Parameter defining a plane
    };

    typedef std::vector<Plane> PlaneVector;

    /**
     *  @brief  SliceFeatureParameters class
     */
    class SliceFeatureParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        SliceFeatureParameters();

        /**
         *  @brief  Set LArTPC geometry information
         *
         *  @param  larTPCMinX the global LArTPC volume minimum x extent
         *  @param  larTPCMaxX the global LArTPC volume maximum x extent
         *  @param  larTPCMinY the global LArTPC volume minimum y extent
         *  @param  larTPCMaxY the global LArTPC volume maximum y extent
         *  @param  larTPCMinZ the global LArTPC volume minimum z extent
         *  @param  larTPCMaxZ the global LArTPC volume maximum z extent
         */
        void Initialize(const float larTPCMinX, const float larTPCMaxX, const float larTPCMinY, const float larTPCMaxY,
            const float larTPCMinZ, const float larTPCMaxZ);

        /**
         *  @brief  Set m_beamLArTPCIntersection
         *
         *  @param  beamLArTPCIntersection intersection of beam and TPC
         */
        void SetBeamLArTPCIntersection(const pandora::CartesianVector &beamLArTPCIntersection);

        /**
         *  @brief  Set m_beamDirection
         *
         *  @param  beamDirection direction of beam
         */
        void SetBeamDirection(const pandora::CartesianVector &beamDirection);

        /**
         *  @brief  Set m_selectedFraction
         *
         *  @param  selectedFraction fraction of hits to use in 3D cluster fits
         */
        void SetSelectedFraction(const float selectedFraction);

        /**
         *  @brief  Set m_nSelectedHits
         *
         *  @param  nSelectedHits minimum number of hits to use in 3D cluster fits
         */
        void SetNSelectedHits(const unsigned int nSelectedHits);

        /**
         *  @brief  Set m_containmentLimit
         *
         *  @param  containmentLimit limit used in is contained definition
         */
        void SetContainmentLimit(const float containmentLimit);

        /**
         *  @brief  Get m_larTPCMinX
         */
        float GetLArTPCMinX() const;

        /**
         *  @brief  Get m_larTPCMaxX
         */
        float GetLArTPCMaxX() const;

        /**
         *  @brief  Get m_larTPCMinY
         */
        float GetLArTPCMinY() const;

        /**
         *  @brief  Get m_larTPCMaxY
         */
        float GetLArTPCMaxY() const;

        /**
         *  @brief  Get m_larTPCMinZ
         */
        float GetLArTPCMinZ() const;

        /**
         *  @brief  Get m_larTPCMaxZ
         */
        float GetLArTPCMaxZ() const;

        /**
         *  @brief  Get vector of planes
         */
        const PlaneVector &GetPlanes() const;

        /**
         *  @brief  Get the beam LArTPC intersection
         */
        const pandora::CartesianVector &GetBeamLArTPCIntersection() const;

        /**
         *  @brief  Get the beam direction
         */
        const pandora::CartesianVector &GetBeamDirection() const;

        /**
         *  @brief  Get m_selectedFraction
         */
        float GetSelectedFraction() const;

        /**
         *  @brief  Get m_nSelectedHits
         */
        unsigned int GetNSelectedHits() const;

        /**
         *  @brief  Get m_containmentLimit
         */
        float GetContainmentLimit() const;

    private:
        float m_larTPCMinX;                                ///< Global LArTPC volume minimum x extent
        float m_larTPCMaxX;                                ///< Global LArTPC volume maximum x extent
        float m_larTPCMinY;                                ///< Global LArTPC volume minimum y extent
        float m_larTPCMaxY;                                ///< Global LArTPC volume maximum y extent
        float m_larTPCMinZ;                                ///< Global LArTPC volume minimum z extent
        float m_larTPCMaxZ;                                ///< Global LArTPC volume maximum z extent
        PlaneVector m_larTPCPlanes;                        ///< Vector of all planes making up global LArTPC volume
        pandora::CartesianVector m_beamLArTPCIntersection; ///< Intersection of beam and global LArTPC volume
        pandora::CartesianVector m_beamDirection;          ///< Beam direction
        float m_selectedFraction;                          ///< Fraction of hits to use in 3D cluster fits
        unsigned int m_nSelectedHits;                      ///< Minimum number of hits to use in 3D cluster fits
        float m_containmentLimit;                          ///< Limit applied in is contained definition
    };

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
         *  @param  geometryInfo geometry information block
         */
        SliceFeatures(const pandora::PfoList &nuPfos, const pandora::PfoList &crPfos, const SliceFeatureParameters &sliceFeatureParameters);

        /**
         *  @brief  Copy constructor
         *
         *  @param  The SliceFeatures to copy
         */
        SliceFeatures(const SliceFeatures &) = default;

        /**
         *  @brief  Destructor
         */
        ~SliceFeatures() = default;

        /**
         *  @brief  Check if all features were calculable
         *
         *  @return true if the feature vector is available
         */
        bool IsFeatureVectorAvailable() const;

        /**
         *  @brief  Get the feature vector for the SVM
         *
         *  @param  featuresVector empty feature vector to populate
         */
        void FillFeatureVector(LArMvaHelper::MvaFeatureVector &featureVector) const;

        /**
         *  @brief  Get the probability that this slice contains a beam particle
         *
         *  @param  adaBoostDecisionTree the adaptive boost decision tree used to calculate the probability
         *
         *  @return the probability that the slice contains a beam particle
         */
        float GetAdaBoostDecisionTreeScore(const AdaBoostDecisionTree &adaBoostDecisionTree) const;

    private:
        /**
         *  @brief  Select a given fraction of a slice's calo hits that are closest to the beam spot
         *
         *  @param  inputCaloHitList all calo hits in slice
         *  @param  outputCaloHitList to receive the list of selected calo hits
         *  @param  closestHitToFaceDistance to receive the distance of closest hit to beam spot
         */
        void GetLeadingCaloHits(
            const pandora::CaloHitList &inputCaloHitList, pandora::CaloHitList &outputCaloHitList, double &closestHitToFaceDistance) const;

        /**
         *  @brief  Find the intercepts of a line with the protoDUNE detector
         *
         *  @param  a0 a point on the line in question
         *  @param  majorAxis the direction of the line in question
         *  @param  interceptOne to receive the first intersection between line and protoDUNE detector
         *  @param  interceptTwo to receive the second intersection between line and protoDUNE detector
         */
        void GetLArTPCIntercepts(const pandora::CartesianVector &a0, const pandora::CartesianVector &majorAxis,
            pandora::CartesianVector &interceptOne, pandora::CartesianVector &interceptTwo) const;

        /**
         *  @brief  Check if a given 3D spacepoint is inside the global LArTPC volume
         *
         *  @param  spacePoint
         */
        bool IsContained(const pandora::CartesianVector &spacePoint, const float limit) const;

        bool m_isAvailable;                                    ///< Is the feature vector available
        const SliceFeatureParameters m_sliceFeatureParameters; ///< Geometry information block
        LArMvaHelper::MvaFeatureVector m_featureVector;        ///< The MVA feature vector
    };

    typedef std::vector<SliceFeatures> SliceFeaturesVector;
    typedef std::pair<unsigned int, float> UintFloatPair;
    typedef std::unordered_map<const pandora::MCParticle *, int> MCParticleToIntMap;

    pandora::StatusCode Initialize();

    /**
     *  @brief  Get the features of each slice
     *
     *  @param  nuSliceHypotheses the input neutrino slice hypotheses
     *  @param  crSliceHypotheses the input cosmic slice hypotheses
     *  @param  sliceFeaturesVector vector to hold the slice features
     */
    void GetSliceFeatures(
        const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, SliceFeaturesVector &sliceFeaturesVector) const;

    /**
     *  @brief  Select all pfos under the same hypothesis
     *
     *  @param  pAlgorithm address of the master algorithm
     *  @param  hypotheses the lists of slices under a certain hypothesis
     *  @param  selectedPfos the list of pfos to populate
     */
    void SelectAllPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &hypotheses, pandora::PfoList &selectedPfos) const;

    /**
     *  @brief  Add the given pfos to the selected Pfo list
     *
     *  @param  pfos the pfos to select
     *  @param  selectedPfos the list of pfos to populate
     */
    void SelectPfos(const pandora::PfoList &pfos, pandora::PfoList &selectedPfos) const;

    /**
     *  @brief  Get the slice with the most neutrino induced hits using Monte-Carlo information
     *
     *  @param  pAlgorithm address of the master algorithm
     *  @param  nuSliceHypotheses the input neutrino slice hypotheses
     *  @param  crSliceHypotheses the input cosmic slice hypotheses
     *  @param  bestSliceIndices vector of slice indices passing quality cuts
     */
    void GetBestMCSliceIndices(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, pandora::IntVector &bestSliceIndices) const;

    /**
     *  @brief  Fill mc particle to nHits map from calo hit list
     *
     *  @param  mcParticleToIntMap map to fill
     *  @param  caloHitList the input calo hits
     */
    void PopulateMCParticleToHitsMap(MCParticleToIntMap &mcParticleToIntMap, const pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Collect all 2D hits in a supplied list of Pfos and push them on to an existing hit list, check so not to double count
     *
     *  @param  pfos input list of pfos
     *  @param  caloHitList output list of all 2d hits in the input pfos
     *  @param  reconstructableCaloHitSet to check if part of before adding
     */
    void Collect2DHits(const pandora::PfoList &pfos, pandora::CaloHitList &caloHitList, const pandora::CaloHitSet &reconstructableCaloHitSet) const;

    /**
     *  @brief  Determine if the event passes the selection cuts for training
     *
     *  @param  purity purity of best slice
     *  @param  completeness completeness of best slice
     *
     *  @return does the evenr pass the quality cuts on purity and completeness
     */
    bool PassesQualityCuts(const float purity, const float completeness) const;

    /**
     *  @brief  Select pfos based on the AdaBDT score that the slice contains a beam particle interaction
     *
     *  @param  pAlgorithm address of the master algorithm
     *  @param  nuSliceHypotheses the input neutrino slice hypotheses
     *  @param  crSliceHypotheses the input cosmic slice hypotheses
     *  @param  sliceFeaturesVector vector holding the slice features
     *  @param  selectedPfos the list of pfos to populate
     */
    void SelectPfosByAdaBDTScore(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, const SliceFeaturesVector &sliceFeaturesVector, pandora::PfoList &selectedPfos) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Training
    bool m_useTrainingMode;           ///< Should use training mode. If true, training examples will be written to the output file
    std::string m_trainingOutputFile; ///< Output file name for training examples
    std::string m_caloHitListName;    ///< Name of input calo hit list
    std::string m_mcParticleListName; ///< Name of input MC particle list
    float m_minPurity;                ///< Minimum purity of the best slice to use event for training
    float m_minCompleteness;          ///< Minimum completeness of the best slice to use event for training

    // Classification
    AdaBoostDecisionTree m_adaBoostDecisionTree;     ///< The adaptive boost decision tree
    std::string m_filePathEnvironmentVariable;       ///< The environment variable providing a list of paths to bdt files
    unsigned int m_maxNeutrinos;                     ///< The maximum number of neutrinos to select in any one event
    float m_minAdaBDTScore;                          ///< Minimum score required to classify a slice as a beam particle
    SliceFeatureParameters m_sliceFeatureParameters; ///< Geometry information block
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParameters::SetBeamLArTPCIntersection(const pandora::CartesianVector &beamLArTPCIntersection)
{
    m_beamLArTPCIntersection = beamLArTPCIntersection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParameters::SetBeamDirection(const pandora::CartesianVector &beamDirection)
{
    m_beamDirection = beamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParameters::SetSelectedFraction(const float selectedFraction)
{
    m_selectedFraction = selectedFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParameters::SetNSelectedHits(const unsigned int nSelectedHits)
{
    m_nSelectedHits = nSelectedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParameters::SetContainmentLimit(const float containmentLimit)
{
    m_containmentLimit = containmentLimit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParameters::GetLArTPCMinX() const
{
    return m_larTPCMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParameters::GetLArTPCMaxX() const
{
    return m_larTPCMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParameters::GetLArTPCMinY() const
{
    return m_larTPCMinY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParameters::GetLArTPCMaxY() const
{
    return m_larTPCMaxY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParameters::GetLArTPCMinZ() const
{
    return m_larTPCMinZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParameters::GetLArTPCMaxZ() const
{
    return m_larTPCMaxZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const BdtBeamParticleIdTool::PlaneVector &BdtBeamParticleIdTool::SliceFeatureParameters::GetPlanes() const
{
    return m_larTPCPlanes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &BdtBeamParticleIdTool::SliceFeatureParameters::GetBeamLArTPCIntersection() const
{
    return m_beamLArTPCIntersection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &BdtBeamParticleIdTool::SliceFeatureParameters::GetBeamDirection() const
{
    return m_beamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParameters::GetSelectedFraction() const
{
    return m_selectedFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int BdtBeamParticleIdTool::SliceFeatureParameters::GetNSelectedHits() const
{
    return m_nSelectedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParameters::GetContainmentLimit() const
{
    return m_containmentLimit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool BdtBeamParticleIdTool::SliceFeatures::IsFeatureVectorAvailable() const
{
    return m_isAvailable;
}

} // namespace lar_content

#endif // #ifndef LAR_BDT_BEAM_PARTICLE_ID_TOOL_H
