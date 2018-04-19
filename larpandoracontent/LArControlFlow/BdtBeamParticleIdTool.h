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

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"

namespace lar_content
{

class AdaBoostDecisionTree;

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
     *  @param  rhs the BdtBeamParticleIdTool to copy
     */
    BdtBeamParticleIdTool(const BdtBeamParticleIdTool &rhs);

    /**
     *  @brief  Assignment operator
     *
     *  @param  rhs the BdtBeamParticleIdTool to assign
     */
    BdtBeamParticleIdTool &operator=(const BdtBeamParticleIdTool &rhs);

    /**
     *  @brief  Destructor
     */
    ~BdtBeamParticleIdTool();

    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &beamSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

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
         *  @param  a0 point on the line
         *  @param  a vector pointing along the line
         */
        pandora::CartesianVector GetLineIntersection(const pandora::CartesianVector &a0, const pandora::CartesianVector &a) const;

    private:
        pandora::CartesianVector      m_unitNormal;                         ///< Unit normal to plane
        pandora::CartesianVector      m_point;                              ///< A point on the plane
        float                         m_d;                                  ///< Parameter defining a plane
    };

    typedef std::vector<Plane> PlaneVector;

    /**
     *  @brief  SliceFeatureParamters class
     */ 
    class SliceFeatureParamters
    {
    public:
        /**
         *  @brief  Constructor
         */
        SliceFeatureParamters();

        /**
         *  @brief  Set TPC geometry information
         *
         *  @param  tpcMinX the global TPC volume minimum x extent
         *  @param  tpcMaxX the global TPC volume maximum x extent
         *  @param  tpcMinY the global TPC volume minimum y extent
         *  @param  tpcMaxY the global TPC volume maximum y extent
         *  @param  tpcMinZ the global TPC volume minimum z extent
         *  @param  tpcMaxZ the global TPC volume maximum z extent
         */
       void SetTPCGeometryInformation(float &tpcMinX, float &tpcMaxX, float &tpcMinY, float &tpcMaxY, float &tpcMinZ, float &tpcMaxZ);

       /**
        *  @brief  Set m_beamTPCIntersection
        *
        *  @param  beamTPCIntersection intersection of beam and TPC
        */
       void SetBeamTPCIntersection(pandora::CartesianVector &beamTPCIntersection);

       /**
        *  @brief  Set m_beamDirection
        *
        *  @param  beamDirection direction of beam
        */
       void SetBeamDirection(pandora::CartesianVector &beamDirection);

       /**
        *  @brief  Set m_selectedFraction
        *
        *  @param  selectedFraction fraction of hits to use in 3D cluster fits
        */
       void SetSelectedFraction(float &selectedFraction);

       /**
        *  @brief  Set m_nSelectedHits
        *
        *  @param  nSelectedHits minimum number of hits to use in 3D cluster fits
        */
       void SetNSelectedHits(unsigned int &nSelectedHits);

       /**
        *  @brief  Get m_tpcMinX
        */
       float GetTPCMinX() const;

       /**
        *  @brief  Get m_tpcMaxX
        */
       float GetTPCMaxX() const;

       /**
        *  @brief  Get m_tpcMinY
        */
       float GetTPCMinY() const;

       /**
        *  @brief  Get m_tpcMaxY
        */
       float GetTPCMaxY() const;

       /**
        *  @brief  Get m_tpcMinZ
        */
       float GetTPCMinZ() const;

       /**
        *  @brief  Get m_tpcMaxZ
        */
       float GetTPCMaxZ() const;

       /**
        *  @brief  Get vector of planes
        */
       PlaneVector GetPlanes() const;

       /**
        *  @brief  Get the beam TPC intersection 
        */
       pandora::CartesianVector GetBeamTPCIntersection() const;

       /**
        *  @brief  Get the beam direction
        */
       pandora::CartesianVector GetBeamDirection() const;

       /**
        *  @brief  Get m_selectedFraction
        */
       float GetSelectedFraction() const;

       /**
        *  @brief  Set m_nSelectedHits
        *
        *  @param  nSelectedHits minimum number of hits to use in 3D cluster fits
        */
       unsigned int GetNSelectedHits() const;

    private:
        float                       m_tpcMinX;                ///< Global TPC volume minimum x extent 
        float                       m_tpcMaxX;                ///< Global TPC volume maximum x extent
        float                       m_tpcMinY;                ///< Global TPC volume minimum y extent
        float                       m_tpcMaxY;                ///< Global TPC volume maximum y extent
        float                       m_tpcMinZ;                ///< Global TPC volume minimum z extent
        float                       m_tpcMaxZ;                ///< Global TPC volume maximum z extent
        PlaneVector                 m_tpcPlanes;              ///< Vector of all planes making up global TPC volume
        pandora::CartesianVector    m_beamTPCIntersection;    ///< Intersection of beam and global TPC volume
        pandora::CartesianVector    m_beamDirection;          ///< Beam direction
        float                       m_selectedFraction;       ///< Fraction of hits to use in 3D cluster fits
        unsigned int                m_nSelectedHits;          ///< Minimum number of hits to use in 3D cluster fits
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
         *  @param  pTool address of the tool using this class
         *  @param  geometryInfo geometry information block
         */
         SliceFeatures(const pandora::PfoList &nuPfos, const pandora::PfoList &crPfos, const BdtBeamParticleIdTool *const pTool, const SliceFeatureParamters *pSliceFeatureParamters);

        /**
         *  @brief  Copy constructor
         *
         *  @param  rhs the SliceFeatures to copy
         */
        SliceFeatures(const SliceFeatures &rhs);

        /**
         *  @brief  Assignment operator
         *
         *  @param  rhs the SliceFeatures to assign
         */
        SliceFeatures &operator=(const SliceFeatures &rhs);

        /**
         *  @brief  Destructor
         */
        ~SliceFeatures();

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
        void GetFeatureVector(LArMvaHelper::MvaFeatureVector &featureVector) const;
        
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
        void GetSelectedCaloHits(const pandora::CaloHitList &inputCaloHitList, pandora::CaloHitList &outputCaloHitList,
            double &closestHitToFaceDistance) const;

        /** 
         *  @brief  Find the intercepts of a line with the protoDUNE detector
         *
         *  @param  a0 a point on the line in question
         *  @param  majorAxis the direction of the line in question
         *  @param  interceptOne to receive the first intersection between line and protoDUNE detector 
         *  @param  interceptTwo to receive the second intersection between line and protoDUNE detector 
         */
        void GetTPCIntercepts(const pandora::CartesianVector &a0, const pandora::CartesianVector &majorAxis,
            pandora::CartesianVector &interceptOne, pandora::CartesianVector &interceptTwo) const;

        /**
         *  @brief  Check if a given 3D spacepoint is inside the global TPC volume
         *
         *  @param  spacePoint
         */
        bool IsContained(const pandora::CartesianVector &spacePoint) const;

        bool                            m_isAvailable;               ///< Is the feature vector available
        const SliceFeatureParamters    *m_pSliceFeatureParamters;    ///< Geometry information block
        LArMvaHelper::MvaFeatureVector  m_featureVector;             ///< The MVA feature vector
        const BdtBeamParticleIdTool    *m_pTool;                     ///< The tool that owns this
    };

    typedef std::vector<SliceFeatures> SliceFeaturesVector;
    typedef std::pair<unsigned int, float> UintFloatPair;
    typedef std::map<const pandora::MCParticle*, int> MCParticleToIntMap;

    /**
     *  @brief  Initialize TPC geometry 
     */
    pandora::StatusCode Initialize();

    /**
     *  @brief  Get the features of each slice
     *
     *  @param  pTool the address of the this NeutrinoId tool
     *  @param  nuSliceHypotheses the input neutrino slice hypotheses
     *  @param  crSliceHypotheses the input cosmic slice hypotheses
     *  @param  sliceFeaturesVector vector to hold the slice features
     */
    void GetSliceFeatures(const BdtBeamParticleIdTool *const pTool, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, SliceFeaturesVector &sliceFeaturesVector) const;

    /**
     *  @brief  Select all pfos under the same hypothesis
     *
     *  @param  hypotheses the lists of slices under a certain hypothesis
     *  @param  selectedPfos the list of pfos to populate
     */
    void SelectAllPfos(const SliceHypotheses &hypotheses, pandora::PfoList &selectedPfos) const;

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
     *  @param  bestSliceIndicies vector of slice indicies passing quality cuts 
     */
    void GetBestMCSliceIndicies(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::IntVector &bestSliceIndicies) const;

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
     *  @param  hitList output list of all 2d hits in the input pfos
     *  @param  reconstructableHits to check if part of before adding
     */
    void Collect2DHits(const pandora::PfoList &pfos, pandora::CaloHitList &hitList, const pandora::CaloHitList &reconstructableHits) const;

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
     *  @param  nuSliceHypotheses the input neutrino slice hypotheses
     *  @param  crSliceHypotheses the input cosmic slice hypotheses
     *  @param  sliceFeaturesVector vector holding the slice features
     *  @param  selectedPfos the list of pfos to populate
     */
    void SelectPfosByAdaBDTScore(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, const SliceFeaturesVector &sliceFeaturesVector, pandora::PfoList &selectedPfos) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Training
    bool                            m_useTrainingMode;                      ///< Should use training mode. If true, training examples will be written to the output file
    std::string                     m_trainingOutputFile;                   ///< Output file name for training examples
    float                           m_minPurity;                            ///< Minimum purity of the best slice to use event for training
    float                           m_minCompleteness;                      ///< Minimum completeness of the best slice to use event for training

    // Classification
    AdaBoostDecisionTree            m_adaBoostDecisionTree;                 ///< The adaptive boost decision tree
    std::string                     m_filePathEnvironmentVariable;          ///< The environment variable providing a list of paths to bdt files
    unsigned int                    m_maxNeutrinos;                         ///< The maximum number of neutrinos to select in any one event 
    float                           m_minAdaBDTScore;                       ///< Minimum score required to classify a slice as a beam particle
    SliceFeatureParamters          *m_pSliceFeatureParamters;               ///< Geometry information block
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParamters::SetBeamTPCIntersection(pandora::CartesianVector &beamTPCIntersection)
{
    m_beamTPCIntersection = beamTPCIntersection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParamters::SetBeamDirection(pandora::CartesianVector &beamDirection)
{
    m_beamDirection = beamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParamters::SetSelectedFraction(float &selectedFraction)
{
    m_selectedFraction = selectedFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void BdtBeamParticleIdTool::SliceFeatureParamters::SetNSelectedHits(unsigned int &nSelectedHits)
{
    m_nSelectedHits = nSelectedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParamters::GetTPCMinX() const 
{
    return m_tpcMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParamters::GetTPCMaxX() const 
{
    return m_tpcMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParamters::GetTPCMinY() const 
{
    return m_tpcMinY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParamters::GetTPCMaxY() const 
{
    return m_tpcMaxY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParamters::GetTPCMinZ() const 
{
    return m_tpcMinZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParamters::GetTPCMaxZ() const 
{
    return m_tpcMaxZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline BdtBeamParticleIdTool::PlaneVector BdtBeamParticleIdTool::SliceFeatureParamters::GetPlanes() const
{
    return m_tpcPlanes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector BdtBeamParticleIdTool::SliceFeatureParamters::GetBeamTPCIntersection() const
{
    return m_beamTPCIntersection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector BdtBeamParticleIdTool::SliceFeatureParamters::GetBeamDirection() const
{
    return m_beamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BdtBeamParticleIdTool::SliceFeatureParamters::GetSelectedFraction() const
{
    return m_selectedFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int BdtBeamParticleIdTool::SliceFeatureParamters::GetNSelectedHits() const
{
    return m_nSelectedHits;
}

} // namespace lar_content

#endif // #ifndef LAR_BDT_BEAM_PARTICLE_ID_TOOL_H
