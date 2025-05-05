/**
 *  @file   larpandoracontent/LArTrackShowerId/MvaPfoCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the mva pfo characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MVA_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_MVA_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"
#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h"

#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
 *  @brief  MvaPfoCharacterisationAlgorithm class
 */
template <typename T>
class MvaPfoCharacterisationAlgorithm : public PfoCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MvaPfoCharacterisationAlgorithm();

    /**
     *  @brief  VertexComparator class for comparison of two points wrt neutrino vertex position
     */
    class VertexComparator
    {
    public:
        /**
         *  @brief  Constructor
         */
        VertexComparator(const pandora::CartesianVector vertexPosition2D);

        /**
         *  @brief  operator <
         *
         *  @param  rhs object for comparison
         *
         *  @return boolean
         */
        bool operator()(const pandora::CaloHit *const left, const pandora::CaloHit *const right) const;

        pandora::CartesianVector m_neutrinoVertex; //The neutrino vertex used to sort
    };

    /**                                                                                                                                                                                                    
     *  @brief  DistanceToVertexComparator class for comparison of two points with different hit type (left and right) wrt neutrino vertex position                                                        
     */
    class DistanceToVertexComparator
    {
    public:
      /**                                                                                                                                                                                                  
       * @brief Constructor                                                                                                                                                                                
       */
      DistanceToVertexComparator(const pandora::CartesianVector vertexPosition2D_A, const pandora::CartesianVector vertexPosition2D_B, 
        const pandora::HitType hitType_A, const pandora::HitType hitType_B);

      /**                                                                                                                                                                                                  
         @brief operator <                                                                                                                                                                                 
         *                                                                                                                                                                                                 
         * @param  rhs object for comparison                                                                                                                                                               
         *                                                                                                                                                                                                 
         * @return boolean                                                                                                                                                                                 
         */
      bool operator()(const pandora::CaloHit *const left, const pandora::CaloHit *const right) const;

      pandora::CartesianVector m_neutrinoVertex_A; // The neutrino vertex used to sort for hit type A, i.e. position is projected onto the correct view based on the hit type                              
      pandora::CartesianVector m_neutrinoVertex_B; // same for hit type B                                                                                                                                  
      pandora::HitType m_hitType_A;
      pandora::HitType m_hitType_B;
    };
    
protected:
    virtual bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    ClusterCharacterisationFeatureTool::FeatureToolMap m_featureToolMap; ///< The feature tool map

    PfoCharacterisationFeatureTool::FeatureToolMap m_featureToolMapThreeD;       ///< FeatureToolMap as a map for 3D info
    PfoCharacterisationFeatureTool::FeatureToolMap m_featureToolMapNoChargeInfo; ///< FeatureToolMap as a map for missing W view

    pandora::StringVector m_algorithmToolNames; ///< Vector of strings saving feature tool order for use in feature calculation
    pandora::StringVector m_algorithmToolNamesNoChargeInfo; ///< Vector of strings saving feature tool order for use in feature calculation (missing W view)

    T m_mva;             ///< The mva
    T m_mvaNoChargeInfo; ///< The mva for missing W view

    bool m_persistFeatures;               ///< Whether to write the features to the properties map
    bool m_trainingSetMode;               ///< Whether to train
    bool m_testBeamMode;                  ///< Whether the training set is from a test beam experiment
    bool m_enableProbability;             ///< Whether to use probabilities instead of binary classification
    bool m_useThreeDInformation;          ///< Whether to use 3D information
    float m_minProbabilityCut;            ///< The minimum probability to label a cluster as track-like
    unsigned int m_minCaloHitsCut;        ///< The minimum number of calo hits to qualify as a track
    bool m_applyFiducialCut;              ///< Whether to apply a fiducial volume cut during training
    float m_fiducialMinX;                 ///< Fiducial volume minimum x
    float m_fiducialMaxX;                 ///< Fiducial volume maximum x
    float m_fiducialMinY;                 ///< Fiducial volume minimum y
    float m_fiducialMaxY;                 ///< Fiducial volume maximum y
    float m_fiducialMinZ;                 ///< Fiducial volume minimum z
    float m_fiducialMaxZ;                 ///< Fiducial volume maximum z
    bool m_applyReconstructabilityChecks; ///< Whether to apply reconstructability checks during training

    std::string m_caloHitListName;    ///< Name of input calo hit list
    std::string m_mcParticleListName; ///< Name of input MC particle list

    std::string m_trainingOutputFile;          ///< The training output file
    std::string m_filePathEnvironmentVariable; ///< The environment variable providing a list of paths to mva files
    std::string m_mvaFileName;                 ///< The mva input file
    std::string m_mvaName;                     ///< The name of the mva to find
    std::string m_mvaFileNameNoChargeInfo;     ///< The mva input file for PFOs missing the W view, and thus charge info
    std::string m_mvaNameNoChargeInfo;         ///< The name of the mva to find for PFOs missing the W view, and thus charge info

    LArMCParticleHelper::PrimaryParameters m_primaryParameters; ///< The mc particle primary selection parameters

private:
    /**
     *  @brief  Checks if the interaction vertex is within the fiducial volume
     *
     *  @param  vertex The coordinates of the vertex
     */
    bool PassesFiducialCut(const pandora::CartesianVector &vertex) const;

    // void OrderCaloHitsByDistanceToVertex(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, 
    //     pandora::CaloHitList &caloHitList);

    void OrderCaloHitsByDistanceToVertex(const pandora::Cluster *const pCluster, pandora::CaloHitList &caloHitList) const;

    void CombineCaloHitListsToHaveCollection(const pandora::CaloHitList &orderedCaloHitList1, const pandora::CaloHitList &orderedCaloHitList2, 
        pandora::CaloHitList &mergedCaloHitList) const;

};

typedef MvaPfoCharacterisationAlgorithm<AdaBoostDecisionTree> BdtPfoCharacterisationAlgorithm;
typedef MvaPfoCharacterisationAlgorithm<SupportVectorMachine> SvmPfoCharacterisationAlgorithm;

} // namespace lar_content

#endif // #ifndef LAR_MVA_PFO_CHARACTERISATION_ALGORITHM_H
