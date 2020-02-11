/**
 *  @file   larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the svm pfo characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SVM_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_SVM_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h"

namespace lar_content
{

/**
 *  @brief  SvmPfoCharacterisationAlgorithm class
 */
class SvmPfoCharacterisationAlgorithm : public PfoCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SvmPfoCharacterisationAlgorithm();

    /**
     *  @brief Destructor
     */
    ~SvmPfoCharacterisationAlgorithm();
class RecoParameters
{
public:
        /**
         *  @brief  Constructor
         */
     RecoParameters();

     unsigned int  m_minPrimaryGoodHits;       ///< the minimum number of primary good Hits
     unsigned int  m_minHitsForGoodView;       ///< the minimum number of Hits for a good view
     unsigned int  m_minPrimaryGoodViews;      ///< the minimum number of primary good views
	 bool          m_foldToPrimaries;          ///< whether to fold all hits to primary pfos and MC particles
     float         m_minHitSharingFraction;    ///< the minimum Hit sharing fraction
    };
protected:

    pandora::StatusCode Run();
    virtual bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    void FillMCToRecoHitsMap(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) const;
    void SelectParticlesByHitCount(const pandora::MCParticleVector &candidateTargets, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap, LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap) const;

    void SelectGoodCaloHits(const pandora::CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap, pandora::CaloHitList &selectedGoodCaloHitList) const;
    void GetMCToSelfMap(const pandora::MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToSelfMap) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    ClusterCharacterisationFeatureTool::FeatureToolVector   m_featureToolVector;         ///< The feature tool map
    PfoCharacterisationFeatureTool::FeatureToolVector       m_featureToolVectorThreeD;   ///< The feature tool map for 3D info
    PfoCharacterisationFeatureTool::FeatureToolVector       m_featureToolVectorNoChargeInfo; ///< The feature tool map for missing W view

    AdaBoostDecisionTree    m_adaBoostDecisionTree;              ///< The Boosted Decision Tree
    AdaBoostDecisionTree    m_adaBoostDecisionTreeNoChargeInfo;  ///< The Boosted Decision Tree for missing W view

    bool                    m_trainingSetMode;              ///< Whether to train
    bool                    m_enableProbability;            ///< Whether to use probabilities instead of binary classification
    bool                    m_useThreeDInformation;         ///< Whether to use 3D information
    float                   m_minProbabilityCut;            ///< The minimum probability to label a cluster as track-like
    unsigned int            m_minCaloHitsCut;               ///< The minimum number of calo hits to qualify as a track

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list

    std::string             m_trainingOutputFile;           ///< The training output file
    std::string             m_filePathEnvironmentVariable;  ///< The environment variable providing a list of paths to svm files
    std::string             m_svmFileName;                  ///< The svm input file
    std::string             m_svmName;                      ///< The name of the svm to find
    std::string             m_svmFileNameNoChargeInfo;      ///< The svm input file for PFOs missing the W view, and thus charge info
    std::string             m_svmNameNoChargeInfo;          ///< The name of the svm to find for PFOs missing the W view, and thus charge info
	bool 					m_writeToTree;///
    RecoParameters m_recoParameters;
};

} // namespace lar_content

#endif // #ifndef LAR_SVM_PFO_CHARACTERISATION_ALGORITHM_H
