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

namespace lar_content
{

/**
 *  @brief  MvaPfoCharacterisationAlgorithm class
 */
template<typename T>
class MvaPfoCharacterisationAlgorithm : public PfoCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MvaPfoCharacterisationAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~MvaPfoCharacterisationAlgorithm();

protected:
    virtual bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    ClusterCharacterisationFeatureTool::FeatureToolVector   m_featureToolVector;                ///< The feature tool map
    PfoCharacterisationFeatureTool::FeatureToolVector       m_featureToolVectorThreeD;          ///< The feature tool map for 3D info
    PfoCharacterisationFeatureTool::FeatureToolVector       m_featureToolVectorNoChargeInfo;    ///< The feature tool map for missing W view

    T                       m_mva;                          ///< The mva
    T                       m_mvaNoChargeInfo;              ///< The mva for missing W view

    bool                    m_trainingSetMode;              ///< Whether to train
    bool                    m_enableProbability;            ///< Whether to use probabilities instead of binary classification
    bool                    m_useThreeDInformation;         ///< Whether to use 3D information
    float                   m_minProbabilityCut;            ///< The minimum probability to label a cluster as track-like
    unsigned int            m_minCaloHitsCut;               ///< The minimum number of calo hits to qualify as a track

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list

    std::string             m_trainingOutputFile;           ///< The training output file
    std::string             m_filePathEnvironmentVariable;  ///< The environment variable providing a list of paths to mva files
    std::string             m_mvaFileName;                  ///< The mva input file
    std::string             m_mvaName;                      ///< The name of the mva to find
    std::string             m_mvaFileNameNoChargeInfo;      ///< The mva input file for PFOs missing the W view, and thus charge info
    std::string             m_mvaNameNoChargeInfo;          ///< The name of the mva to find for PFOs missing the W view, and thus charge info
    bool                    m_writeToTree;
    std::string             m_treeName;
    std::string             m_fileName;

    LArMCParticleHelper::PrimaryParameters  m_primaryParameters;    ///< The mc particle primary selection parameters
};

typedef MvaPfoCharacterisationAlgorithm<AdaBoostDecisionTree> BdtPfoCharacterisationAlgorithm;
typedef MvaPfoCharacterisationAlgorithm<SupportVectorMachine> SvmPfoCharacterisationAlgorithm;

} // namespace lar_content

#endif // #ifndef LAR_MVA_PFO_CHARACTERISATION_ALGORITHM_H
