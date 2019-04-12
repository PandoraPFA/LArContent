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

protected:
    virtual bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    ClusterCharacterisationFeatureTool::FeatureToolVector   m_featureToolVector;         ///< The feature tool map
    PfoCharacterisationFeatureTool::FeatureToolVector       m_featureToolVectorThreeD;   ///< The feature tool map for 3D info
    PfoCharacterisationFeatureTool::FeatureToolVector       m_featureToolVectorNoChargeInfo; ///< The feature tool map for missing W view

    SupportVectorMachine    m_supportVectorMachine;              ///< The support vector machine
    SupportVectorMachine    m_supportVectorMachineNoChargeInfo;  ///< The support vector machine for missing W view

    bool                    m_trainingSetMode;              ///< Whether to train
    bool                    m_enableProbability;            ///< Whether to use probabilities instead of binary classification
    bool                    m_useThreeDInformation;         ///< Whether to use 3D information
    float                   m_minProbabilityCut;            ///< The minimum probability to label a cluster as track-like
    unsigned int            m_minCaloHitsCut;               ///< The minimum number of calo hits to qualify as a track

    std::string             m_trainingOutputFile;           ///< The training output file
    std::string             m_filePathEnvironmentVariable;  ///< The environment variable providing a list of paths to svm files
    std::string             m_svmFileName;                  ///< The svm input file
    std::string             m_svmName;                      ///< The name of the svm to find
    std::string             m_svmFileNameNoChargeInfo;      ///< The svm input file for PFOs missing the W view, and thus charge info
    std::string             m_svmNameNoChargeInfo;          ///< The name of the svm to find for PFOs missing the W view, and thus charge info

};

} // namespace lar_content

#endif // #ifndef LAR_SVM_PFO_CHARACTERISATION_ALGORITHM_H
