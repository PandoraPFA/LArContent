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
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    ClusterCharacterisationFeatureTool::FeatureToolVector   m_featureToolVector;    ///< The feature tool map
    SupportVectorMachine    m_supportVectorMachine;         ///< The support vector machine

    bool                    m_trainingSetMode;              ///< Whether to train
    bool                    m_ratioVariables;               ///< Whether to divide all variables by the straight line length
    bool                    m_enableProbability;            ///< Whether to use probabilities instead of binary classification
    unsigned int            m_minCaloHitsCut;               ///< The minimum number of calo hits to qualify as a track

    std::string             m_trainingOutputFile;           ///< The training output file
    std::string             m_svmFileName;                  ///< The svm input file
    std::string             m_svmName;                      ///< The name of the svm to find
};

} // namespace lar_content

#endif // #ifndef LAR_SVM_PFO_CHARACTERISATION_ALGORITHM_H
