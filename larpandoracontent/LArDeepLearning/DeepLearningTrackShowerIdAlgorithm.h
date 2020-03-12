/**
 *  @file   larpandoracontent/LArMonitoring/DeepLearningTrackShowerIdAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H
#define LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  DeepLearningTrackShowerIdAlgorithm class
 */
class DeepLearningTrackShowerIdAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DeepLearningTrackShowerIdAlgorithm();
    virtual ~DeepLearningTrackShowerIdAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode Train();
    pandora::StatusCode Infer();

    pandora::StringVector     m_caloHitListNames;    ///< Name of input calo hit list
    std::string               m_modelFileName;       ///< Model file name
    float                     m_xMin;                ///< Min X used in training model
    float                     m_xMax;                ///< Max X used in training model
    float                     m_zMinU;               ///< Min Z used in training model in U view
    float                     m_zMaxU;               ///< Max Z used in training model in U view
    float                     m_zMinV;               ///< Min Z used in training model in V view
    float                     m_zMaxV;               ///< Max Z used in training model in V view
    float                     m_zMinW;               ///< Min Z used in training model in W view
    float                     m_zMaxW;               ///< Max Z used in training model in W view
    int                       m_nBins;               ///< Number of bins used in training model (assumption is same number in X and Z)
    bool                      m_visualize;           ///< Whether to visualize the track shower ID scores
    bool                      m_useTrainingMode;     ///< Training mode
    bool                      m_profile;
    std::string               m_trainingOutputFile;  ///< Output file name for training examples
};

} // namespace lar_content

#endif // LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H
