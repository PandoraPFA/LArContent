/**
 *  @file   larpandoradlcontent/LArMonitoring/DeepLearningTrackShowerIdAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H
#define LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_dl_content
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
    bool                      m_visualize;           ///< Whether to visualize the track shower ID scores
    bool                      m_useTrainingMode;     ///< Training mode
    bool                      m_profile;
    std::string               m_trainingOutputFile;  ///< Output file name for training examples
};

} // namespace lar_dl_content

#endif // LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H
