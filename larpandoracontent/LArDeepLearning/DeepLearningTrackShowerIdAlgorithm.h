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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector     m_caloHitListNames;    ///< Name of input calo hit list
    std::string               m_modelFileName;       ///< Model file name
    float                     m_xMin;                ///< Min X used in training model
    float                     m_xMax;                ///< Max X used in training model
    float                     m_zMin;                ///< Min Z used in training model
    float                     m_zMax;                ///< Max Z used in training model
    int                       m_nBins;               ///< Number of bins used in training model (assumption is same number in X and Z)
    bool                      m_visualize;           ///< Whether to visualize the track shower ID scores
};

} // namespace lar_content

#endif // LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H
