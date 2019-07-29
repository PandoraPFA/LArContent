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

//    std::string                                 m_caloHitListName;          ///< Name of input calo hit list
    std::string m_modelFileName; ///< Model file name
};

} // namespace lar_content

#endif // LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H
