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
    typedef std::map<const pandora::CaloHit*, std::tuple<int, int, int, int>> CaloHitToPixelMap;
    typedef std::map<int, int> PixelToTileMap;

    pandora::StatusCode Run();
    pandora::StatusCode Train();
    pandora::StatusCode Infer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetHitRegion(const pandora::CaloHitList& caloHitList, float& xMin, float& xMax, float& zMin, float& zMax);
    void GetSparseTileMap(const pandora::CaloHitList& caloHitList, const float xMin, const float zMin, const float tileSize,
        const int nTilesX, PixelToTileMap& sparseMap);

    pandora::StringVector     m_caloHitListNames;    ///< Name of input calo hit list
    std::string               m_modelFileName;       ///< Model file name
    bool                      m_visualize;           ///< Whether to visualize the track shower ID scores
    bool                      m_useTrainingMode;     ///< Training mode
    bool                      m_profile;
    std::string               m_trainingOutputFile;  ///< Output file name for training examples
};

} // namespace lar_content

#endif // LAR_DEEP_LEARNING_TRACK_SHOWER_ID_ALGORITHM_H
