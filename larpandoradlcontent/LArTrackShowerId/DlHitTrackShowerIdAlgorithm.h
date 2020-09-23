/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
#define LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_dl_content
{

/**
 *  @brief  DlHitTrackShowerIdAlgorithm class
 */
class DlHitTrackShowerIdAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlHitTrackShowerIdAlgorithm();

    virtual ~DlHitTrackShowerIdAlgorithm();

private:
    typedef std::map<const pandora::CaloHit*, std::tuple<int, int, int, int>> CaloHitToPixelMap;
    typedef std::map<int, int> PixelToTileMap;

    pandora::StatusCode Run();

    /**
     *  @brief  Produce files that act as inputs to network training
     */
    pandora::StatusCode Train();

    /**
     *  @brief  Run network inference
     */
    pandora::StatusCode Infer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Identify the XZ range containing the hits for an event
     *
     *  @param  caloHitList The list of CaloHits for which the range is to be found
     *  @param  xMin The output minimum x-coordinate
     *  @param  xMax The output maximum x-coordinate
     *  @param  zMin The output minimum z-coordinate
     *  @param  zMax The output maximum z-coordinate
     */
    void GetHitRegion(const pandora::CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax);

    /**
     *  @brief  Populate a map between pixels and tiles
     *
     *  @param  caloHitList The list of CaloHits for which the map is to be populated
     *  @param  xMin The minimum x-coordinate
     *  @param  zMin The minimum z-coordinate
     *  @param  nTilesX The number of tiles in the x direction
     *  @param  sparseMap The output map between pixels and tiles
     */
    void GetSparseTileMap(const pandora::CaloHitList &caloHitList, const float xMin, const float zMin, const int nTilesX, PixelToTileMap &sparseMap);

    pandora::StringVector     m_caloHitListNames;    ///< Name of input calo hit list
    std::string               m_modelFileNameU;      ///< Model file name for U view
    std::string               m_modelFileNameV;      ///< Model file name for V view
    std::string               m_modelFileNameW;      ///< Model file name for W view
    LArDLHelper::TorchModel   m_modelU;              ///< Model for the U view
    LArDLHelper::TorchModel   m_modelV;              ///< Model for the V view
    LArDLHelper::TorchModel   m_modelW;              ///< Model for the W view
    int                       m_imageHeight;         ///< Height of images in pixels
    int                       m_imageWidth;          ///< Width of images in pixels
    float                     m_tileSize;            ///< Size of tile in cm
    bool                      m_visualize;           ///< Whether to visualize the track shower ID scores
    bool                      m_useTrainingMode;     ///< Training mode
    std::string               m_trainingOutputFile;  ///< Output file name for training examples
};

} // namespace lar_dl_content

#endif // LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
