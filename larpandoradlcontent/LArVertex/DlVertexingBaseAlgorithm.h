/**
 *  @file   larpandoradlcontent/LArVertexing/DlVertexingBaseAlgorithm.h
 *
 *  @brief  Header file for the deep learning vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_VERTEXING_BASE_ALGORITHM_H
#define LAR_DL_VERTEXING_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArObjects/VertexTuple.h"

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningTrackShowerIdAlgorithm class
 */
class DlVertexingBaseAlgorithm : public pandora::Algorithm
{
public:
    typedef std::map<std::pair<int, int>, std::vector<const pandora::CaloHit *>> PixelToCaloHitsMap;

    /**
     *  @brief Default constructor
     */
    DlVertexingBaseAlgorithm();

    ~DlVertexingBaseAlgorithm();

protected:
    typedef std::pair<int, int> Pixel; // A Pixel is a row, column pair
    typedef std::vector<Pixel> PixelVector;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    /**
     *  @brief  Determines the parameters of the canvas for extracting the vertex location.
     *          The network predicts the distance that each pixel associated with a hit is located from the vertex, but says nothing about
     *          the direction. As a result, the ring describing the potential vertices associated with that hit can extend beyond the
     *          original canvas size. This function returns the size of the required canvas and the offset for the bottom left corner.
     *
     *  @param  networkOutput The TorchOutput object populated by the network inference step
     *  @param  pixelVector The vector of populated pixels
     *  @param  columnOffset The output column offset for the canvas
     *  @param  rowOffset The output row offset for the canvas
     *  @param  width The output width for the canvas
     *  @param  height The output height for the canvas
     */
    void GetCanvasParameters(const LArDLHelper::TorchOutput &networkOutput, const PixelVector &pixelVector, int &columnOffset,
        int &rowOffset, int &width, int &height) const;

protected:
    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingOutputFile;         ///< Output file name for training examples
    std::string m_inputVertexListName;        ///< Input vertex list name if 2nd pass
    std::string m_outputVertexListName;       ///< Output vertex list name
    pandora::StringVector m_caloHitListNames; ///< Names of input calo hit lists
    LArDLHelper::TorchModel m_modelU;         ///< The model for the U view
    LArDLHelper::TorchModel m_modelV;         ///< The model for the V view
    LArDLHelper::TorchModel m_modelW;         ///< The model for the W view
    int m_pass;                               ///< The pass of the train/infer step
    int m_nClasses;                           ///< The number of distance classes
    int m_height;                             ///< The height of the images
    int m_width;                              ///< The width of the images
    float m_driftStep;                        ///< The size of a pixel in the drift direction in cm (most relevant in pass 2)
    std::vector<double> m_thresholds;         ///< Distance class thresholds
    std::string m_volumeType;                 ///< The name of the fiducial volume type for the monitoring output
};

} // namespace lar_dl_content

#endif // LAR_DL_VERTEXING_BASE_ALGORITHM_H
