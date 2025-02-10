/**
 *  @file   larpandoradlcontent/LArVertexing/DlVertexingAlgorithm.h
 *
 *  @brief  Header file for the deep learning vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_VERTEXING_ALGORITHM_H
#define LAR_DL_VERTEXING_ALGORITHM_H 1

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
class DlVertexingAlgorithm : public pandora::Algorithm
{
public:
    typedef std::map<std::pair<int, int>, std::vector<const pandora::CaloHit *>> PixelToCaloHitsMap;

    /**
     *  @brief Default constructor
     */
    DlVertexingAlgorithm();

    virtual ~DlVertexingAlgorithm();

private:
    typedef std::pair<int, int> Pixel; // A Pixel is a row, column pair
    typedef std::vector<Pixel> PixelVector;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    /*
     *  @brief  Create input for the network from a calo hit list
     *
     *  @param  caloHits The CaloHitList from which the input should be made
     *  @param  view The wire plane view
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     *  @param  networkInput The TorchInput object to populate
     *  @param  pixelVector The output vector of populated pixels
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeNetworkInputFromHits(const pandora::CaloHitList &caloHits, const pandora::HitType view, const float xMin,
        const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector) const;

    /*
     *  @brief  Create a list of wire plane-space coordinates from a canvas
     *
     *  @param  canvas The input canvas
     *  @param  canvasWidth The width of the canvas
     *  @param  canvasHeight The height of the canvas
     *  @param  columnOffset The column offset used when populating the canvas
     *  @param  rowOffset The row offset used when populating the canvas
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     *  @param  positionVector The output vector of wire plane positions
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeWirePlaneCoordinatesFromCanvas(float **canvas, const int canvasWidth, const int canvasHeight,
        const int columnOffset, const int rowOffset, const pandora::HitType view, const float xMin, const float xMax, const float zMin,
        const float zMax, pandora::CartesianPointVector &positionVector) const;

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

    /*
     *  @brief  Retrieve the map from MC to calo hits for reconstructable particles
     *
     *  @param  mcToHitsMap The map to populate
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    /*
     *  @brief  Construct a list of the MC particles from the MC to calo hits map, completing the interaction hierarchy with the invisible
     *          upstream particles.
     *
     *  @param  mcToHitsMap The map of reconstructible MC particles to calo hits
     *  @param  mcHierarchy The output list of MC particles representing the interaction
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, pandora::MCParticleList &mcHierarchy) const;

    /*
     *  @brief  Determine the physical bounds associated with a CaloHitList.
     *
     *  @param  caloHitList The CaloHitList for which bounds should be determined
     *  @param  xMin The output minimum x value
     *  @param  xMax The output maximum x value
     *  @param  zMin The output minimum z value
     *  @param  zMax The output maximum z value
     *
     *  @return The StatusCode resulting from the function
     */
    void GetHitRegion(const pandora::CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const;

    /**
     *  @brief Create a vertex list from the candidate vertices.
     *
     *  @param  candidates The candidate positions with which to create the list.
     *
     *  @return The StatusCode resulting from the function
     */
    pandora::StatusCode MakeCandidateVertexList(const pandora::CartesianPointVector &positions);

    /**
     *  @brief  Retrieve the true neutrino vertex position.
     */
    void GetTrueVertexPosition(float &x, float &y, float &z) const;

    /**
     *  @brief  Retrieve the true neutrino vertex position.
     */
    void GetTrueVertexPosition(float &x, float &u, float &v, float &w) const;

    /**
     *  @brief  Retrieve the true neutrino vertex.
     */
    const pandora::CartesianVector &GetTrueVertex() const;

#ifdef MONITORING
    /**
     *  @brief  Populate a root true with vertex information.
     */
    void PopulateRootTree(const std::vector<VertexTuple> &vertexTuples, const pandora::CartesianPointVector &vertexCandidatesU,
        const pandora::CartesianPointVector &vertexCandidatesV, const pandora::CartesianPointVector &vertexCandidatesW) const;
#endif

    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingOutputFile;         ///< Output file name for training examples
    std::string m_inputVertexListName;        ///< Input vertex list name if 2nd pass
    std::string m_outputVertexListName;       ///< Output vertex list name
    pandora::StringVector m_caloHitListNames; ///< Names of input calo hit lists
    LArDLHelper::TorchModel m_modelU;         ///< The model for the U view
    LArDLHelper::TorchModel m_modelV;         ///< The model for the V view
    LArDLHelper::TorchModel m_modelW;         ///< The model for the W view
    int m_event;                              ///< The current event number
    int m_pass;                               ///< The pass of the train/infer step
    int m_nClasses;                           ///< The number of distance classes
    int m_height;                             ///< The height of the images
    int m_width;                              ///< The width of the images
    float m_driftStep;                        ///< The size of a pixel in the drift direction in cm (most relevant in pass 2)
    bool m_visualise;                         ///< Whether or not to visualise the candidate vertices
    bool m_writeTree;                         ///< Whether or not to write validation details to a ROOT tree
    std::string m_rootTreeName;               ///< The ROOT tree name
    std::string m_rootFileName;               ///< The ROOT file name
    std::mt19937 m_rng;                       ///< The random number generator
    std::vector<double> m_thresholds;         ///< Distance class thresholds
    std::string m_volumeType;                 ///< The name of the fiducial volume type for the monitoring output
};

} // namespace lar_dl_content

#endif // LAR_DL_VERTEXING_ALGORITHM_H
