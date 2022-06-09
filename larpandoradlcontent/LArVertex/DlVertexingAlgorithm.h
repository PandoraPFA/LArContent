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

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningTrackShowerIdAlgorithm class
 */
class DlVertexingAlgorithm: public pandora::Algorithm
{
public:
    typedef std::map<std::pair<int, int>, std::vector<const pandora::CaloHit*>> PixelToCaloHitsMap;

    /**
     *  @brief Default constructor
     */
    DlVertexingAlgorithm();

    virtual ~DlVertexingAlgorithm();

private:
    class VertexTuple
    {
    public:
        VertexTuple(const pandora::Pandora &pandora, const pandora::CartesianVector &vertexU, const pandora::CartesianVector &vertexV,
            const pandora::CartesianVector &vertexW);

        VertexTuple(const pandora::Pandora &pandora, const pandora::CartesianVector &vertex1, const pandora::CartesianVector &vertex2,
            const pandora::HitType view1, const pandora::HitType view2);

        const pandora::CartesianVector &GetPosition() const;
        float GetChi2() const;
        std::string ToString() const;

    private:
        pandora::CartesianVector m_pos;         ///< Calculated 3D position
        float                   m_chi2;         ///< Chi squared of calculated position
    };

    typedef std::pair<int, int> Pixel;          // A Pixel is a row, column pair
    typedef std::vector<Pixel> PixelVector;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    /*
     *  @brief  Create input for the network from a calo hit list
     *
     *  @param  caloHits The CaloHitList from which the input should be made
     *  @param  networkInput The TorchInput object to populate
     *  @param  pixelVector The output vector of populated pixels
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeNetworkInputFromHits(const pandora::CaloHitList &caloHits, LArDLHelper::TorchInput &networkInput,
        PixelVector &pixelVector) const;

    /*
     *  @brief  Create a list of wire plane-space coordinates from a list of Pixels
     *
     *  @param  caloHits The CaloHitList from which the pixel vector was originally constructed
     *  @param  pixelVector The input vector of populated pixels
     *  @param  positionVector The output vector of wire plane positions
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeWirePlaneCoordinatesFromPixels(const pandora::CaloHitList &caloHits, const PixelVector &pixelVector,
        pandora::CartesianPointVector &positionVector) const;

    /*
     *  @brief  Create a list of wire plane-space coordinates from a canvas
     *
     *  @param  caloHits The CaloHitList from which the pixel vector was originally constructed
     *  @param  canvas The input canvas
     *  @param  canvasWidth The width of the canvas
     *  @param  canvasHeight The height of the canvas
     *  @param  columnOffset The column offset used when populating the canvas
     *  @param  rowOffset The row offset used when populating the canvas
     *  @param  positionVector The output vector of wire plane positions
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeWirePlaneCoordinatesFromCanvas(const pandora::CaloHitList &caloHits, float **canvas, const int canvasWidth,
        const int canvasHeight, const int columnOffset, const int rowOffset, pandora::CartesianPointVector &positionVector) const;

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
    void GetCanvasParameters(const LArDLHelper::TorchOutput &networkOutput, const PixelVector &pixelVector, int &columnOffset, int &rowOffset,
        int &width, int &height) const;

    /**
     *  @brief  Add a filled ring to the specified canvas.
     *          The ring has an inner radius based on the minimum predicted distance to the vertex and an outer radius based on the maximum
     *          predicted distance to the vertex. The centre of the ring is the location of the hit used to predict the distance to the
     *          vertex. Each pixel to be filled is augmented by the specified weight. In this way, once all hits have been considered, a
     *          consensus view emerges of the likely vertex location based on the overlap of various rings centred at different locations.
     *
     *          The underlying implementation is a variant of the Bresenham midpoint circle algorithm and therefore only computes pixel
     *          coordinates for one octant of each circle (one of radius 'inner', one of radius 'outer') and interpolates the fill between
     *          points using integer arithmetic, guaranteeing each pixel of the ring is filled once and only once, and then mirrored to the
     *          remaining seven octants.
     *
     *  @param  networkOutput The TorchOutput object populated by the network inference step
     *  @param  pixelVector The vector of populated pixels
     *  @param  thresholds The fractional distance thresholds representing the classes predicted by the network
     *  @param  columnOffset The output column offset for the canvas
     *  @param  rowOffset The output row offset for the canvas
     *  @param  width The output width for the canvas
     *  @param  height The output height for the canvas
     */
    void DrawRing(float **canvas, const int row, const int col, const int inner, const int outer, const float weight) const;

    /**
     *  @brief  Update the coordinates along the loci of a circle.
     *          When drawing the ring we need an efficient means to determine the next pixel defining the inner and outer loci of the ring.
     *          This update function uses the Bresenham midpoint circle update function to determine this location. The row position is
     *          always incremented by 1 pixel, the column position is left unchanged, or decremented by 1 pixel to best follow the arc of
     *          the true underlying circle.
     *
     *  @param  radius2 The squared radius of the circle under consideration
     *  @param  col The input/output column position to (potentially) update
     *  @param  row The input/output row position to update
     */
    void Update(const int radius, int &col, int &row) const;

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
    pandora::StatusCode CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, pandora::MCParticleList& mcHierarchy)
        const;

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
    void GetHitRegion(const pandora::CaloHitList& caloHitList, float& xMin, float& xMax, float& zMin, float& zMax) const;

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

    /**
     *  @brief  Populate a root true with vertex information.
     */
    void PopulateRootTree(const std::vector<VertexTuple> &vertexTuples, const pandora::CartesianPointVector &vertexCandidatesU,
        const pandora::CartesianPointVector &vertexCandidatesV, const pandora::CartesianPointVector &vertexCandidatesW) const;

    bool                    m_trainingMode;             ///< Training mode
    std::string             m_trainingOutputFile;       ///< Output file name for training examples
    std::string             m_inputVertexListName;      ///< Input vertex list name if 2nd pass
    std::string             m_outputVertexListName;     ///< Output vertex list name
    pandora::StringVector   m_caloHitListNames;         ///< Names of input calo hit lists
    LArDLHelper::TorchModel m_modelU;                   ///< The model for the U view
    LArDLHelper::TorchModel m_modelV;                   ///< The model for the V view
    LArDLHelper::TorchModel m_modelW;                   ///< The model for the W view
    int                     m_event;                    ///< The current event number
    int                     m_pass;                     ///< The pass of the train/infer step
    int                     m_nClasses;                 ///< The number of distance classes
    int                     m_height;                   ///< The height of the images
    int                     m_width;                    ///< The width of the images
    float                   m_maxHitAdc;                ///< Maximum ADC value to allow
    float                   m_regionSize;               ///< The half width/height of the event region to consider in cm
    bool                    m_visualise;                ///< Whether or not to visualise the candidate vertices
    bool                    m_writeTree;                ///< Whether or not to write validation details to a ROOT tree
    std::string             m_rootTreeName;             ///< The ROOT tree name
    std::string             m_rootFileName;             ///< The ROOT file name
    std::mt19937            m_rng;                      ///< The random number generator
    std::vector<double>     m_thresholds;               ///< Distance class thresholds
    const double            PY_EPSILON{1.1920929e-7};   ///< The value of epsilon in Python
};

} // namespace lar_dl_content

#endif // LAR_DL_VERTEXING_ALGORITHM_H

