/**
 *  @file   larpandoradlcontent/LArEventClassification/CNNTrackShowerCountingAlgorithm.h
 *
 *  @brief  Header file for the deep learning track/shower counting algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CNN_TRACK_SHOWER_COUNTING_ALGORITHM_H
#define LAR_CNN_TRACK_SHOWER_COUNTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include <torch/torch.h>

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  CNNTrackShowerCountingAlgorithm class
 */
class CNNTrackShowerCountingAlgorithm : public pandora::Algorithm
{
public:
    typedef std::pair<unsigned int, unsigned int> Pixel; // A Pixel is a row, column pair
    typedef std::map<Pixel, float> PixelMap;             // Sparse representation of the pixel map
    typedef std::map<Pixel, std::vector<const pandora::CaloHit *>> PixelToCaloHitsMap;
    typedef std::vector<Pixel> PixelVector;

    /**
     *  @brief A simple class for storing vertex position information
     */
    class VertexPosition
    {
    public:
        /**
        *  @brief  Default constructor
        */
        VertexPosition();

        /**
        *  @brief  Constructor with cartesian vector
        *
        *  @param  pos the vertex position
        */
        VertexPosition(const pandora::CartesianVector &pos);

        /**
        *  @brief  Add information about the vertex position in the pixel map
        *
        *  @param  isDrift to determine if this is the wire or x coordiate
        *  @param  bin the image row or column containing the vertex 
        */
        void AddVertexBin(const bool isDrift, const unsigned int bin);

        /**
        *  @brief  Query the vertex position
        *
        *  @return the vertex position
        */
        pandora::CartesianVector GetPosition() const;

        /**
        *  @brief  Query the vertex position in the drift dimension of the image
        *
        *  @return the vertex drift bin
        */
        unsigned int GetDriftBin() const;

        /**
        *  @brief  Query the vertex position in the wire dimension of the image
        *
        *  @return the vertex wire bin
        */
        unsigned int GetWireBin() const;

    private:
        pandora::CartesianVector m_position;  ///< The vertex position
        unsigned int m_driftBin;              ///< The vertex position in the drift dimension of the image
        unsigned int m_wireBin;               ///< The vertex position in the wire dimension of the image
    };
    typedef std::map<pandora::HitType, VertexPosition> ViewToVertexPositionMap;

    /**
     *  @brief A simple class for storing results from the track and shower counting
     */
    class TrackShowerCountingResults
    {
    public:
        /**
        *  @brief  Add the predicted scores for each network output
        *
        *  @param  view the readout view
        *  @param  nuScores the scores from the neutrino id output
        *  @param  trackScores the scores from the track counting output
        *  @param  showerScores the scores from the shower counting output
        */
        void AddScoresFromView(const pandora::HitType &view, const torch::Tensor &nuScores, const torch::Tensor &trackScores,
            const torch::Tensor &showerScores);

        /**
        *  @brief  Get the predicted scores for neutrino id output
        *
        *  @return the output scores
        */
        const std::vector<float>& GetNuScoresFromView(const pandora::HitType &view) const;

        /**
        *  @brief  Get the neutrino class with the highest score
        *
        *  @return the neutrino class with the highest score
        */
        unsigned int GetNuClassPredictionFromView(const pandora::HitType &view) const;

        /**
        *  @brief  Get the predicted scores for track counting output
        *
        *  @return the output scores
        */
        const std::vector<float>& GetTrackScoresFromView(const pandora::HitType &view) const;

        /**
        *  @brief  Get the track counting class with the highest score
        *
        *  @return the track counting class with the highest score
        */
        unsigned int GetTrackClassPredictionFromView(const pandora::HitType &view) const;

        /**
        *  @brief  Get the predicted scores for shower counting output
        *
        *  @return the output scores
        */
        const std::vector<float>& GetShowerScoresFromView(const pandora::HitType &view) const;

        /**
        *  @brief  Get the shower counting class with the highest score
        *
        *  @return the shower counting class with the highest score
        */
        unsigned int GetShowerClassPredictionFromView(const pandora::HitType &view) const;

    private:
        std::vector<float> TensorToVector(const torch::Tensor &scores) const;

        std::map<pandora::HitType, std::vector<float>> m_nuScores;      ///< Vector of neutrino class scores
        std::map<pandora::HitType, std::vector<float>> m_trackScores;   ///< Vector of track class scores
        std::map<pandora::HitType, std::vector<float>> m_showerScores;  ///< Vector of shower class scores
    };

    /**
     *  @brief  Default constructor
     */
    CNNTrackShowerCountingAlgorithm();

    /**
     *  @brief  Default destructor
     */
    virtual ~CNNTrackShowerCountingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create the input files for training the network
     */
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief  Perform inference to get the event classification
     */
    pandora::StatusCode Infer();

    /**
     *  @brief  Persist the results of the track and shower counting
     *
     *  @param  result the result of the CNN inference
     */
    pandora::StatusCode StorePredictions(const TrackShowerCountingResults &result);

    /**
     *  @brief  Find which particles are visible in the final state
     *
     *  @param  pMCParticleList the pointer to the list of the input MCParticles
     *  @param  pCaloHitList the pointer to the input CaloHitList
     *  @param  mcToHitsMap to show the CaloHits associated to each MCParticle
     */
    void GetVisibleParticles(const pandora::MCParticleList *const pMCParticleList, const pandora::MCParticle *const pMCNeutrino,
        const pandora::CaloHitList *const pCaloHitList, LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    /**
     *  @brief  Create a pixel map from a CaloHitList object
     *
     *  @param  pCaloHitList the pointer to the input CaloHitList
     *  @param  pixelMap the PixelMap to be filled
     *  @param  vertex the reconstructed vertex position
     */
    void CreatePixelMap(const pandora::CaloHitList *const pCaloHitList, PixelMap &pixelMap, VertexPosition &vertex) const;

    /*
     *  @brief  Create input for the network from a calo hit list
     *
     *  @param  pixelMap the PixelMap object from which the input should be made
     *  @param  networkInput The TorchInput object to populate
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeNetworkInputFromPixelMap(const PixelMap &pixelMap, LArDLHelper::TorchInput &networkInput) const;

    /**
     *  @brief  Calculate the best cropped region around the event
     *
     *  @param  pCaloHitList the pointer to the input CaloHitList
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     **/
    void GetHitRegion(const pandora::CaloHitList *const pCaloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const;

    /**
     *  @brief  Find the min and max values for x or z that to create an image of the requested size
     *
     *  @param  pCaloHitList the pointer to the input CaloHitList
     *  @param  vertexPosition the position of the reconstructed interaction vertex
     *  @param  min to take the minimum coordinate value after cropping
     *  @param  max to take the maximum coordinate value after cropping
     *  @param  span the span in coordinate values required
     *  @param  isDrift whether this is the drift coordinate or not
     **/
    void GetCrop1D(const pandora::CaloHitList *const pCaloHitList, VertexPosition &vertexPosition, float &min, float &max,
        const float &span, const bool &isDrift) const;

    /*
     *  @brief  Get a list of all the primary MC particles
     *
     *  @param  pMCParticleList the pointer to the MCParticle list
     *  @param  mcPrimaries to take the list of primary MCParticles
     **/
    void GetMCPrimaries(const pandora::MCParticleList *pMCParticleList, pandora::MCParticleList &mcPrimaries) const;

    /*
     *  @brief  Determine if an MCParticle is tracklike
     *
     *  @param  pMCParticle the pointer to the MCParticle in question
     *
     *  @return if the object is a tracklike particle
     **/
    bool IsTracklike(const pandora::MCParticle *const pMCParticle) const;

    /*
     *  @brief  Determine if an MCParticle is showerlike
     *
     *  @param  pMCParticle the pointer to the MCParticle in question
     *
     *  @return if the object is a showerlike particle
     **/
    bool IsShowerlike(const pandora::MCParticle *const pMCParticle) const;

    /*
     *  @brief  Count the number of primary tracks and showers
     *
     *  @param  mcToHitsMap the MC particles and associated calo hits
     *  @param  nTracks to take the value of the number of primary tracks
     *  @param  nShowers to take the value of the number of primary showers
     **/
    void CountMCPrimaries(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, unsigned int &nTracks, unsigned int &nShowers) const;

    /*
     *  @brief  Get the reconstructed neutrino vertex for each 2D view
     *
     *  @param  vertices to take the vertex position for each view
     **/
    void GetVertexPositions(ViewToVertexPositionMap &vertices) const;

    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingTreeName;           ///< Tree name for training examples
    std::string m_trainingFileName;           ///< Output file name for training examples
    pandora::StringVector m_caloHitListNames; ///< Names of input calo hit lists
    LArDLHelper::TorchModel m_model;          ///< The network model
    unsigned int m_height;                    ///< The height of the images in pixels
    unsigned int m_width;                     ///< The width of the images in pixels
    unsigned int m_wiresPerPixel;             ///< The number of wires per pixel
    float m_driftStep;                        ///< The size of a pixel in the drift direction in cm
    std::string m_vertexListName;             ///< The vertex list name
    unsigned int m_goodMCTrackHits;           ///< The number of hits an MC primary track needs to be considered visible per view
    unsigned int m_goodMCShowerHits;          ///< The number of hits an MC primary shower needs to be considered visible per view
    float m_mcHitWeightThreshold;             ///< Fraction above which a given MCParticle is considered to have been responsible for a hit
    float m_secondaryDistanceThreshold;       ///< Distance from the neutrino vertex below which a secondary is considered primary if the primary is not visible
    unsigned int m_minHits;                   ///< Minimum number of hits to create a training example
    float m_maxChargeThreshold;               ///< Value at which to truncate the charge of each pixel
    std::string m_outputPfoListName;          ///< Name of the output PfoList containing the dummy event pfo
};

} // namespace lar_dl_content

#endif // LAR_CNN_TRACK_SHOWER_COUNTING_ALGORITHM_H
