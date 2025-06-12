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
    typedef std::map<Pixel, float> PixelMap; // Sparse representation of the pixel map
    typedef std::map<Pixel, std::vector<const pandora::CaloHit *>> PixelToCaloHitsMap;
    typedef std::vector<Pixel> PixelVector;
    typedef std::map<pandora::HitType, pandora::CartesianVector> ViewToVertexPositionMap;

    class TrackShowerCountingResults
    {
        public:

        void AddScoresFromView(const pandora::HitType &view, const torch::Tensor &trackScores, const torch::Tensor &showerScores);

        private:

        std::map<pandora::HitType, std::vector<float>> m_trackScores;
        std::map<pandora::HitType, std::vector<float>> m_showerScores;
    };

    /**
     *  @brief Default constructor
     */
    CNNTrackShowerCountingAlgorithm();

    /**
     *  @brief Default destructor
     */
    virtual ~CNNTrackShowerCountingAlgorithm();

private:
  
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Create the input files for training the network
     */ 
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief Perform inference to get the event classification
     */ 
    pandora::StatusCode Infer();

    /**
     *  @brief Create a pixel map from a CaloHitList object
     *
     *  @param  pCaloHitList the pointer to the input CaloHitList
     *  @param  pixelMap the PixelMap to be filled
     */ 
    void CreatePixelMap(const pandora::CaloHitList *const pCaloHitList, PixelMap &pixelMap) const;

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
     *  @brief Calculate the best cropped region around the event
     *
     *  @param  pCaloHitList the pointer to the input CaloHitList
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     **/
    void GetHitRegion(const pandora::CaloHitList *const pCaloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const;

    /**
     *  @brief Find the min and max values for x or z that to create an image of the requested size
     *
     *  @param pCaloHitList the pointer to the input CaloHitList
     *  @param vertexPosition the position of the interaction vertex
     *  @param min to take the minimum coordinate value after cropping
     *  @param max to take the maximum coordinate value after cropping
     *  @param span the span in coordinate values required
     *  @param isDrift whether this is the drift coordinate or not
     **/
    void GetCrop1D(const pandora::CaloHitList *const pCaloHitList, const pandora::CartesianVector &vertexPosition,
        float &min, float &max, const float &span, const bool &isDrift) const;

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
     *  @brief  Get a list of all the primary MC particles
     *
     *  @param  mcHierarchy the hierarchy of MCParticle objects
     *  @param  mcPrimaries to take the lst of primary MCParticles
     **/
    void GetMCPrimaries(const pandora::MCParticleList &mcHierarchy, pandora::MCParticleList &mcPrimaries) const;

    /*
     *  @brief  Count the number of primary tracks and showers
     *
     *  @param  mcToHitsMap the MC particles and associated calo hits
     *  @param  nTracks to take the value of the number of primary tracks
     *  @param  nShowers to take the value of the number of primary showers
     **/
    void CountMCPrimaries(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, unsigned int &nTracks, unsigned int &nShowers) const;

    /*
     *  @brief  Get the reconstructed neutrino vertex (if it exists) for each 2D view
     *
     *  @param  vertices to take the vertex position for each view
     **/
    void GetVertexPositions(ViewToVertexPositionMap &vertices) const;

    bool m_trainingMode;                        ///< Training mode
    bool m_validationMode;                      ///< Gets the truth information to compare to the predictions
    std::string m_trainingTreeName;             ///< Tree name for training examples
    std::string m_trainingFileName;             ///< Output file name for training examples
    pandora::StringVector m_caloHitListNames;   ///< Names of input calo hit lists
    LArDLHelper::TorchModel m_model;            ///< The network model
    unsigned int m_height;                      ///< The height of the images in pixels
    unsigned int m_width;                       ///< The width of the images in pixels
    unsigned int m_wiresPerPixel;               ///< The number of wires per pixel
    float m_driftStep;                          ///< The size of a pixel in the drift direction in cm
    bool m_useVertexForCrops;                   ///< Make use of the reconstructed vertex to perform the cropping
    std::string m_vertexListName;               ///< The vertex list name
    unsigned int m_goodMCPrimaryHits;           ///< The number of hits an MC primary needs to be considered reconstructable per view
    unsigned int m_minHits;                     ///< Minimum number of hits to create a training example
};

} // namespace lar_dl_content

#endif // LAR_CNN_TRACK_SHOWER_COUNTING_ALGORITHM_H
