/**
 *  @file   larpandoradlcontent/LArSignalID/DlSNSignalAlgorithm.h
 *
 *  @brief  Header file for the deep learning signal algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SIGNAL_ALGORITHM_H
#define LAR_DL_SIGNAL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningSignalIdAlgorithm class
 */
class DlSNSignalAlgorithm : public pandora::Algorithm
{
public:
    typedef std::map<std::pair<int, int>, std::vector<const pandora::CaloHit *>> PixelToCaloHitsMap;

    /**
     *  @brief Default constructor
     */
    DlSNSignalAlgorithm();

    virtual ~DlSNSignalAlgorithm();

private:

    typedef std::pair<int, int> Pixel; // A Pixel is a row, column pair
    typedef std::map<const pandora::CaloHit*, Pixel> PixelMap; //A mapping between a CaloHit and a Pixel in the canvas

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();
    pandora::StatusCode CheatedSeparation();
    pandora::StatusCode SignalZoomRandomised();

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
        const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelMap &pixelMap) const;

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

    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingOutputFile;         ///< Output file name for training examples
    std::string m_inputSignalListName;        ///< Input vertex list name if 2nd pass
    pandora::StringVector m_caloHitListNames; ///< Names of input calo hit lists
    LArDLHelper::TorchModel m_modelU;         ///< The model for the U view
    LArDLHelper::TorchModel m_modelV;         ///< The model for the V view
    LArDLHelper::TorchModel m_modelW;         ///< The model for the W view
    int m_event;                              ///< The current event number
    int m_pass;                               ///< The pass of the train/infer step
    int m_height;                             ///< The height of the images
    int m_width;                              ///< The width of the images
    float m_driftStep;                        ///< The size of a pixel in the drift direction in cm (most relevant in pass 2)
    bool m_visualise;                         ///< Whether or not to visualise the candidate vertices
    bool m_visualiseZoom;                     ///< Whether or not to visulaise the zoomed version of detector
    bool m_writeTree;                         ///< Whether or not to write validation details to a ROOT tree
    std::string m_rootTreeName;               ///< The ROOT tree name
    std::string m_rootFileName;               ///< The ROOT file name
    std::mt19937 m_rng;                       ///< The random number generator
    bool m_printOut;                          ///< Whether or not to print out network outputs of CaloHitList names and sizes
    std::string m_signalListNameU;            ///< Output signal CaloHitListU name
    std::string m_signalListNameV;            ///< Output signal CaloHitListV name
    std::string m_signalListNameW;            ///< Output signal CaloHitListW name
    std::string m_signalListName2D;           ///< Output signal CaloHitList2D name
    std::string m_caloHitListName2D;          ///< Input CaloHitList2D name
    std::string m_zoomListNameU;              ///< Output zoom CaloHitListU name
    std::string m_zoomListNameV;              ///< Output zoom CaloHitListV name
    std::string m_zoomListNameW ;             ///< Output zoom CaloHitListW name    
    pandora::StringVector m_inputCaloHitListNames; ///< Names of input calo hit lists, passed from Pass 1 of DLSignalAlg
    std::string m_backgroundListName;         ///< Input Background CaloHitList name
    bool m_applyCheatedSeparation;            ///< Whether cheating to separate background and signal hits 
    bool m_simpleZoom;                        ///< Decide whethere to run a simple loop to find highest adc hit or run network
    long unsigned int m_passOneTrustThreshold;///< Number of pixels in pass one required to trust the wire finding ability, below this                                                           threshold, the algorithm will use highest ADC within Drift Min/Max to set wire limits
    bool m_applySignalZoomRandomised;         ///< Whether to apply a function which reduced combinatorics of background, by reducing size of 							 detctor volume and just containing signal in a smaller, randomised region.
    const int PHOTON_CLASS{2};                ///< Constant for network classification for photons
    const int ELECTRON_CLASS{3};              ///< Constant for network classification for electrons
    const int SIGNAL_CLASS{2};                ///< Constant for network classification for signal
};

} // namespace lar_dl_content

#endif // LAR_DL_SIGNAL_ALGORITHM_H
