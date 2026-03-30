/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlTrackCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the deep learning track characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_TRACK_CHARACTERISATION_ALGORITHM_H
#define LAR_DL_TRACK_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_dl_content
{

/**
 *  @brief  DlTrackCharacterisationAlgorithm class
 */
class DlTrackCharacterisationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlTrackCharacterisationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~DlTrackCharacterisationAlgorithm();

    typedef std::vector<std::vector<float>> TrackHitFeatures;
    typedef std::map<pandora::HitType, TrackHitFeatures> ViewToTrackHitFeaturesMap;

    enum class AuxInputs
    {
        kNTrkChildren = 0,
        kNShwChildren,
        kNDescendants,
        kNDescendantHits,
        kRecoTier,
        kNHitsU,
        kNHitsV,
        kNHitsW
    };

    /**
    *  @brief  Simple class to store track features
    */
    class TrackFeatures
    {
    public:
        /**
         *  @brief  Add information from a CaloHit
         *
         *  @param  pCaloHit the hit to extract features from
         *  @param  view the readout view that the hit is from
         *  @param  chargeThreshold to truncate the hit charge
         */
        void AddHit(const pandora::CaloHit *const pCaloHit, const pandora::HitType view, const float chargeThreshold);

        /**
         *  @brief  Add an auxillary feature
         *
         *  @param  value the value of the feature to add
         */
        void AddAuxillaryValue(const float value);

        /**
         *  @brief  Get the map of view to track hit features
         *
         *  @return the map of view to track hit features
         */
        const ViewToTrackHitFeaturesMap &GetViewToTrackHitFeaturesMap() const;

        /**
         *  @brief  Get the features for all hits of a given type
         *
         *  @return the features of all hits in the requested view
         */
        const TrackHitFeatures &GetAllHitFeatures(const pandora::HitType view) const;

        /**
         *  @brief  Get the features for a given hit
         *
         *  @param  h the index of the hit
         *
         *  @return the features of all hits
         */
        const std::vector<float> &GetHit(const unsigned int h, const pandora::HitType view) const;

        /**
         *  @brief  Get the number of hits
         *
         *  @return the number of hits
         */
        unsigned int GetNHits(const pandora::HitType view) const;

        /**
         *  @brief  Get the auxillary features
         *
         *  @return the auxillary features
         */
        const std::vector<float> &GetAuxillaryFeatures() const;

        /**
         *  @brief  Get the number of auxillary features
         *
         *  @return the number of auxillary features
         */
        unsigned int GetNAuxillaryFeatures() const;

        /**
         *  @brief  Set the truth target values for training
         *
         *  @param  pdg the PDG code of the track
         */
        void SetTruePdg(const int pdg);

        /**
         *  @brief  Get the true PDG code
         *
         *  @return the PDG code
         */
        int GetTruePdg() const;

        /**
         *  @brief  Set whether or not the track exits the detector
         *
         *  @param  isExiting if the track exits or not
         */
        void SetExitingStatus(bool isExiting);

        /**
         *  @brief  Get whether or not the track exits the detector
         *
         *  @return if the track exits or not
         */
        bool GetExitingStatus() const;

        /**
         *  @brief  Add an empty TrackFeatures object to view without hits
         */
        void AddFeaturesToMissingViews();

    private:
        ViewToTrackHitFeaturesMap m_viewToTrackHitFeatures; ///< Map of view to track hit features
        std::vector<float> m_auxillaryFeatures;             ///< Vector of auxillary track features
        bool m_isExiting;                                   ///< Whether the track exits the detector or not
        int m_truePdg;                                      ///< True pdg code for the track
    };

    typedef std::map<const pandora::ParticleFlowObject *, TrackFeatures> PfoToTrackFeaturesMap;

protected:
    pandora::StatusCode Run();

    /**
     *  @brief  Create the training sample of tracks
     *
     *  @return status code
     */
    pandora::StatusCode PrepareTrainingSample() const;

    /**
     *  @brief  Perform the inference to get the predicted track properties
     *
     *  @return status code
     */
    pandora::StatusCode Infer();

    /**
     *  @brief  Store the inference results in the pfo metadata
     *
     *  @param  pPfo the pfo for which the results were obtained
     *  @param  scores the tensor containing the raw results (before softmax)
     */
    void SaveResultsToPfoMetadata(const pandora::ParticleFlowObject *pPfo, const LArDLHelper::TorchOutput &scores);

    /**
     *  @brief  Get the track features for all track-like pfos
     *
     *  @param  pfoToTrackFeaturesMap to store the track features for the pfos
     *  @return status code
     */
    pandora::StatusCode GetAllTrackFeatures(PfoToTrackFeaturesMap &pfoToTrackFeaturesMap) const;

    /**
     *  @brief  Get the track features for a given pfo
     *
     *  @param  caloHits the list of 3D calo hits associated to the pfo
     *
     *  @return the track features for this pfo
     */
    TrackFeatures GetTrackFeatures(const pandora::CaloHitList &caloHits) const;

    /**
     *  @brief  Get the track auxillary feature variables
     *
     *  @param  pPfo the pointer to the pfo we want the track features for
     *  @param  trackFeatures the object in which to store the auxillary variables
     *
     */
    void GetTrackAuxillaryInfo(const pandora::ParticleFlowObject *pPfo, TrackFeatures &trackFeatures) const;

    /**
     *  @brief  Check if the track leaves the detector
     *
     *  @param  caloHits the list of 3D calo hits associated to the pfo
     *
     *  @return if the track exits the detector
     */
    bool IsTrackExiting(const pandora::CaloHitList &caloHits) const;

    /**
     *  @brief  Create the sequence input tensor for the network {x, y, z, q}
     *
     *  @param  trackHitFeatures the features to extract to the tensor
     *  @param  sequenceInput the torch tensor to take the sequence input
     */
    void CreateSequenceInput(const TrackHitFeatures &trackHitFeatures, LArDLHelper::TorchInput &sequenceInput) const;

    /**
     *  @brief  Create the auxillary input tensor for the network
     *
     *  @param  trackFeatures the features to extract to the tensor
     *  @param  auxillaryInput the torch tensor to take the auxillary input
     */
    void CreateAuxillaryInput(const TrackFeatures &trackFeatures, LArDLHelper::TorchInput &auxillaryInput) const;

    /**
     *  @brief  Calculate the detector limits
     */
    void GetDetectorLimits();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_trainingMode; ///< Flag to set whether we are in training or inference modes

    unsigned int m_minTrackHits;    ///< Minimum number of hits for a track to be considered
    unsigned int m_sequenceLength;  ///< Maximum number of hits considered in the sequence input
    float m_maxChargeThreshold;     ///< Threshold to which hit charges are set if they exceed the threshold
    float m_detEdgeTolerance;       ///< Tolerance on determining if a track exits the detector
    std::string m_trackPfoListName; ///< Name of the track-like pfo list

    std::string m_trainingFileName; ///< Name of the training file to create
    std::string m_trainingTreeName; ///< Name of the ROOT tree within the training file

    std::string m_networkFileName;            ///< Name of the network model
    LArDLHelper::TorchModel m_exitingModel;   ///< The network model used for inference of exiting tracks
    LArDLHelper::TorchModel m_containedModel; ///< The network model used for inference of contained tracks

    float m_detLimitMinX; ///< Minimum detector x value
    float m_detLimitMaxX; ///< Maximum detector x value
    float m_detLimitMinY; ///< Minimum detector y value
    float m_detLimitMaxY; ///< Maximum detector y value
    float m_detLimitMinZ; ///< Minimum detector z value
    float m_detLimitMaxZ; ///< Maximum detector z value
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_TRACK_CHARACTERISATION_ALGORITHM_H
