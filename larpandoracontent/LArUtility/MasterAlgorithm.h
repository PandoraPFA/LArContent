/**
 *  @file   larpandoracontent/LArUtility/MasterAlgorithm.h
 *
 *  @brief  Header file for the master algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MASTER_ALGORITHM_H
#define LAR_MASTER_ALGORITHM_H 1

#include "Pandora/AlgorithmTool.h"
#include "Pandora/ExternallyConfiguredAlgorithm.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

namespace lar_content
{

class CosmicRayTaggingBaseTool;
class NeutrinoIdBaseTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MasterAlgorithm class
 */
class MasterAlgorithm : public pandora::ExternallyConfiguredAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MasterAlgorithm();

    /**
     *  @brief  External steering parameters class
     */
    class ExternalSteeringParameters : public pandora::ExternalParameters
    {
    public:
        pandora::InputBool      m_shouldRunAllHitsCosmicReco;       ///< Whether to run all hits cosmic-ray reconstruction
        pandora::InputBool      m_shouldRunStitching;               ///< Whether to stitch cosmic-ray muons crossing between volumes
        pandora::InputBool      m_shouldRunCosmicHitRemoval;        ///< Whether to remove hits from tagged cosmic-rays
        pandora::InputBool      m_shouldRunSlicing;                 ///< Whether to slice events into separate regions for processing
        pandora::InputBool      m_shouldRunNeutrinoRecoOption;      ///< Whether to run neutrino reconstruction for each slice
        pandora::InputBool      m_shouldRunCosmicRecoOption;        ///< Whether to run cosmic-ray reconstruction for each slice
        pandora::InputBool      m_shouldIdentifyNeutrinoSlice;      ///< Whether to identify most appropriate neutrino slice
        pandora::InputBool      m_printOverallRecoStatus;           ///< Whether to print current operation status messages
    };

private:
    pandora::StatusCode Initialize();

    /**
     *  @brief  Create a pandora worker instance to handle a single LArTPC
     *
     *  @param  larTPC the lar tpc
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPC &larTPC, const pandora::DetectorGapList &gapList, const std::string &settingsFile) const;

    /**
     *  @brief  Create a pandora worker instance to handle a number of LArTPCs
     *
     *  @param  larTPCMap the lar tpc map
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPCMap &larTPCMap, const pandora::DetectorGapList &gapList, const std::string &settingsFile) const;

    pandora::StatusCode Run();

    /**
     *  @brief  Copy a specified calo hit to the provided pandora instance
     *
     *  @param  pPandora the address of the target pandora instance
     *  @param  pCaloHit the address of the calo hit
     */
    pandora::StatusCode Copy(const pandora::Pandora *const pPandora, const pandora::CaloHit *const pCaloHit) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Read settings from external steering parameters block, if present, otherwise from xml as standard
     *
     *  @param  pExternalParameters the address of the external parameters (if present)
     *  @param  inputBool the input boolean value, from the external parameters (if present)
     *  @param  xmlHandle the xml handle
     *  @param  xmlTag the xml tag
     *  @param  outputBool to receive the output boolean value
     */
    pandora::StatusCode ReadExternalSettings(const ExternalSteeringParameters *const pExternalParameters, const pandora::InputBool inputBool,
        const pandora::TiXmlHandle xmlHandle, const std::string &xmlTag, bool &outputBool);

    bool                        m_shouldRunAllHitsCosmicReco;       ///< Whether to run all hits cosmic-ray reconstruction
    bool                        m_shouldRunStitching;               ///< Whether to stitch cosmic-ray muons crossing between volumes
    bool                        m_shouldRunCosmicHitRemoval;        ///< Whether to remove hits from tagged cosmic-rays
    bool                        m_shouldRunSlicing;                 ///< Whether to slice events into separate regions for processing
    bool                        m_shouldRunNeutrinoRecoOption;      ///< Whether to run neutrino reconstruction for each slice
    bool                        m_shouldRunCosmicRecoOption;        ///< Whether to run cosmic-ray reconstruction for each slice
    bool                        m_shouldIdentifyNeutrinoSlice;      ///< Whether to identify most appropriate neutrino slice
    bool                        m_printOverallRecoStatus;           ///< Whether to print current operation status messages

    CosmicRayTaggingBaseTool   *m_pCosmicRayTaggingTool;            ///< The address of the cosmic-ray tagging tool
    NeutrinoIdBaseTool         *m_pNeutrinoIdTool;                  ///< The address of the neutrino id tool

    PandoraInstanceList         m_crWorkerInstances;                ///< The list of cosmic-ray reconstruction worker instances
    const pandora::Pandora     *m_pSlicingWorkerInstance;           ///< The slicing worker instance
    const pandora::Pandora     *m_pSliceNuWorkerInstance;           ///< The per-slice neutrino reconstruction worker instance
    const pandora::Pandora     *m_pSliceCRWorkerInstance;           ///< The per-slice cosmic-ray reconstruction worker instance

    std::string                 m_crSettingsFile;                   ///< The cosmic-ray reconstruction settings file
    std::string                 m_nuSettingsFile;                   ///< The neutrino reconstruction settings file
    std::string                 m_slicingSettingsFile;              ///< The slicing settings file
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CosmicRayTaggingBaseTool class
 */
class CosmicRayTaggingBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Find the list of ambiguous pfos (could represent cosmic-ray muons or neutrinos)
     *
     *  @param  parentCosmicRayPfos the list of parent cosmic-ray pfos
     *  @param  ambiguousPfos to receive the list of ambiguous pfos
     */
    virtual void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  NeutrinoIdBaseTool class
 */
class NeutrinoIdBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief SliceProperties class
     */
    class SliceProperties
    {
    public:
        /**
         *  @brief  Default constructor
         */
        SliceProperties();

        float       m_weight;               ///< Generic weight property
        float       m_nuVtxY;               ///< The neutrino vertex x coordinate
        float       m_nuVtxZ;               ///< The neutrino vertex y coordinate
        float       m_nuCosWeightedDir;     ///< Cosine of angle between z-axis and (hit weighted) mean direction of primary particles
    };

    typedef std::map<unsigned int, pandora::PfoList> SliceIndexToPfoListMap;
    typedef std::map<unsigned int, SliceProperties> SliceIndexToPropertiesMap;

    /**
     *  @brief  Fill slice properties given provided neutrino pfos
     *
     *  @param  pPfoList the address of the parent neutrino pfo list
     *  @param  sliceProperties to receive the populated neutrino slice properties
     */
    virtual void FillNeutrinoProperties(const pandora::PfoList *const pPfoList, SliceProperties &sliceProperties) const = 0;

    /**
     *  @brief  Fill slice properties given provided cosmic-ray pfos
     *
     *  @param  pPfoList the address of the parent cosmic-ray pfo list
     *  @param  sliceProperties to receive the populated cosmic-ray slice properties
     */
    virtual void FillCosmicRayProperties(const pandora::PfoList *const pPfoList, SliceProperties &sliceProperties) const = 0;

    /**
     *  @brief  Try to identify a neutrino slice. If such a slice is found, return true and provide the neutrino slice index
     *
     *  @param  sliceIndexToPropertiesMap the slice index to slice properties mapping
     *  @param  neutrinoSliceIndex to receive the neutrino slice index, if a neutrino slice is identified
     *
     *  @return whether a neutrino slice is identified
     */
    virtual bool GetNeutrinoSliceIndex(const SliceIndexToPropertiesMap &sliceIndexToPropertiesMap, unsigned int &neutrinoSliceIndex) const = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline NeutrinoIdBaseTool::SliceProperties::SliceProperties() :
    m_weight(0.f),
    m_nuVtxY(std::numeric_limits<float>::max()),
    m_nuVtxZ(std::numeric_limits<float>::max()),
    m_nuCosWeightedDir(std::numeric_limits<float>::max())
{
}

} // namespace lar_content

#endif // #ifndef LAR_MASTER_ALGORITHM_H
