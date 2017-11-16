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

    typedef std::unordered_map<const pandora::ParticleFlowObject*, const pandora::LArTPC*> PfoToLArTPCMap;

    /**
     *  @brief  StitchingInfo class
     */
    class StitchingInfo
    {
    public:
        PfoToLArTPCMap      m_pfoToLArTPCMap;         ///< Mapping between Pfos and LArTPCs
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

    /**
     *  @brief  Recreate a specified list of pfos in the current pandora instance
     *
     *  @param  inputPfoList the input pfo list
     *  @param  newPfoList to receive the list of new pfos
     */
    pandora::StatusCode Recreate(const pandora::PfoList &inputPfoList, pandora::PfoList &newPfoList) const;

    /**
     *  @brief  Recreate a specified pfo in the current pandora instance
     *
     *  @param  pInputPfo the input pfo
     *  @param  pNewParentPfo the new parent of the new output pfo (nullptr if none)
     *  @param  newPfoList to receive the list of new pfos
     */
    pandora::StatusCode Recreate(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ParticleFlowObject *const pNewParentPfo,
        pandora::PfoList &newPfoList) const;

    /**
     *  @brief  Create a new calo hit in the current pandora instance, based upon the provided input calo hit
     *
     *  @param  pInputCaloHit the address of the input calo hit
     *  @param  pParentCaloHit the address of the parent calo hit
     *
     *  @return the address of the new calo hit
     */
    const pandora::CaloHit *CreateCaloHit(const pandora::CaloHit *const pInputCaloHit, const pandora::CaloHit *const pParentCaloHit) const;

    /**
     *  @brief  Create a new cluster in the current pandora instance, based upon the provided input cluster
     *
     *  @param  pInputCluster the address of the input cluster
     *  @param  newCaloHitList the list of calo hits for the new cluster
     *  @param  newIsolatedCaloHitList the list of isolated calo hits for the new cluster
     *
     *  @return the address of the new cluster
     */
    const pandora::Cluster *CreateCluster(const pandora::Cluster *const pInputCluster, const pandora::CaloHitList &newCaloHitList,
        const pandora::CaloHitList &newIsolatedCaloHitList) const;

    /**
     *  @brief  Create a new vertex in the current pandora instance, based upon the provided input vertex
     *
     *  @param  pInputVertex the address of the input vertex
     *
     *  @return the address of the new vertex
     */
    const pandora::Vertex *CreateVertex(const pandora::Vertex *const pInputVertex) const;

    /**
     *  @brief  Create a new pfo in the current pandora instance, based upon the provided input pfo
     *
     *  @param  pInputPfo the address of the input pfo
     *  @param  newClusterList the list of clusters for the new pfo
     *  @param  newVertexList the list of vertices for the new pfo
     *
     *  @return the address of the new pfo
     */
    const pandora::ParticleFlowObject *CreatePfo(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ClusterList &newClusterList,
        const pandora::VertexList &newVertexList) const;

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

    PandoraInstanceList         m_crWorkerInstances;                ///< The list of cosmic-ray reconstruction worker instances
    const pandora::Pandora     *m_pSlicingWorkerInstance;           ///< The slicing worker instance
    const pandora::Pandora     *m_pSliceNuWorkerInstance;           ///< The per-slice neutrino reconstruction worker instance
    const pandora::Pandora     *m_pSliceCRWorkerInstance;           ///< The per-slice cosmic-ray reconstruction worker instance

    CosmicRayTaggingBaseTool   *m_pCosmicRayTaggingTool;            ///< The address of the cosmic-ray tagging tool
    NeutrinoIdBaseTool         *m_pNeutrinoIdTool;                  ///< The address of the neutrino id tool

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
    typedef std::vector<pandora::PfoList> SliceHypotheses;

    /**
     *  @brief  Select which reconstruction hypotheses to use; neutrino outcomes or cosmic-ray muon outcomes for each slice
     *
     *  @param  nuSliceHypotheses the parent pfos representing the neutrino outcome for each slice
     *  @param  crSliceHypotheses the parent pfos representing the cosmic-ray muon outcome for each slice
     *  @param  sliceNuPfos to receive the list of selected pfos
     */
    virtual void SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses,
        pandora::PfoList &selectedPfos) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_MASTER_ALGORITHM_H
