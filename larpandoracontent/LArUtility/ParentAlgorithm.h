/**
 *  @file   larpandoracontent/LArUtility/ParentAlgorithm.h
 *
 *  @brief  Header file for the parent algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PARENT_ALGORITHM_H
#define LAR_PARENT_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/ParentSlicingBaseAlgorithm.h"

namespace lar_content
{

class CosmicRayTaggingBaseTool;
class NeutrinoIdBaseTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ParentAlgorithm class
 */
class ParentAlgorithm : public ParentSlicingBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ParentAlgorithm();

    /**
     *  @brief  External steering parameters class
     */
    class ExternalSteeringParameters : public ExternalParameters
    {
    public:
        pandora::InputBool      m_shouldRunAllHitsCosmicReco;       ///< Whether to run all hits cosmic-ray reconstruction
        pandora::InputBool      m_shouldRunCosmicHitRemoval;        ///< Whether to remove hits from tagged cosmic-rays
        pandora::InputBool      m_shouldRunSlicing;                 ///< Whether to slice events into separate regions for processing
        pandora::InputBool      m_shouldRunNeutrinoRecoOption;      ///< Whether to run neutrino reconstruction for each slice
        pandora::InputBool      m_shouldRunCosmicRecoOption;        ///< Whether to run cosmic-ray reconstruction for each slice
        pandora::InputBool      m_shouldIdentifyNeutrinoSlice;      ///< Whether to identify most appropriate neutrino slice
        pandora::InputBool      m_printOverallRecoStatus;           ///< Whether to print current operation status messages
    };

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
    
        float       m_weight;           ///< Generic weight property
    };

    typedef std::map<unsigned int, pandora::PfoList> SliceIndexToPfoListMap;
    typedef std::map<unsigned int, SliceProperties> SliceIndexToPropertiesMap;

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Run the cosmic-ray reconstruction on the list of all input hits
     *
     *  @param  parentCosmicRayPfos to receive the list of parent cosmic-ray pfos
     */
    void RunAllHitsCosmicRayReconstruction(pandora::PfoList &parentCosmicRayPfos) const;

    /**
     *  @brief  Remove ambiguous pfos from the cosmic-ray working lists
     *
     *  @param  ambiguousPfos the list of ambiguous pfos
     */
    void RemoveAmbiguousCosmicRayPfos(const pandora::PfoList &ambiguousPfos) const;

    /**
     *  @brief  Run the fast reconstruction and slicing (if configured to do so)
     *
     *  @param  sliceList to receive the slice list
     */
    void RunFastReconstructionAndSlicing(SliceList &sliceList) const;

    /**
     *  @brief  Run neutrino and cosmic-ray reconstructions for each slice, leaving the cosmic-ray outcomes in the output lists.
     *          Maps of slice index to lists of parent cosmic-ray particles, and slice index to slice properties, are populated.
     *
     *  @param  sliceList the slice list
     *  @param  sliceToCosmicRayPfosMap the slice index to cosmic-ray parent pfo list mapping
     *  @param  sliceIndexToPropertiesMap the slice index to slice properties mapping
     */
    void ReconstructSlices(const SliceList &sliceList, SliceIndexToPfoListMap &sliceToCosmicRayPfosMap, SliceIndexToPropertiesMap &sliceIndexToPropertiesMap) const;

    /**
     *  @brief  Remove the cosmic-ray outcome from the output lists for a given slice index
     *
     *  @param  sliceToCosmicRayPfosMap the slice index to cosmic-ray parent pfo list mapping
     *  @param  neutrinoSliceIndex the relevant slice index
     */
    void RemoveSliceCosmicRayReconstruction(const SliceIndexToPfoListMap &sliceToCosmicRayPfosMap, const unsigned int neutrinoSliceIndex) const;

    /**
     *  @brief  Rerun the neutrino reconstruction for a given slice index and add the results to the output lists
     *
     *  @param  sliceList the slice list
     *  @param  neutrinoSliceIndex the relevant slice index
     */
    void AddSliceNeutrinoReconstruction(const SliceList &sliceList, const unsigned int neutrinoSliceIndex) const;

    void FastReconstruction() const;

    /**
     *  @brief  Perform cosmic-ray reconstruction using the provided slice and its index
     *
     *  @param  slice the slice
     *  @param  sliceIndexString the slice index string/identifier
     */
    void CosmicRayReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const;

    /**
     *  @brief  Perform neutrino reconstruction using the provided slice and its index
     *
     *  @param  slice the slice
     *  @param  sliceIndexString the slice index string/identifier
     */
    void NeutrinoReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const;

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

    bool                        m_shouldRunAllHitsCosmicReco;         ///< Steering: whether to run all hits cosmic-ray reconstruction
    bool                        m_shouldRunCosmicHitRemoval;          ///< Steering: whether to remove hits from tagged cosmic-rays
    bool                        m_shouldRunNeutrinoRecoOption;        ///< Steering: whether to run neutrino reconstruction for each slice
    bool                        m_shouldRunCosmicRecoOption;          ///< Steering: whether to run cosmic-ray reconstruction for each slice
    bool                        m_shouldIdentifyNeutrinoSlice;        ///< Steering: whether to identify most appropriate neutrino slice
    bool                        m_printOverallRecoStatus;             ///< Steering: whether to print current operation status messages

    CosmicRayTaggingBaseTool   *m_pCosmicRayTaggingTool;              ///< The address of the cosmic-ray tagging tool
    NeutrinoIdBaseTool         *m_pNeutrinoIdTool;                    ///< The address of the neutrino id tool

    std::string                 m_crParentListName;                   ///< Output: the name of the cr parent pfo list
    std::string                 m_crDaughterListName;                 ///< Output: the name of the cr daughter pfo list 
    std::string                 m_nuParentListName;                   ///< Output: the name of the nu parent pfo list
    std::string                 m_outputListPrefix;                   ///< Output: the prefix applied to output list names
    std::string                 m_outputListPruningAlgorithm;         ///< Output: the name of the list pruning algorithm    

    std::string                 m_crTrackClusteringAlgorithm;         ///< CR: the name of the two dimensional track clustering algorithm
    std::string                 m_crDeltaRayClusteringAlgorithm;      ///< CR: the name of the two dimensional delta ray clustering algorithm
    pandora::StringVector       m_crTwoDAlgorithms;                   ///< CR: the names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_crThreeDTrackAlgorithms;            ///< CR: the names of the three dimensional track reconstruction algorithms
    pandora::StringVector       m_crDeltaRayAlgorithms;               ///< CR: the names of the delta ray algorithms
    pandora::StringVector       m_crTwoDRemnantAlgorithms;            ///< CR: the names of the two dimensional remnant algorithms
    pandora::StringVector       m_crThreeDRemnantAlgorithms;          ///< CR: the names of the three dimensional remnant algorithms
    pandora::StringVector       m_crThreeDHitAlgorithms;              ///< CR: the names of the three dimensional hit creation algorithms
    pandora::StringVector       m_crVertexAlgorithms;                 ///< CR: the names of the vertex reconstruction algorithms
    std::string                 m_crListPruningAlgorithm;             ///< CR: the name of the list pruning algorithm
    std::string                 m_crListDeletionAlgorithm;            ///< CR: the name of the list deletion algorithm
    std::string                 m_crListMovingAlgorithm;              ///< CR: the name of the list moving algorithm

    std::string                 m_nuClusteringAlgorithm;              ///< Nu: the name of the two dimensional clustering algorithm
    pandora::StringVector       m_nuTwoDAlgorithms;                   ///< Nu: the names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_nuThreeDTrackAlgorithms;            ///< Nu: the names of the three dimensional track reconstruction algorithms
    pandora::StringVector       m_nuThreeDShowerAlgorithms;           ///< Nu: the names of the three dimensional shower reconstruction algorithms
    pandora::StringVector       m_nuThreeDRecoveryAlgorithms;         ///< Nu: the names of the three dimensional recovery reconstruction algorithms
    pandora::StringVector       m_nuThreeDHitAlgorithms;              ///< Nu: the names of the three dimensional hit creation algorithms
    pandora::StringVector       m_nuVertexAlgorithms;                 ///< Nu: the names of the vertex reconstruction algorithms
    pandora::StringVector       m_nuTwoDMopUpAlgorithms;              ///< Nu: the names of the two dimensional mop-up algorithms
    pandora::StringVector       m_nuThreeDMopUpAlgorithms;            ///< Nu: the names of the three dimensional mop-up algorithms
    pandora::StringVector       m_nuNeutrinoAlgorithms;               ///< Nu: the names of the neutrino building algorithms
    std::string                 m_nuListDeletionAlgorithm;            ///< Nu: the name of the list deletion algorithm
    std::string                 m_nuListMovingAlgorithm;              ///< Nu: the name of the list moving algorithm

    bool                        m_nuRepeatThreeDTrackReco;            ///< Nu: whether to repeat 3d track algs in between 3d shower and recovery algs
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline ParentAlgorithm::SliceProperties::SliceProperties() :
    m_weight(0.f)
{
}

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
     *  @brief  Fill slice properties given provided neutrino pfos
     *
     *  @param  pPfoList the address of the parent neutrino pfo list
     *  @param  sliceProperties to receive the populated neutrino slice properties
     */
    virtual void FillNeutrinoProperties(const pandora::PfoList *const pPfoList, ParentAlgorithm::SliceProperties &sliceProperties) const = 0;

    /**
     *  @brief  Fill slice properties given provided cosmic-ray pfos
     *
     *  @param  pPfoList the address of the parent cosmic-ray pfo list
     *  @param  sliceProperties to receive the populated cosmic-ray slice properties
     */
    virtual void FillCosmicRayProperties(const pandora::PfoList *const pPfoList, ParentAlgorithm::SliceProperties &sliceProperties) const = 0;

    /**
     *  @brief  Try to identify a neutrino slice. If such a slice is found, return true and provide the neutrino slice index
     *
     *  @param  sliceIndexToPropertiesMap the slice index to slice properties mapping
     *  @param  neutrinoSliceIndex to receive the neutrino slice index, if a neutrino slice is identified
     *
     *  @return whether a neutrino slice is identified
     */
    virtual bool GetNeutrinoSliceIndex(const ParentAlgorithm::SliceIndexToPropertiesMap &sliceIndexToPropertiesMap, unsigned int &neutrinoSliceIndex) const = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_ALGORITHM_H
