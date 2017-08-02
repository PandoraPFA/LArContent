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

/**
 *  @brief  ParentAlgorithm class
 */
class ParentAlgorithm : public ParentSlicingBaseAlgorithm
{
private:
    pandora::StatusCode Run();

    /**
     *  @brief SliceProperties class
     */
    class SliceProperties
    {
        // Placeholder
    };

    typedef std::map<unsigned int, pandora::PfoList> SliceIndexToPfoListMap;
    typedef std::map<unsigned int, SliceProperties> SliceIndexToPropertiesMap;

    /**
     *  @brief  Run the cosmic-ray reconstruction on the list of all input hits
     *
     *  @param  parentCosmicRayPfos to receive the list of parent cosmic-ray pfos
     */
    void RunAllHitsCosmicRayReconstruction(pandora::PfoList &parentCosmicRayPfos) const;

    /**
     *  @brief  Find the list of ambiguous pfos (could represent cosmic-ray muons or neutrinos)
     *
     *  @param  parentCosmicRayPfos the list of parent cosmic-ray pfos
     *  @param  ambiguousPfos to receive the list of ambiguous pfos
     */
    void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos) const;

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
     *  @brief  Try to identify a neutrino slice. If such a slice is found, return true and provide the neutrino slice index
     *
     *  @param  sliceIndexToPropertiesMap the slice index to slice properties mapping
     *  @param  neutrinoSliceIndex to receive the neutrino slice index, if a neutrino slice is identified
     *
     *  @return whether a neutrino slice is identified
     */
    bool GetNeutrinoSliceIndex(const SliceIndexToPropertiesMap &sliceIndexToPropertiesMap, unsigned int &neutrinoSliceIndex) const;

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

    std::string                 m_crParentListName;                   ///< Output: the name of the cr parent pfo list
    std::string                 m_crDaughterListName;                 ///< Output: the name of the cr daughter pfo list 
    std::string                 m_outputListPrefix;                   ///< Output: the prefix applied to output list names
    std::string                 m_outputListPruningAlgorithm;         ///< Output: the name of the list pruning algorithm    

    std::string                 m_crTrackClusteringAlgorithm;         ///< CR: the name of the two dimensional track clustering algorithm
    std::string                 m_crDeltaRayClusteringAlgorithm;      ///< CR: the name of the two dimensional delta ray clustering algorithm
    pandora::StringVector       m_crTwoDAlgorithms;                   ///< CR: the names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_crThreeDAlgorithms;                 ///< CR: the names of the three dimensional reconstruction algorithms
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
    pandora::StringVector       m_nuThreeDAlgorithms;                 ///< Nu: the names of the three dimensional reconstruction algorithms
    pandora::StringVector       m_nuThreeDHitAlgorithms;              ///< Nu: the names of the three dimensional hit creation algorithms
    pandora::StringVector       m_nuVertexAlgorithms;                 ///< Nu: the names of the vertex reconstruction algorithms
    pandora::StringVector       m_nuTwoDMopUpAlgorithms;              ///< Nu: the names of the two dimensional mop-up algorithms
    pandora::StringVector       m_nuThreeDMopUpAlgorithms;            ///< Nu: the names of the three dimensional mop-up algorithms
    pandora::StringVector       m_nuNeutrinoAlgorithms;               ///< Nu: the names of the neutrino building algorithms
    std::string                 m_nuListDeletionAlgorithm;            ///< Nu: the name of the list deletion algorithm
    std::string                 m_nuListMovingAlgorithm;              ///< Nu: the name of the list moving algorithm
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_ALGORITHM_H
