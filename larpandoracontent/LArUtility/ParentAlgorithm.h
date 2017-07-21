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
