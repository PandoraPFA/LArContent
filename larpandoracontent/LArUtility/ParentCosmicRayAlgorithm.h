/**
 *  @file   larpandoracontent/LArUtility/ParentCosmicRayAlgorithm.h
 *
 *  @brief  Header file for the parent cosmic ray algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PARENT_COSMIC_RAY_ALGORITHM_H
#define LAR_PARENT_COSMIC_RAY_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/ParentSlicingBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ParentCosmicRayAlgorithm class
 */
class ParentCosmicRayAlgorithm : public ParentSlicingBaseAlgorithm
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

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                 m_trackClusteringAlgorithm;         ///< CR: the name of the two dimensional track clustering algorithm
    std::string                 m_deltaRayClusteringAlgorithm;      ///< CR: the name of the two dimensional delta ray clustering algorithm
    pandora::StringVector       m_twoDAlgorithms;                   ///< CR: the names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_threeDAlgorithms;                 ///< CR: the names of the three dimensional reconstruction algorithms
    pandora::StringVector       m_deltaRayAlgorithms;               ///< CR: the names of the delta ray algorithms
    pandora::StringVector       m_twoDRemnantAlgorithms;            ///< CR: the names of the two dimensional remnant algorithms
    pandora::StringVector       m_threeDRemnantAlgorithms;          ///< CR: the names of the three dimensional remnant algorithms
    pandora::StringVector       m_threeDHitAlgorithms;              ///< CR: the names of the three dimensional hit creation algorithms
    pandora::StringVector       m_vertexAlgorithms;                 ///< CR: the names of the vertex reconstruction algorithms
    std::string                 m_listPruningAlgorithm;             ///< CR: the name of the list pruning algorithm
    std::string                 m_listMovingAlgorithm;              ///< CR: the name of the list moving algorithm

    pandora::StringVector       m_nuTwoDAlgorithms;                 ///< Nu: the names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_nuThreeDAlgorithms;               ///< Nu: the names of the three dimensional reconstruction algorithms
    pandora::StringVector       m_nuThreeDHitAlgorithms;            ///< Nu: the names of the three dimensional hit creation algorithms
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_COSMIC_RAY_ALGORITHM_H
