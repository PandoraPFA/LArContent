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

    /**
     *  @brief  Run two dimensional remnant reconstruction using list names provided via algorithm config
     */
    void TwoDRemnantReconstruction() const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                 m_trackClusteringAlgorithm;         ///< The name of the two dimensional track clustering algorithm
    std::string                 m_deltaRayClusteringAlgorithm;      ///< The name of the two dimensional delta ray clustering algorithm
    std::string                 m_listPruningAlgorithm;             ///< The name of the list pruning algorithm

    pandora::StringVector       m_twoDAlgorithms;                   ///< The names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_threeDAlgorithms;                 ///< The names of the three dimensional reconstruction algorithms
    pandora::StringVector       m_deltaRayAlgorithms;               ///< The names of the delta ray algorithms
    pandora::StringVector       m_twoDRemnantAlgorithms;            ///< The names of the two dimensional remnant algorithms
    pandora::StringVector       m_threeDRemnantAlgorithms;          ///< The names of the three dimensional remnant algorithms
    pandora::StringVector       m_threeDHitAlgorithms;              ///< The names of the three dimensional hit creation algorithms
    pandora::StringVector       m_vertexAlgorithms;                 ///< The names of the vertex reconstruction algorithms
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_COSMIC_RAY_ALGORITHM_H
