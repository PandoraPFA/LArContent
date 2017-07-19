/**
 *  @file   larpandoracontent/LArUtility/ParentNeutrinoAlgorithm.h
 *
 *  @brief  Header file for the parent neutrino algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PARENT_NEUTRINO_ALGORITHM_H
#define LAR_PARENT_NEUTRINO_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/ParentSlicingBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ParentNeutrinoAlgorithm class
 */
class ParentNeutrinoAlgorithm : public ParentSlicingBaseAlgorithm
{
private:
    pandora::StatusCode Run();

    /**
     *  @brief  Use first-pass 3D event reconstruction to slice events into separate, distinct interactions for processing
     *
     *  @param  sliceList the slice list to receive the slice list
     */
    void PerformSlicing(SliceList &sliceList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                 m_clusteringAlgorithm;              ///< The name of the two dimensional clustering algorithm
    pandora::StringVector       m_twoDAlgorithms;                   ///< The names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_threeDAlgorithms;                 ///< The names of the three dimensional reconstruction algorithms
    pandora::StringVector       m_threeDHitAlgorithms;              ///< The names of the three dimensional hit creation algorithms
    pandora::StringVector       m_vertexAlgorithms;                 ///< The names of the vertex reconstruction algorithms
    pandora::StringVector       m_twoDMopUpAlgorithms;              ///< The names of the two dimensional mop-up algorithms
    pandora::StringVector       m_threeDMopUpAlgorithms;            ///< The names of the three dimensional mop-up algorithms
    pandora::StringVector       m_neutrinoAlgorithms;               ///< The names of the neutrino building algorithms
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_NEUTRINO_ALGORITHM_H
