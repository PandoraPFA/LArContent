/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the bounded cluster merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode BoundedClusterMergingAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BoundedClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    m_vertexVetoRadius = 2.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexVetoRadius", m_vertexVetoRadius));

    m_minCaloHitsInCluster = 4;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinCaloHitsInCluster", m_minCaloHitsInCluster));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
