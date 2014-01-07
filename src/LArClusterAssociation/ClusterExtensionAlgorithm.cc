/**
 *  @file   LArContent/src/LArClusterAssociation/ClusterExtensionAlgorithm.cc
 * 
 *  @brief  Implementation of the longitudinal extension algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/ClusterExtensionAlgorithm.h"

using namespace pandora;

namespace lar
{

void ClusterExtensionAlgorithm::FillClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterAssociationMatrix clusterAssociationMatrix;
    this->FillClusterAssociationMatrix(clusterVector, clusterAssociationMatrix);
    this->FillClusterMergeMap(clusterAssociationMatrix, clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
