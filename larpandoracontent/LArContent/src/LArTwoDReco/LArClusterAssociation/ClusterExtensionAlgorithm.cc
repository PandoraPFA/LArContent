/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the cluster extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void ClusterExtensionAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterAssociationMatrix clusterAssociationMatrix;
    this->FillClusterAssociationMatrix(clusterVector, clusterAssociationMatrix);
    this->FillClusterMergeMap(clusterAssociationMatrix, clusterMergeMap);
}

} // namespace lar_content
