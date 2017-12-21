/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h
 *
 *  @brief  Header file for the transverse extension algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRANSVERSE_EXTENSION_ALGORITHM_H
#define LAR_TRANSVERSE_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TransverseExtensionAlgorithm class
 */
class TransverseExtensionAlgorithm : public ClusterExtensionAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TransverseExtensionAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FillClusterAssociationMatrix(const pandora::ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const;
    void FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Form associations between parent and daughter clusters
     *
     *  @param  parentCluster the parent pointing cluster
     *  @param  pDaughterCluster the address of the daughter cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillClusterAssociationMatrix(const LArPointingCluster &parentCluster, const pandora::Cluster* const pDaughterCluster, ClusterAssociationMatrix &clusterAssociationMatrix) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minClusterLength;              ///<
    float m_maxLongitudinalDisplacement;   ///<
    float m_maxTransverseDisplacement;     ///<
};

} // namespace lar_content

#endif // #ifndef LAR_TRANSVERSE_EXTENSION_ALGORITHM_H
