/**
 *  @file   LArContent/include/LArClusterAssociation/TransverseExtensionAlgorithm.h
 * 
 *  @brief  Header file for the hello world algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRANSVERSE_EXTENSION_ALGORITHM_H
#define LAR_TRANSVERSE_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArClusterAssociation/ClusterExtensionAlgorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  TransverseExtensionAlgorithm class
 */
class TransverseExtensionAlgorithm : public ClusterExtensionAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FillClusterAssociationMatrix(const pandora::ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const;
    void FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const;



    /**
     *  @brief  Form association between two pointing clusters
     * 
     *  @param  parentCluster the parent pointing cluster
     *  @param  daughterCluster the daughter pointing cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillAssociationMatrix(const LArPointingCluster &parentCluster, const LArPointingCluster &daughterCluster, ClusterAssociationMatrix &clusterAssociationMatrix) const;


    /**
     *  @brief  Form association between a parent pointing cluster and a daughter cluster
     * 
     *  @param  parentCluster the parent pointing cluster
     *  @param  pDaughterCluster the daughter cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillAssociationMatrix(const LArPointingCluster &parentCluster, const pandora::Cluster* const pDaughterCluster, ClusterAssociationMatrix &clusterAssociationMatrix) const;
    


    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    float m_maxClusterLength;

    float m_maxLongitudinalDisplacement;
    float m_maxTransverseDisplacement;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TransverseExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new TransverseExtensionAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRANSVERSE_EXTENSION_ALGORITHM_H
