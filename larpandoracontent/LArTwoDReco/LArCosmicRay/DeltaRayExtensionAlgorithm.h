/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h
 *
 *  @brief  Header file for the delta ray extension algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_EXTENSION_ALGORITHM_H
#define LAR_DELTA_RAY_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayExtensionAlgorithm class
 */
class DeltaRayExtensionAlgorithm : public ClusterExtensionAlgorithm
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

    /**
     *  @brief  Default constructor
     */
    DeltaRayExtensionAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FillClusterAssociationMatrix(const pandora::ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const;
    void FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Form association between two clusters
     *
     *  @param  pParentCluster the parent cluster
     *  @param  pDaughterCluster the daughter cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillClusterAssociationMatrix(const pandora::Cluster *const pParentCluster, const pandora::Cluster *const pDaughterCluster,
        ClusterAssociationMatrix &clusterAssociationMatrix) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float    m_minClusterLength;                ///<
    float    m_maxClusterLength;                ///<

    float    m_maxLongitudinalDisplacement;     ///<
    float    m_maxTransverseDisplacement;       ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRayExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRayExtensionAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_EXTENSION_ALGORITHM_H
