/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h
 *
 *  @brief  Header file for the cosmic-ray extension algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_EXTENSION_ALGORITHM_H
#define LAR_COSMIC_RAY_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayExtensionAlgorithm class
 */
class CosmicRayExtensionAlgorithm : public ClusterExtensionAlgorithm
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
     *  @param  clusterI the first pointing cluster
     *  @param  clusterJ the second pointing cluster
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillClusterAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ,
        ClusterAssociationMatrix &clusterAssociationMatrix) const;

    /**
     *  @brief  Calculate RMS deviation of a cluster with respect to the reference line
     *
     *  @param  pCluster the input cluster
     *  @param  position the intercept of the reference line
     *  @param  direction the direction of the reference line
     */
    float CalculateRms(const pandora::Cluster *const pCluster, const pandora::CartesianVector &position,
        const  pandora::CartesianVector &direction) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float    m_minClusterLength;                ///<
    float    m_minSeedClusterLength;            ///<
    float    m_maxTransverseDisplacement;       ///<
    float    m_maxLongitudinalDisplacement;     ///<
    float    m_minCosRelativeAngle;             ///<
    float    m_maxAverageRms;                   ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayExtensionAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_EXTENSION_ALGORITHM_H
