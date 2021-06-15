/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.h
 *
 *  @brief  Header file for the longitudinal extension algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H
#define LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  LongitudinalExtensionAlgorithm class
 */
class LongitudinalExtensionAlgorithm : public ClusterExtensionAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    LongitudinalExtensionAlgorithm();

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
    void FillClusterAssociationMatrix(
        const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ, ClusterAssociationMatrix &clusterAssociationMatrix) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_clusterMinLength;                    ///<
    float m_clusterMinLayerOccupancy;            ///<
    float m_nodeMaxDisplacement;                 ///<
    float m_nodeMaxCosRelativeAngle;             ///<
    float m_emissionMaxLongitudinalDisplacement; ///<
    float m_emissionMaxTransverseDisplacement;   ///<
    float m_emissionMaxCosRelativeAngle;         ///<
};

} // namespace lar_content

#endif // #ifndef LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H
