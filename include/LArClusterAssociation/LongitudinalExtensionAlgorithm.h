/**
 *  @file   LArContent/include/LArClusterAssociation/LongitudinalExtensionAlgorithm.h
 * 
 *  @brief  Header file for the cluster extension algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H
#define LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H 1

#include "LArClusterAssociation/ClusterExtensionAlgorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  LongitudinalExtensionAlgorithm class
 */
class LongitudinalExtensionAlgorithm : public ClusterExtensionAlgorithm
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
    void FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ, ClusterAssociationMatrix &clusterAssociationMatrix) const;

    /**
     *  @brief  Remove double-counting from cluster association matrix
     * 
     *  @param  inputAssociationMatrix the input matrix of associations
     *  @param  outputAssociationMatrix the output matrix of reduced associations
     */
    void FillReducedAssociationMatrix(const ClusterAssociationMatrix &inputAssociationMatrix, ClusterAssociationMatrix &outputAssociationMatrix) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    float   m_clusterMinLength;                     ///< 
    float   m_clusterMinLayerOccupancy;             ///<  
    float   m_nodeMaxDisplacement;                  ///< 
    float   m_nodeMaxCosRelativeAngle;              ///< 
    float   m_emissionMaxCosRelativeAngle;          ///<
    float   m_emissionMaxLongitudinalDisplacement;  ///< 
    float   m_emissionMaxTransverseDisplacement;    ///< 
   
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *LongitudinalExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new LongitudinalExtensionAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_LONGITUDINAL_EXTENSION_ALGORITHM_H
