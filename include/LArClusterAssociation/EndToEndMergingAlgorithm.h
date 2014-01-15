/**
 *  @file   LArContent/include/LArClusterAssociation/EndToEndMergingAlgorithm.h
 *
 *  @brief  Header file for the cluster extension algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_END_TO_END_MERGING_ALGORITHM_H
#define LAR_END_TO_END_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArClusterAssociation/ClusterExtensionAlgorithm.h"

namespace lar
{

/**
 *  @brief  EndToEndMergingAlgorithm class
 */
class EndToEndMergingAlgorithm : public ClusterExtensionAlgorithm
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



    void FillClusterAssociationMatrix(const pandora::Cluster* const pParentCluster, const pandora::Cluster* const pDaughterCluster,
	ClusterAssociationMatrix &clusterAssociationMatrix) const;



    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    float m_minClusterLength;
    float m_maxLongitudinalDisplacement;
    float m_maxTransverseDisplacement;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EndToEndMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new EndToEndMergingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_END_TO_END_MERGING_ALGORITHM_H
