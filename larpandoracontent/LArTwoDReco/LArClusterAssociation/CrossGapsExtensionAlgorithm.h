/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsExtensionAlgorithm.h
 *
 *  @brief  Header file for the cross gaps extension algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CROSS_GAPS_EXTENSION_ALGORITHM_H
#define LAR_GROSS_GAPS_EXTENSION_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CrossGapsExtensionAlgorithm class
 */
class CrossGapsExtensionAlgorithm : public ClusterExtensionAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CrossGapsExtensionAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FillClusterAssociationMatrix(const pandora::ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const;
    void FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief Build lists of pointing clusters that are in proximity to a detector gap
     *
     *  @param clusterVector the input vector of clusters
     *  @param pointingClusterList the list pointing clusters that are in proximity to a detector gap
     */
    void BuildPointingClusterList(const pandora::ClusterVector &clusterVector, LArPointingClusterList &pointingClusterList) const;

    /**
     *  @brief Use pointing information to determine whether two clusters are associated
     *
     *  @param pointingVertex1 the first pointing vertex
     *  @param pointingVertex2 the second pointing vertex
     */
    bool IsAssociated(const LArPointingCluster::Vertex &pointingVertex1, const LArPointingCluster::Vertex &pointingVertex2) const;

    /**
     *  @brief Determine whether a start and end position sit either side of a gap
     *
     *  @param minZ the start position
     *  @param maxZ the end position
     *  @param hitType the hitType
     */
    bool IsAcrossGap(const float minZ, const float maxZ, const pandora::HitType hitType) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float   m_minClusterLength;               ///<
    float   m_maxGapTolerance;                ///<
    float   m_maxTransverseDisplacement;      ///<
    float   m_maxRelativeAngle;               ///<
    float   m_minCosRelativeAngle;            ///<
};

} // namespace lar_content

#endif // #ifndef LAR_CROSS_GAPS_EXTENSION_ALGORITHM_H
