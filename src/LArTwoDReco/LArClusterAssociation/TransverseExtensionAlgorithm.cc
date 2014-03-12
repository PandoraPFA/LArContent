/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the transverse extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h"

using namespace pandora;

namespace lar
{


void TransverseExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
   for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        clusterVector.push_back(*iter);

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    // Convert long clusters into pointing clusters
    LArPointingClusterList pointingClusterList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster* const pCluster(*iter);

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength )
            continue;

        try
        {
            pointingClusterList.push_back(LArPointingCluster(*iter));
        }
        catch (StatusCodeException &)
        {
        }
    }

    // Form associations between clusters
    for (LArPointingClusterList::const_iterator iter1 = pointingClusterList.begin(), iterEnd1 = pointingClusterList.end(); iter1 != iterEnd1; ++iter1)
    {
        const LArPointingCluster &parentCluster = *iter1;

        for (ClusterVector::const_iterator iter2 = clusterVector.begin(), iterEnd2 = clusterVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pDaughterCluster(*iter2);

            if (parentCluster.GetCluster() == pDaughterCluster)
                continue;

            this->FillClusterAssociationMatrix(parentCluster, pDaughterCluster, clusterAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::FillClusterAssociationMatrix(const LArPointingCluster &pointingCluster, const Cluster* const pDaughterCluster, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    const Cluster *const pParentCluster(pointingCluster.GetCluster());

    if (pParentCluster == pDaughterCluster)
        return;

    for (unsigned int useInner = 0; useInner < 2; ++useInner)
    {
        const LArPointingCluster::Vertex &pointingVertex(useInner==1 ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());
        const ClusterAssociation::VertexType vertexType(useInner==1 ? ClusterAssociation::INNER : ClusterAssociation::OUTER);

        if (pointingVertex.GetRms() > 0.5f)
            continue;

        const float projectedDisplacement(LArClusterHelper::GetClosestDistance(pointingVertex.GetPosition(), pDaughterCluster));

        if (projectedDisplacement > m_maxLongitudinalDisplacement)
            continue;

        ClusterAssociation::AssociationType associationType(ClusterAssociation::WEAK);
        float figureOfMerit(projectedDisplacement);

        CartesianVector firstCoordinate(0.f, 0.f, 0.f);
        CartesianVector secondCoordinate(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinatesXZ(pDaughterCluster, firstCoordinate, secondCoordinate);

        float firstL(0.f), firstT(0.f), secondT(0.f), secondL(0.f);
        LArPointingClusterHelper::GetImpactParameters(pointingVertex, firstCoordinate, firstL, firstT);
        LArPointingClusterHelper::GetImpactParameters(pointingVertex, secondCoordinate, secondL, secondT);

        const float innerL(firstL < secondL ? firstL : secondL);
        const float innerT(firstL < secondL ? firstT : secondT);
        const float outerL(firstL > secondL ? firstL : secondL);
        const float outerT(firstL > secondL ? firstT : secondT);

        if (innerL > 0.f && innerL < 2.5f && outerL < m_maxLongitudinalDisplacement &&
            innerT < m_maxTransverseDisplacement && outerT < 1.5f * m_maxTransverseDisplacement)
        {
            associationType = ClusterAssociation::STRONG;
            figureOfMerit = outerL;
        }

// --- BEGIN DISPLAY ---
// if(associationType != ClusterAssociation::NONE)
// {
// if (associationType == ClusterAssociation::WEAK) std::cout << " WEAK " << std::endl; else std::cout << " STRONG " << std::endl;
// std::cout << " RMS=" << pointingVertex.GetRms() << std::endl;
// ClusterList tempList1, tempList2;
// tempList1.insert((Cluster*)pParentCluster);
// tempList2.insert((Cluster*)pDaughterCluster);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "ParentCluster", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "DaughterCluster", BLUE);
// PandoraMonitoringApi::ViewEvent();
// }
// --- END DISPLAY ---
        (void) clusterAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pDaughterCluster,
               ClusterAssociation(vertexType, vertexType, associationType, figureOfMerit)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &parentToDaughterMatrix, ClusterMergeMap &clusterMergeMap) const
{
    ClusterAssociationMatrix daughterToParentMatrix;

    // Loop over parent clusters and select nearby daughter clusters that are closer than another parent cluster
    for (ClusterAssociationMatrix::const_iterator iter1 = parentToDaughterMatrix.begin(), iterEnd1 = parentToDaughterMatrix.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster* pParentCluster(iter1->first);
        const ClusterAssociationMap &parentToDaughterMap(iter1->second);

        float maxDisplacementInner(std::numeric_limits<float>::max());
        float maxDisplacementOuter(std::numeric_limits<float>::max());

        // Find the nearest parent cluster
        for (ClusterAssociationMap::const_iterator iter2 = parentToDaughterMap.begin(), iterEnd2 = parentToDaughterMap.end(); iter2 != iterEnd2; ++iter2)
        {
            const ClusterAssociation &clusterAssociation(iter2->second);

            if (clusterAssociation.GetAssociation() == ClusterAssociation::WEAK)
            {
                if (clusterAssociation.GetParent() == ClusterAssociation::INNER && clusterAssociation.GetFigureOfMerit() < maxDisplacementInner)
                    maxDisplacementInner = clusterAssociation.GetFigureOfMerit();

                if (clusterAssociation.GetParent() == ClusterAssociation::OUTER && clusterAssociation.GetFigureOfMerit() < maxDisplacementOuter)
                    maxDisplacementOuter = clusterAssociation.GetFigureOfMerit();
            }
        }

        // Select daughter clusters that are closer than the nearest parent cluster
        for (ClusterAssociationMap::const_iterator iter2 = parentToDaughterMap.begin(), iterEnd2 = parentToDaughterMap.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pDaughterCluster(iter2->first);
            const ClusterAssociation &clusterAssociation(iter2->second);

            if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
            {
                if (clusterAssociation.GetParent() == ClusterAssociation::INNER && clusterAssociation.GetFigureOfMerit() < maxDisplacementInner)
                    (void) daughterToParentMatrix[pDaughterCluster].insert(ClusterAssociationMap::value_type(pParentCluster, clusterAssociation));

                if (clusterAssociation.GetParent() == ClusterAssociation::OUTER && clusterAssociation.GetFigureOfMerit() < maxDisplacementOuter)
                    (void) daughterToParentMatrix[pDaughterCluster].insert(ClusterAssociationMap::value_type(pParentCluster, clusterAssociation));
            }
        }
    }

    // Loop over daughter clusters and select the nearest parent clusters
    for (ClusterAssociationMatrix::const_iterator iter1 = daughterToParentMatrix.begin(), iterEnd1 = daughterToParentMatrix.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster* pDaughterCluster(iter1->first);
        const ClusterAssociationMap &daughterToParentMap(iter1->second);

        Cluster* pParentCluster = NULL;
        float minDisplacement(std::numeric_limits<float>::max());

        for (ClusterAssociationMap::const_iterator iter2 = daughterToParentMap.begin(), iterEnd2 = daughterToParentMap.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pCandidateCluster(iter2->first);
            const ClusterAssociation &clusterAssociation(iter2->second);

            if (clusterAssociation.GetFigureOfMerit() < minDisplacement)
            {
            minDisplacement = clusterAssociation.GetFigureOfMerit();
            pParentCluster = (Cluster*)pCandidateCluster;
            }
        }

        if (pParentCluster)
        {
            clusterMergeMap[pParentCluster].insert((Cluster*)pDaughterCluster);
            clusterMergeMap[pDaughterCluster].insert(pParentCluster);

// --- BEGIN DISPLAY ---
// ClusterList tempList1, tempList2;
// tempList1.insert((Cluster*)pParentCluster);
// tempList2.insert((Cluster*)pDaughterCluster);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "ParentCluster", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "DaughterCluster", BLUE);
// PandoraMonitoringApi::ViewEvent();
// --- END DISPLAY ---
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterLength = 5.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "MinClusterLength", m_minClusterLength));

    m_maxLongitudinalDisplacement = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    m_maxTransverseDisplacement = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
