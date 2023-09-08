/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the transverse extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TransverseExtensionAlgorithm::TransverseExtensionAlgorithm() :
    m_minClusterLength(5.f),
    m_maxLongitudinalDisplacement(10.f),
    m_maxTransverseDisplacement(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

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
        const Cluster *const pCluster(*iter);

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
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
            const Cluster *const pDaughterCluster(*iter2);

            if (parentCluster.GetCluster() == pDaughterCluster)
                continue;

            this->FillClusterAssociationMatrix(parentCluster, pDaughterCluster, clusterAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::FillClusterAssociationMatrix(const LArPointingCluster &pointingCluster,
    const Cluster *const pDaughterCluster, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    const Cluster *const pParentCluster(pointingCluster.GetCluster());

    if (pParentCluster == pDaughterCluster)
        return;

    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), LArClusterHelper::GetClusterHitType(pParentCluster))};
    const float maxLongitudinalDisplacementAdjusted{ratio * m_maxLongitudinalDisplacement};
    const float maxTransverseDisplacementAdjusted{ratio * m_maxTransverseDisplacement};

    for (unsigned int useInner = 0; useInner < 2; ++useInner)
    {
        const LArPointingCluster::Vertex &pointingVertex(useInner == 1 ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());
        const ClusterAssociation::VertexType vertexType(useInner == 1 ? ClusterAssociation::INNER : ClusterAssociation::OUTER);

        if (pointingVertex.GetRms() > 0.5f)
            continue;

        const float projectedDisplacement(LArClusterHelper::GetClosestDistance(pointingVertex.GetPosition(), pDaughterCluster));

        if (projectedDisplacement > maxLongitudinalDisplacementAdjusted)
            continue;

        ClusterAssociation::AssociationType associationType(ClusterAssociation::WEAK);
        float figureOfMerit(projectedDisplacement);

        CartesianVector firstCoordinate(0.f, 0.f, 0.f);
        CartesianVector secondCoordinate(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pDaughterCluster, firstCoordinate, secondCoordinate);

        float firstL(0.f), firstT(0.f), secondT(0.f), secondL(0.f);
        LArPointingClusterHelper::GetImpactParameters(pointingVertex, firstCoordinate, firstL, firstT);
        LArPointingClusterHelper::GetImpactParameters(pointingVertex, secondCoordinate, secondL, secondT);

        const float innerL(firstL < secondL ? firstL : secondL);
        const float innerT(firstL < secondL ? firstT : secondT);
        const float outerL(firstL > secondL ? firstL : secondL);
        const float outerT(firstL > secondL ? firstT : secondT);

        if (innerL > 0.f && innerL < 2.5f && outerL < maxLongitudinalDisplacementAdjusted && innerT < maxTransverseDisplacementAdjusted &&
            outerT < 1.5f * maxTransverseDisplacementAdjusted)
        {
            associationType = ClusterAssociation::STRONG;
            figureOfMerit = outerL;
        }

        (void)clusterAssociationMatrix[pParentCluster].insert(
            ClusterAssociationMap::value_type(pDaughterCluster, ClusterAssociation(vertexType, vertexType, associationType, figureOfMerit)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &parentToDaughterMatrix, ClusterMergeMap &clusterMergeMap) const
{
    ClusterAssociationMatrix daughterToParentMatrix;

    // Loop over parent clusters and select nearby daughter clusters that are closer than another parent cluster
    ClusterVector sortedParentClusters;
    for (const auto &mapEntry : parentToDaughterMatrix)
        sortedParentClusters.push_back(mapEntry.first);
    std::sort(sortedParentClusters.begin(), sortedParentClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : sortedParentClusters)
    {
        const ClusterAssociationMap &daughterToAssociationMap(parentToDaughterMatrix.at(pParentCluster));

        float maxDisplacementInner(std::numeric_limits<float>::max());
        float maxDisplacementOuter(std::numeric_limits<float>::max());

        // Find the nearest parent cluster
        ClusterVector sortedLocalDaughterClusters;
        for (const auto &mapEntry : daughterToAssociationMap)
            sortedLocalDaughterClusters.push_back(mapEntry.first);
        std::sort(sortedLocalDaughterClusters.begin(), sortedLocalDaughterClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedLocalDaughterClusters)
        {
            const ClusterAssociation &clusterAssociation(daughterToAssociationMap.at(pDaughterCluster));

            if (clusterAssociation.GetAssociation() == ClusterAssociation::WEAK)
            {
                if (clusterAssociation.GetParent() == ClusterAssociation::INNER && clusterAssociation.GetFigureOfMerit() < maxDisplacementInner)
                    maxDisplacementInner = clusterAssociation.GetFigureOfMerit();

                if (clusterAssociation.GetParent() == ClusterAssociation::OUTER && clusterAssociation.GetFigureOfMerit() < maxDisplacementOuter)
                    maxDisplacementOuter = clusterAssociation.GetFigureOfMerit();
            }
        }

        // Select daughter clusters that are closer than the nearest parent cluster
        for (const Cluster *const pDaughterCluster : sortedLocalDaughterClusters)
        {
            const ClusterAssociation &clusterAssociation(daughterToAssociationMap.at(pDaughterCluster));

            if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
            {
                if (clusterAssociation.GetParent() == ClusterAssociation::INNER && clusterAssociation.GetFigureOfMerit() < maxDisplacementInner)
                    (void)daughterToParentMatrix[pDaughterCluster].insert(ClusterAssociationMap::value_type(pParentCluster, clusterAssociation));

                if (clusterAssociation.GetParent() == ClusterAssociation::OUTER && clusterAssociation.GetFigureOfMerit() < maxDisplacementOuter)
                    (void)daughterToParentMatrix[pDaughterCluster].insert(ClusterAssociationMap::value_type(pParentCluster, clusterAssociation));
            }
        }
    }

    // Loop over daughter clusters and select the nearest parent clusters
    ClusterVector sortedDaughterClusters;
    for (const auto &mapEntry : daughterToParentMatrix)
        sortedDaughterClusters.push_back(mapEntry.first);
    std::sort(sortedDaughterClusters.begin(), sortedDaughterClusters.end(), LArClusterHelper::SortByNHits);

    // Loop over parent clusters and select nearby daughter clusters that are closer than another parent cluster
    for (const Cluster *const pDaughterCluster : sortedDaughterClusters)
    {
        const ClusterAssociationMap &parentToAssociationMap(daughterToParentMatrix.at(pDaughterCluster));

        const Cluster *pParentCluster(nullptr);
        float minDisplacement(std::numeric_limits<float>::max());

        ClusterVector sortedLocalParentClusters;
        for (const auto &mapEntry : parentToAssociationMap)
            sortedLocalParentClusters.push_back(mapEntry.first);
        std::sort(sortedLocalParentClusters.begin(), sortedLocalParentClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pCandidateParentCluster : sortedLocalParentClusters)
        {
            const ClusterAssociation &clusterAssociation(parentToAssociationMap.at(pCandidateParentCluster));

            if (clusterAssociation.GetFigureOfMerit() < minDisplacement)
            {
                minDisplacement = clusterAssociation.GetFigureOfMerit();
                pParentCluster = pCandidateParentCluster;
            }
        }

        if (pParentCluster)
        {
            ClusterList &parentList(clusterMergeMap[pParentCluster]);

            if (parentList.end() == std::find(parentList.begin(), parentList.end(), pDaughterCluster))
                parentList.push_back(pDaughterCluster);

            ClusterList &daughterList(clusterMergeMap[pDaughterCluster]);

            if (daughterList.end() == std::find(daughterList.begin(), daughterList.end(), pParentCluster))
                daughterList.push_back(pParentCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
