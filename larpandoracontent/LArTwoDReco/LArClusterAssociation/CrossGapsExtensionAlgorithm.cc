/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the cross gaps extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CrossGapsExtensionAlgorithm::CrossGapsExtensionAlgorithm() :
    m_minClusterLength(10.f),
    m_maxGapTolerance(2.5f),
    m_maxTransverseDisplacement(2.5f),
    m_maxRelativeAngle(15.f),
    m_minCosRelativeAngle(0.966)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // ATTN - Opt out completely if there is no gap information available
    if (PandoraContentApi::GetGeometry(*this)->GetDetectorGapList().empty())
         return;

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    // Build and loop over list of pointing clusters
    LArPointingClusterList pointingClusterList;
    this->BuildPointingClusterList(clusterVector, pointingClusterList);

    for (LArPointingClusterList::const_iterator iter1 = pointingClusterList.begin(), iterEnd1 = pointingClusterList.end();
        iter1 != iterEnd1; ++iter1)
    {
        const LArPointingCluster &pointingCluster1 = *iter1;
        const Cluster *const pCluster1(pointingCluster1.GetCluster());

        for (LArPointingClusterList::const_iterator iter2 = iter1, iterEnd2 = pointingClusterList.end(); iter2 != iterEnd2; ++iter2)
        {
            const LArPointingCluster &pointingCluster2 = *iter2;
            const Cluster *const pCluster2(pointingCluster2.GetCluster());

            if (pCluster1 == pCluster2)
                continue;

            // Get hit types and check they're the same
            const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

            if (hitType1 != hitType2)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            // Get closest pair of vertices from pointing clusters
            LArPointingCluster::Vertex closestVertex1, closestVertex2;

            try
            {
                LArPointingClusterHelper::GetClosestVertices(pointingCluster1, pointingCluster2, closestVertex1, closestVertex2);
            }
            catch (StatusCodeException &)
            {
                continue;
            }

            // Check that these vertices have associated proximity and pointing information
            if (!this->IsAssociated(closestVertex1, closestVertex2))
                continue;

            // Check that these vertices lie across a registered gap
            if (!(LArGeometryHelper::IsInGap(this->GetPandora(), closestVertex1.GetPosition(), closestVertex2.GetPosition(),
                hitType1, m_maxGapTolerance)))
                continue;

            // Calculate figure of merit for this pair of clusters and require a positive answer
            const float clusterLength1(LArClusterHelper::GetLength(pCluster1));
            const float clusterLength2(LArClusterHelper::GetLength(pCluster2));
            const float clusterSeparation((closestVertex1.GetPosition() - closestVertex2.GetPosition()).GetMagnitude());

            const float cosTheta(-1.f * closestVertex1.GetDirection().GetDotProduct(closestVertex2.GetDirection()));
            const float logTheta(-1.f * std::log10(1.0 - cosTheta + std::numeric_limits<float>::epsilon()));

            const float figureOfMerit12(logTheta * (2.f * clusterLength2 - clusterSeparation));
            const float figureOfMerit21(logTheta * (2.f * clusterLength1 - clusterSeparation));

            if (figureOfMerit12 + figureOfMerit21 < std::numeric_limits<float>::epsilon())
                continue;

            const ClusterAssociation::VertexType vertexType1(closestVertex1.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
            const ClusterAssociation::VertexType vertexType2(closestVertex2.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);

            (void) clusterAssociationMatrix[pCluster1].insert(ClusterAssociationMap::value_type(pCluster2,
                ClusterAssociation(vertexType1, vertexType2, ClusterAssociation::STRONG, figureOfMerit12)));
            (void) clusterAssociationMatrix[pCluster2].insert(ClusterAssociationMap::value_type(pCluster1,
                ClusterAssociation(vertexType2, vertexType1, ClusterAssociation::STRONG, figureOfMerit21)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsExtensionAlgorithm::BuildPointingClusterList(const ClusterVector &clusterVector, LArPointingClusterList &pointingClusterList) const
{
    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            pointingClusterList.push_back(LArPointingCluster(pCluster));
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsExtensionAlgorithm::IsAssociated(const LArPointingCluster::Vertex &pointingVertex1, const LArPointingCluster::Vertex &pointingVertex2) const
{
    const float maxLongitudinalDisplacement((pointingVertex2.GetPosition() - pointingVertex1.GetPosition()).GetMagnitude());

    const bool isAssociated1(LArPointingClusterHelper::IsEmission(pointingVertex1.GetPosition(), pointingVertex2, -1.f,
        maxLongitudinalDisplacement + 1.f, m_maxTransverseDisplacement, m_maxRelativeAngle));
    const bool isAssociated2(LArPointingClusterHelper::IsEmission(pointingVertex2.GetPosition(), pointingVertex1, -1.f,
        maxLongitudinalDisplacement + 1.f, m_maxTransverseDisplacement, m_maxRelativeAngle));
    const bool isAssociated3(-1.f * pointingVertex1.GetDirection().GetDotProduct(pointingVertex2.GetDirection()) > m_minCosRelativeAngle);

    return (isAssociated1 && isAssociated2 && isAssociated3);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &inputAssociationMatrix, ClusterMergeMap &clusterMergeMap) const
{
    // Decide which associations will become merges
    // To make the merge A <-> B, both A -> B and B -> A must be strong associations
    // with the largest figures of merit of all the A -> X and B -> Y associations

    // First step: remove double-counting from the map of associations
    // i.e. if the map has A <-> B, B <-> C, A <-> C, then remove A <-> C
    ClusterAssociationMatrix clusterAssociationMatrix;

    ClusterVector sortedInputClusters;
    for (const auto &mapEntry : inputAssociationMatrix) sortedInputClusters.push_back(mapEntry.first);
    std::sort(sortedInputClusters.begin(), sortedInputClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster1 : sortedInputClusters)
    {
        const ClusterAssociationMap &associationMap1(inputAssociationMatrix.at(pCluster1));

        for (const Cluster *const pCluster2 : sortedInputClusters)
        {
            if (pCluster1 == pCluster2)
                continue;

            const ClusterAssociationMap &associationMap2(inputAssociationMatrix.at(pCluster2));

            ClusterAssociationMap::const_iterator iter12 = associationMap1.find(pCluster2);
            if (associationMap1.end() == iter12)
                continue;

            ClusterAssociationMap::const_iterator iter21 = associationMap2.find(pCluster1);
            if (associationMap2.end() == iter21)
                continue;

            const ClusterAssociation &association12(iter12->second);
            const ClusterAssociation &association21(iter21->second);

            bool isAssociated(true);

            ClusterVector sortedAssociationClusters;
            for (const auto &mapEntry : associationMap1) sortedAssociationClusters.push_back(mapEntry.first);
            std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

            for (const Cluster *const pCluster3 : sortedAssociationClusters)
            {
                const ClusterAssociation &association13(associationMap1.at(pCluster3));

                ClusterAssociationMap::const_iterator iter23 = associationMap2.find(pCluster3);
                if (associationMap2.end() == iter23)
                    continue;

                const ClusterAssociation &association23(iter23->second);

                if (association12.GetParent() == association13.GetParent() &&
                    association23.GetParent() == association21.GetParent() &&
                    association13.GetDaughter() != association23.GetDaughter())
                {
                    isAssociated = false;
                    break;
                }
            }

            if (isAssociated)
            {
                (void) clusterAssociationMatrix[pCluster1].insert(ClusterAssociationMap::value_type(pCluster2, association12));
                (void) clusterAssociationMatrix[pCluster2].insert(ClusterAssociationMap::value_type(pCluster1, association21));
            }
        }
    }


    // Second step: find the best associations A -> X and B -> Y
    ClusterAssociationMatrix intermediateAssociationMatrix;

    ClusterVector sortedClusters;
    for (const auto &mapEntry : clusterAssociationMatrix) sortedClusters.push_back(mapEntry.first);
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : sortedClusters)
    {
        const ClusterAssociationMap &clusterAssociationMap(clusterAssociationMatrix.at(pParentCluster));

        const Cluster *pBestClusterInner(nullptr);
        ClusterAssociation bestAssociationInner(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        const Cluster *pBestClusterOuter(nullptr);
        ClusterAssociation bestAssociationOuter(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        ClusterVector sortedAssociationClusters;
        for (const auto &mapEntry : clusterAssociationMap) sortedAssociationClusters.push_back(mapEntry.first);
        std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedAssociationClusters)
        {
            const ClusterAssociation &clusterAssociation(clusterAssociationMap.at(pDaughterCluster));

            // Inner associations
            if (clusterAssociation.GetParent() == ClusterAssociation::INNER)
            {
                if (clusterAssociation.GetFigureOfMerit() > bestAssociationInner.GetFigureOfMerit())
                {
                    bestAssociationInner = clusterAssociation;
                    pBestClusterInner = pDaughterCluster;
                }
            }

            // Outer associations
            if (clusterAssociation.GetParent() == ClusterAssociation::OUTER)
            {
                if (clusterAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit())
                {
                    bestAssociationOuter = clusterAssociation;
                    pBestClusterOuter = pDaughterCluster;
                }
            }
        }

        if (pBestClusterInner)
            (void) intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterInner, bestAssociationInner));

        if (pBestClusterOuter)
            (void) intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterOuter, bestAssociationOuter));
    }


    // Third step: make the merge if A -> X and B -> Y is in fact A -> B and B -> A
    ClusterVector intermediateSortedClusters;
    for (const auto &mapEntry : intermediateAssociationMatrix) intermediateSortedClusters.push_back(mapEntry.first);
    std::sort(intermediateSortedClusters.begin(), intermediateSortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : intermediateSortedClusters)
    {
        const ClusterAssociationMap &parentAssociationMap(intermediateAssociationMatrix.at(pParentCluster));

        ClusterVector sortedAssociationClusters;
        for (const auto &mapEntry : parentAssociationMap) sortedAssociationClusters.push_back(mapEntry.first);
        std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedAssociationClusters)
        {
            const ClusterAssociation &parentToDaughterAssociation(parentAssociationMap.at(pDaughterCluster));

            ClusterAssociationMatrix::const_iterator iter5 = intermediateAssociationMatrix.find(pDaughterCluster);

            if (intermediateAssociationMatrix.end() == iter5)
                continue;

            const ClusterAssociationMap &daughterAssociationMap(iter5->second);

            ClusterAssociationMap::const_iterator iter6 = daughterAssociationMap.find(pParentCluster);

            if (daughterAssociationMap.end() == iter6)
                continue;

            const ClusterAssociation &daughterToParentAssociation(iter6->second);

            if (parentToDaughterAssociation.GetParent() == daughterToParentAssociation.GetDaughter() &&
                parentToDaughterAssociation.GetDaughter() == daughterToParentAssociation.GetParent())
            {
                ClusterList &parentList(clusterMergeMap[pParentCluster]);

                if (parentList.end() == std::find(parentList.begin(), parentList.end(), pDaughterCluster))
                    parentList.push_back(pDaughterCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CrossGapsExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGapTolerance", m_maxGapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAngularAllowance", m_maxRelativeAngle));
    m_minCosRelativeAngle = std::cos(m_maxRelativeAngle * M_PI / 180.0);

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
