/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the cross gaps extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CrossGapsExtensionAlgorithm::CrossGapsExtensionAlgorithm() :
    m_minClusterLength(5.f), m_minGapFraction(0.5f), m_maxGapTolerance(2.f), m_maxTransverseDisplacement(2.5f), m_maxRelativeAngle(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // ATTN May want to opt-out completely if no gap information available
    // if (PandoraContentApi::GetGeometry(*this)->GetDetectorGapList().empty())
    //     return;

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
    // Build lists of pointing clusters in proximity to gaps
    LArPointingClusterList innerPointingClusterList, outerPointingClusterList;
    this->BuildPointingClusterList(clusterVector, innerPointingClusterList, outerPointingClusterList);

    // Form associations between pairs of pointing clusters
    for (const LArPointingCluster &pointingClusterInner : innerPointingClusterList)
    {
        const LArPointingCluster::Vertex &pointingVertexInner(pointingClusterInner.GetInnerVertex());
        const float zInner(pointingVertexInner.GetPosition().GetZ());

        for (const LArPointingCluster &pointingClusterOuter : outerPointingClusterList)
        {
            const LArPointingCluster::Vertex &pointingVertexOuter(pointingClusterOuter.GetOuterVertex());
            const float zOuter(pointingVertexOuter.GetPosition().GetZ());

            if (!this->IsAcrossGap(zOuter, zInner, LArClusterHelper::GetClusterHitType(pointingClusterInner.GetCluster())))
                continue;

            if (!this->IsAssociated(pointingVertexInner, pointingVertexOuter))
                continue;

            const Cluster *const pClusterInner(pointingClusterInner.GetCluster());
            const Cluster *const pClusterOuter(pointingClusterOuter.GetCluster());

            const float lengthSquaredInner(LArClusterHelper::GetLengthSquared(pClusterInner));
            const float lengthSquaredOuter(LArClusterHelper::GetLengthSquared(pClusterOuter));

            (void)clusterAssociationMatrix[pClusterInner].insert(ClusterAssociationMap::value_type(pClusterOuter,
                ClusterAssociation(ClusterAssociation::INNER, ClusterAssociation::OUTER, ClusterAssociation::STRONG, lengthSquaredOuter)));
            (void)clusterAssociationMatrix[pClusterOuter].insert(ClusterAssociationMap::value_type(pClusterInner,
                ClusterAssociation(ClusterAssociation::OUTER, ClusterAssociation::INNER, ClusterAssociation::STRONG, lengthSquaredInner)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsExtensionAlgorithm::BuildPointingClusterList(const ClusterVector &clusterVector,
    LArPointingClusterList &innerPointingClusterList, LArPointingClusterList &outerPointingClusterList) const
{
    // Convert each input cluster into a pointing cluster
    LArPointingClusterList pointingClusterList;

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

    // Identify clusters adjacent to detector gaps
    this->BuildPointingClusterList(true, pointingClusterList, innerPointingClusterList);
    this->BuildPointingClusterList(false, pointingClusterList, outerPointingClusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsExtensionAlgorithm::BuildPointingClusterList(
    const bool useInner, const LArPointingClusterList &inputPointingClusterList, LArPointingClusterList &outputPointingClusterList) const
{
    for (const LArPointingCluster &pointingCluster : inputPointingClusterList)
    {
        const LArPointingCluster::Vertex &pointingVertex(useInner ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());
        const HitType hitType(LArClusterHelper::GetClusterHitType(pointingCluster.GetCluster()));

        if (LArGeometryHelper::IsInGap(this->GetPandora(), pointingVertex.GetPosition(), hitType, m_maxGapTolerance))
            outputPointingClusterList.push_back(pointingCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsExtensionAlgorithm::IsAssociated(const LArPointingCluster::Vertex &pointingVertex1, const LArPointingCluster::Vertex &pointingVertex2) const
{
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), LArClusterHelper::GetClusterHitType(pointingVertex1.GetCluster()))};
    const float maxTransverseDisplacementAdjusted{ratio * m_maxTransverseDisplacement};
    const float maxLongitudinalDisplacement((pointingVertex2.GetPosition() - pointingVertex1.GetPosition()).GetMagnitude());

    const bool isAssociated1(LArPointingClusterHelper::IsEmission(pointingVertex1.GetPosition(), pointingVertex2, -1.f,
        maxLongitudinalDisplacement + 1.f, maxTransverseDisplacementAdjusted, m_maxRelativeAngle));
    const bool isAssociated2(LArPointingClusterHelper::IsEmission(pointingVertex2.GetPosition(), pointingVertex1, -1.f,
        maxLongitudinalDisplacement + 1.f, maxTransverseDisplacementAdjusted, m_maxRelativeAngle));

    return (isAssociated1 && isAssociated2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsExtensionAlgorithm::IsAcrossGap(const float minZ, const float maxZ, const HitType hitType) const
{
    if (maxZ - minZ < std::numeric_limits<float>::epsilon())
        return false;

    const float gapDeltaZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, hitType));

    if (gapDeltaZ / (maxZ - minZ) < m_minGapFraction)
        return false;

    return true;
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
    for (const auto &mapEntry : inputAssociationMatrix)
        sortedInputClusters.push_back(mapEntry.first);
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
            for (const auto &mapEntry : associationMap1)
                sortedAssociationClusters.push_back(mapEntry.first);
            std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

            for (const Cluster *const pCluster3 : sortedAssociationClusters)
            {
                const ClusterAssociation &association13(associationMap1.at(pCluster3));

                ClusterAssociationMap::const_iterator iter23 = associationMap2.find(pCluster3);
                if (associationMap2.end() == iter23)
                    continue;

                const ClusterAssociation &association23(iter23->second);

                if (association12.GetParent() == association13.GetParent() && association23.GetParent() == association21.GetParent() &&
                    association13.GetDaughter() != association23.GetDaughter())
                {
                    isAssociated = false;
                    break;
                }
            }

            if (isAssociated)
            {
                (void)clusterAssociationMatrix[pCluster1].insert(ClusterAssociationMap::value_type(pCluster2, association12));
                (void)clusterAssociationMatrix[pCluster2].insert(ClusterAssociationMap::value_type(pCluster1, association21));
            }
        }
    }

    // Second step: find the best associations A -> X and B -> Y
    ClusterAssociationMatrix intermediateAssociationMatrix;

    ClusterVector sortedClusters;
    for (const auto &mapEntry : clusterAssociationMatrix)
        sortedClusters.push_back(mapEntry.first);
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : sortedClusters)
    {
        const ClusterAssociationMap &clusterAssociationMap(clusterAssociationMatrix.at(pParentCluster));

        const Cluster *pBestClusterInner(nullptr);
        ClusterAssociation bestAssociationInner(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        const Cluster *pBestClusterOuter(nullptr);
        ClusterAssociation bestAssociationOuter(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        ClusterVector sortedAssociationClusters;
        for (const auto &mapEntry : clusterAssociationMap)
            sortedAssociationClusters.push_back(mapEntry.first);
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
            (void)intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterInner, bestAssociationInner));

        if (pBestClusterOuter)
            (void)intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterOuter, bestAssociationOuter));
    }

    // Third step: make the merge if A -> X and B -> Y is in fact A -> B and B -> A
    ClusterVector intermediateSortedClusters;
    for (const auto &mapEntry : intermediateAssociationMatrix)
        intermediateSortedClusters.push_back(mapEntry.first);
    std::sort(intermediateSortedClusters.begin(), intermediateSortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : intermediateSortedClusters)
    {
        const ClusterAssociationMap &parentAssociationMap(intermediateAssociationMatrix.at(pParentCluster));

        ClusterVector sortedAssociationClusters;
        for (const auto &mapEntry : parentAssociationMap)
            sortedAssociationClusters.push_back(mapEntry.first);
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
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinGapFraction", m_minGapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapTolerance", m_maxGapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    float maxCosRelativeAngle(std::cos(m_maxRelativeAngle * M_PI / 180.0));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCosRelativeAngle", maxCosRelativeAngle));
    m_maxRelativeAngle = (180.0 / M_PI) * std::acos(maxCosRelativeAngle);

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
