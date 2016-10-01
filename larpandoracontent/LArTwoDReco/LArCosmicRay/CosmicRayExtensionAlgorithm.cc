/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic-ray extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CosmicRayExtensionAlgorithm::CosmicRayExtensionAlgorithm() :
    m_minClusterLength(3.f),
    m_minSeedClusterLength(6.f),
    m_maxLongitudinalDisplacement(30.f),
    m_maxTransverseDisplacement(2.f),
    m_minCosRelativeAngle(0.966f),
    m_maxAverageRms(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    // Convert each input cluster into a pointing cluster
    LArPointingClusterList pointingClusterList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        try
        {
            pointingClusterList.push_back(LArPointingCluster(*iter));
        }
        catch (StatusCodeException &)
        {
        }
    }

    // Form associations between pairs of pointing clusters
    for (LArPointingClusterList::const_iterator iterI = pointingClusterList.begin(), iterEndI = pointingClusterList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster &clusterI = *iterI;

        for (LArPointingClusterList::const_iterator iterJ = iterI, iterEndJ = pointingClusterList.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster &clusterJ = *iterJ;

            if (clusterI.GetCluster() == clusterJ.GetCluster())
            continue;

            this->FillClusterAssociationMatrix(clusterI, clusterJ, clusterAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayExtensionAlgorithm::FillClusterAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ,
    ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    const Cluster *const pClusterI(clusterI.GetCluster());
    const Cluster *const pClusterJ(clusterJ.GetCluster());

    if (pClusterI == pClusterJ)
        return;

    // Get closest pair of vertices
    LArPointingCluster::Vertex clusterVertexI, clusterVertexJ;

    try
    {
        LArPointingClusterHelper::GetClosestVertices(clusterI, clusterJ, clusterVertexI, clusterVertexJ);
    }
    catch (StatusCodeException &)
    {
        return;
    }

    // (Just in case...)
    if (!(clusterVertexI.IsInitialized() && clusterVertexJ.IsInitialized()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const CartesianVector vertexI(clusterVertexI.GetPosition());
    const CartesianVector vertexJ(clusterVertexJ.GetPosition());

    const CartesianVector endI(clusterVertexI.IsInnerVertex() ? clusterI.GetOuterVertex().GetPosition() : clusterI.GetInnerVertex().GetPosition());
    const CartesianVector endJ(clusterVertexJ.IsInnerVertex() ? clusterJ.GetOuterVertex().GetPosition() : clusterJ.GetInnerVertex().GetPosition());

    // Requirements on length
    const float lengthSquaredI(LArPointingClusterHelper::GetLengthSquared(clusterI));
    const float lengthSquaredJ(LArPointingClusterHelper::GetLengthSquared(clusterJ));
    
    if (std::max(lengthSquaredI, lengthSquaredJ) < m_minSeedClusterLength * m_minSeedClusterLength)
        return;

    // Requirements on proximity
    const float distanceSquaredIJ((vertexI - vertexJ).GetMagnitudeSquared());

    if (distanceSquaredIJ > m_maxLongitudinalDisplacement * m_maxLongitudinalDisplacement)
        return;

    // Requirements on pointing information
    const CartesianVector directionI((endI - vertexI).GetUnitVector());
    const CartesianVector directionJ((endJ - vertexJ).GetUnitVector());

    const float cosTheta(-directionI.GetDotProduct(directionJ));
    const float cosThetaCut(-1.f + 2.f * m_minCosRelativeAngle);

    if (cosTheta < cosThetaCut)
        return;

    // Requirements on overlap between clusters
    const CartesianVector directionIJ((endJ - endI).GetUnitVector());
    const CartesianVector directionJI((endI - endJ).GetUnitVector());

    const float overlapL(directionIJ.GetDotProduct(vertexJ - vertexI));
    const float overlapT(directionIJ.GetCrossProduct(vertexJ - vertexI).GetMagnitude());

    if (overlapL < -1.f || overlapL * overlapL > 2.f * std::min(lengthSquaredI, lengthSquaredJ) ||
        overlapT > m_maxTransverseDisplacement + std::fabs(overlapL) * std::tan(5.f * M_PI / 180.f))
        return;

    // Calculate RMS deviations on composite cluster
    const float rms1(this->CalculateRms(pClusterI, endI, directionIJ));
    const float rms2(this->CalculateRms(pClusterJ, endJ, directionJI));

    const float rms(0.5f * (rms1 + rms2));
    const float rmsCut(2.f * m_maxAverageRms * (cosTheta - cosThetaCut) / (1.0 - cosThetaCut));

    if (rms > rmsCut)
        return;

    // Record the association
    const ClusterAssociation::VertexType vertexTypeI(clusterVertexI.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
    const ClusterAssociation::VertexType vertexTypeJ(clusterVertexJ.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
    const ClusterAssociation::AssociationType associationType(ClusterAssociation::STRONG);

    (void) clusterAssociationMatrix[pClusterI].insert(ClusterAssociationMap::value_type(pClusterJ,
        ClusterAssociation(vertexTypeI, vertexTypeJ, associationType, lengthSquaredJ)));

    (void) clusterAssociationMatrix[pClusterJ].insert(ClusterAssociationMap::value_type(pClusterI,
        ClusterAssociation(vertexTypeJ, vertexTypeI, associationType, lengthSquaredI)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &inputAssociationMatrix, ClusterMergeMap &clusterMergeMap) const
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

float CosmicRayExtensionAlgorithm::CalculateRms(const Cluster *const pCluster, const CartesianVector &position, const CartesianVector &direction) const
{
    float totalChi2(0.f);
    float totalHits(0.f);

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);

    for (CaloHitList::const_iterator iter = caloHitList.begin(), iterEnd = caloHitList.end(); iter != iterEnd; ++iter)
    {
        const CaloHit *const pCaloHit = *iter;

        const CartesianVector hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector predictedPosition(position + direction * direction.GetDotProduct(hitPosition - position));

        totalChi2 += (predictedPosition - hitPosition).GetMagnitudeSquared();
        totalHits += 1.f;
    }

    if (totalHits > 0.f)
        return std::sqrt(totalChi2/totalHits);

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeedClusterLength", m_minSeedClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAverageRms", m_maxAverageRms));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
