/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic-ray extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void CosmicRayExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }
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
    // First step: remove double-counting from the map of associations
    // i.e. if the map has A <-> B, B <-> C, A <-> C, then remove A <-> C
    ClusterAssociationMatrix clusterAssociationMatrix;

    for (ClusterAssociationMatrix::const_iterator iter1 = inputAssociationMatrix.begin(), iterEnd1 = inputAssociationMatrix.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster* pCluster1(iter1->first);
        const ClusterAssociationMap &associationMap1(iter1->second);

        for (ClusterAssociationMatrix::const_iterator iter2 = iter1, iterEnd2 = inputAssociationMatrix.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pCluster2(iter2->first);
            const ClusterAssociationMap &associationMap2(iter2->second);

            if (pCluster1 == pCluster2)
            continue;

            ClusterAssociationMap::const_iterator iter12 = associationMap1.find(pCluster2);
            if (associationMap1.end() == iter12)
            continue;

            ClusterAssociationMap::const_iterator iter21 = associationMap2.find(pCluster1);
            if (associationMap2.end() == iter21)
            continue;

            const ClusterAssociation &association12(iter12->second);
            const ClusterAssociation &association21(iter21->second);

            bool isAssociated(true);

            for (ClusterAssociationMap::const_iterator iter13 = associationMap1.begin(), iterEnd13 = associationMap1.end(); iter13 != iterEnd13; ++iter13)
            {
                const Cluster* pCluster3(iter13->first);

                ClusterAssociationMap::const_iterator iter23 = associationMap2.find(pCluster3);
                if (associationMap2.end() == iter23)
                    continue;

                const ClusterAssociation &association13(iter13->second);
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

    for (ClusterAssociationMatrix::const_iterator iter1 = clusterAssociationMatrix.begin(), iterEnd1 = clusterAssociationMatrix.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster* pParentCluster(iter1->first);
        const ClusterAssociationMap &clusterAssociationMap(iter1->second);

        Cluster* pBestClusterInner = NULL;
        ClusterAssociation bestAssociationInner(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        Cluster* pBestClusterOuter = NULL;
        ClusterAssociation bestAssociationOuter(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        for (ClusterAssociationMap::const_iterator iter2 = clusterAssociationMap.begin(), iterEnd2 = clusterAssociationMap.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pDaughterCluster(iter2->first);
            const ClusterAssociation &clusterAssociation(iter2->second);

            // Inner associations
            if (clusterAssociation.GetParent() == ClusterAssociation::INNER)
            {
                if (clusterAssociation.GetFigureOfMerit() > bestAssociationInner.GetFigureOfMerit())
                {
                    bestAssociationInner = clusterAssociation;
                    pBestClusterInner = (Cluster*)pDaughterCluster;
                }
            }

            // Outer associations
            if (clusterAssociation.GetParent() == ClusterAssociation::OUTER)
            {
                if (clusterAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit())
                {
                    bestAssociationOuter = clusterAssociation;
                    pBestClusterOuter = (Cluster*)pDaughterCluster;
                }
            }
        }

        if (pBestClusterInner)
            (void) intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterInner, bestAssociationInner));

        if (pBestClusterOuter)
            (void) intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterOuter, bestAssociationOuter));
    }


    // Third step: make the merge if A -> B and B -> A are mutually the best associations
    for (ClusterAssociationMatrix::const_iterator iter3 = intermediateAssociationMatrix.begin(), iterEnd3 = intermediateAssociationMatrix.end(); iter3 != iterEnd3; ++iter3)
    {
        const Cluster *pParentCluster(iter3->first);
        const ClusterAssociationMap &parentAssociationMap(iter3->second);

        for (ClusterAssociationMap::const_iterator iter4 = parentAssociationMap.begin(), iterEnd4 = parentAssociationMap.end(); iter4 != iterEnd4; ++iter4)
        {
            const Cluster *pDaughterCluster(iter4->first);
            const ClusterAssociation &parentToDaughterAssociation(iter4->second);

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
                clusterMergeMap[pParentCluster].insert((Cluster*)pDaughterCluster);

// ---- BEGIN DISPLAY ----
// ClusterList tempList1, tempList2;
// tempList1.insert((Cluster*)pParentCluster);
// tempList2.insert((Cluster*)pDaughterCluster);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "ParentCluster", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "DaughterCluster", BLUE);
// PandoraMonitoringApi::ViewEvent();
// ---- END DISPLAY ----
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
        CaloHit *pCaloHit = *iter;

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
    m_minClusterLength = 3.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    m_minSeedClusterLength = 6.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeedClusterLength", m_minSeedClusterLength));

    m_maxLongitudinalDisplacement = 30.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    m_maxTransverseDisplacement = 2.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_minCosRelativeAngle = 0.966f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    m_maxAverageRms = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAverageRms", m_maxAverageRms));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
