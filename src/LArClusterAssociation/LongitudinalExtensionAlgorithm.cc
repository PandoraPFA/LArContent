/**
 *  @file   LArContent/src/LArClusterAssociation/LongitudinalExtensionAlgorithm.cc
 * 
 *  @brief  Implementation of the longitudinal extension algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/LongitudinalExtensionAlgorithm.h"

#include "LArHelpers/LArVertexHelper.h"
#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

void LongitudinalExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < 15)
            continue;

        const CartesianVector innerVertex(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
        const CartesianVector outerVertex(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));

        if ((outerVertex - innerVertex).GetMagnitudeSquared() < 25.f)
            continue;

        if (LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75f)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    LArPointingClusterList pointingClusterList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        pointingClusterList.push_back(LArPointingCluster(*iter));
    }

    for (LArPointingClusterList::const_iterator iterI = pointingClusterList.begin(), iterEndI = pointingClusterList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster &clusterI = *iterI;

        for (LArPointingClusterList::const_iterator iterJ = iterI, iterEndJ = pointingClusterList.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster &clusterJ = *iterJ;

            if (clusterI.GetCluster() == clusterJ.GetCluster())
                continue;

            this->FillAssociationMatrix(clusterI, clusterJ, clusterAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ,
    ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    if (clusterI.GetCluster() == clusterJ.GetCluster())
        return;

    // Check association of closest two end points
    const LArPointingCluster::Vertex &innerVertexI(clusterI.GetInnerVertex());
    const LArPointingCluster::Vertex &outerVertexI(clusterI.GetOuterVertex());
    const LArPointingCluster::Vertex &innerVertexJ(clusterJ.GetInnerVertex());
    const LArPointingCluster::Vertex &outerVertexJ(clusterJ.GetOuterVertex());

    const float distSquared_innerI_to_innerJ((innerVertexI.GetPosition()-innerVertexJ.GetPosition()).GetMagnitudeSquared());
    const float distSquared_innerI_to_outerJ((innerVertexI.GetPosition()-outerVertexJ.GetPosition()).GetMagnitudeSquared());
    const float distSquared_outerI_to_innerJ((outerVertexI.GetPosition()-innerVertexJ.GetPosition()).GetMagnitudeSquared());
    const float distSquared_outerI_to_outerJ((outerVertexI.GetPosition()-outerVertexJ.GetPosition()).GetMagnitudeSquared());

    // (a) Check association of inner I <-> inner J
    if ((distSquared_innerI_to_innerJ < std::min(distSquared_outerI_to_outerJ, std::min(distSquared_innerI_to_outerJ, distSquared_outerI_to_innerJ))) &&
        (distSquared_outerI_to_outerJ > std::max(distSquared_innerI_to_innerJ, std::max(distSquared_innerI_to_outerJ, distSquared_outerI_to_innerJ))) )
    {
        return this->FillAssociationMatrix(innerVertexI, innerVertexJ, clusterAssociationMatrix);
    }

    // (b) Check association of inner I <-> outer J
    if ((distSquared_innerI_to_outerJ < std::min(distSquared_outerI_to_innerJ, std::min(distSquared_outerI_to_outerJ, distSquared_innerI_to_innerJ))) &&
        (distSquared_outerI_to_innerJ > std::max(distSquared_innerI_to_outerJ, std::max(distSquared_outerI_to_outerJ, distSquared_innerI_to_innerJ))) )
    {
        return this->FillAssociationMatrix(innerVertexI, outerVertexJ, clusterAssociationMatrix);
    }

    // (c) Check association of outer I <-> inner J
    if ((distSquared_outerI_to_innerJ < std::min(distSquared_innerI_to_outerJ, std::min(distSquared_innerI_to_innerJ, distSquared_outerI_to_outerJ))) &&
        (distSquared_innerI_to_outerJ > std::max(distSquared_outerI_to_innerJ, std::max(distSquared_innerI_to_innerJ, distSquared_outerI_to_outerJ))) )
    {
        return this->FillAssociationMatrix(outerVertexI, innerVertexJ, clusterAssociationMatrix);
    }

    // (d) Check association of outer I <-> outer J
    if ((distSquared_outerI_to_outerJ < std::min(distSquared_innerI_to_innerJ, std::min(distSquared_innerI_to_outerJ, distSquared_outerI_to_innerJ))) &&
        (distSquared_innerI_to_innerJ > std::max(distSquared_outerI_to_outerJ, std::max(distSquared_innerI_to_outerJ, distSquared_outerI_to_innerJ))) )
    {
        return this->FillAssociationMatrix(outerVertexI, outerVertexJ, clusterAssociationMatrix);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const LArPointingCluster::Vertex &clusterVertexI, const LArPointingCluster::Vertex &clusterVertexJ,
    ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    const Cluster *const pClusterI(clusterVertexI.GetCluster());
    const Cluster *const pClusterJ(clusterVertexJ.GetCluster());

    if (pClusterI == pClusterJ)
        return;

    if ((clusterVertexI.GetRms() > 0.5f) || (clusterVertexJ.GetRms() > 0.5f))
        return;

    // Check that new layer occupancy would be reasonable
    if (LArClusterHelper::GetLayerOccupancy(pClusterI, pClusterJ) < 0.75)
        return;

    // Association types
    ClusterAssociation::AssociationType associationType(ClusterAssociation::NOTHING);
    ClusterAssociation::StrengthType strengthType(ClusterAssociation::UNASSOCIATED);

    // Requirements on cluster end points
    const CartesianVector &vertexPositionI(clusterVertexI.GetPosition());
    const CartesianVector &vertexPositionJ(clusterVertexJ.GetPosition());
    const CartesianVector &vertexDirectionI(clusterVertexI.GetDirection());
    const CartesianVector &vertexDirectionJ(clusterVertexJ.GetDirection());

    const float distanceSquared((vertexPositionI - vertexPositionJ).GetMagnitudeSquared());

    if (distanceSquared < 2.f * m_spatialResolution * m_spatialResolution)
    {
        associationType = ClusterAssociation::NODE;
        strengthType = ClusterAssociation::WEAK;

        if (distanceSquared < m_spatialResolution * m_spatialResolution)
        {
            const float cosTheta(-vertexDirectionI.GetDotProduct(vertexDirectionJ));

            if (cosTheta > 0.906f)
            {
                strengthType = ClusterAssociation::STANDARD;

                if (cosTheta > 0.966f)
                    strengthType = ClusterAssociation::STRONG;
            }
        }
    }

    // Calculate lengths of clusters
    const CartesianVector innerCentroidI(pClusterI->GetCentroid(pClusterI->GetInnerPseudoLayer()));
    const CartesianVector outerCentroidI(pClusterI->GetCentroid(pClusterI->GetOuterPseudoLayer()));
    const CartesianVector innerCentroidJ(pClusterJ->GetCentroid(pClusterJ->GetInnerPseudoLayer()));
    const CartesianVector outerCentroidJ(pClusterJ->GetCentroid(pClusterJ->GetOuterPseudoLayer()));

    const float clusterOccupancyI(LArClusterHelper::GetLayerOccupancy(pClusterI));
    const float clusterOccupancyJ(LArClusterHelper::GetLayerOccupancy(pClusterJ));
    const float clusterLengthI(clusterOccupancyI * (outerCentroidI - innerCentroidI).GetMagnitude());
    const float clusterLengthJ(clusterOccupancyJ * (outerCentroidJ - innerCentroidJ).GetMagnitude());

    if ((clusterLengthI > 10.f) && (clusterLengthJ > 10.f))
    {
        // Calculate impact parameters
        if (strengthType < ClusterAssociation::STRONG)
        {
            const float cosTheta(-vertexDirectionI.GetDotProduct(vertexDirectionJ));
            const float cosThetaI((vertexPositionI - vertexPositionJ).GetUnitVector().GetDotProduct(vertexDirectionI));
            const float cosThetaJ((vertexPositionJ - vertexPositionI).GetUnitVector().GetDotProduct(vertexDirectionJ));

            float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
            LArVertexHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, vertexPositionJ, rL1, rT1);
            LArVertexHelper::GetImpactParameters(vertexPositionJ, vertexDirectionJ, vertexPositionI, rL2, rT2);

            if ((cosTheta > 0.985f) && (std::fabs(cosThetaI) > 0.25f) && (std::fabs(cosThetaJ) > 0.25f) && (std::fabs(rL1) < 10.f) &&
                (std::fabs(rL2) < 10.f) && (rT1 < 2.f * m_spatialResolution) && (rT2 < 2.f * m_spatialResolution))
            {
                associationType = ClusterAssociation::EMISSION;
                strengthType = ClusterAssociation::STRONG;
            }
        }

        // Remove 10 layers and re-calculate impact parameters
        if (strengthType < ClusterAssociation::STRONG)
        {
            const LArPointingCluster::Vertex truncatedVertexI(clusterVertexI.GetCluster(), clusterVertexI.IsInnerVertex(), 10);
            const float cosTheta(-truncatedVertexI.GetDirection().GetDotProduct(vertexDirectionJ));

            float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
            LArVertexHelper::GetImpactParameters(truncatedVertexI.GetPosition(), truncatedVertexI.GetDirection(), vertexPositionJ, rL1, rT1);
            LArVertexHelper::GetImpactParameters(vertexPositionJ, vertexDirectionJ, truncatedVertexI.GetPosition(), rL2, rT2);

            if ((cosTheta > 0.985f) && (truncatedVertexI.GetRms() < 0.5f) && (std::fabs(rL1) < 10.f) && (std::fabs(rL2) < 10.f) &&
                (rT1 < m_spatialResolution) && (rT2 < m_spatialResolution))
            {
                associationType = ClusterAssociation::EMISSION;
                strengthType = ClusterAssociation::STRONG;
            }
        }

        if (strengthType < ClusterAssociation::STRONG)
        {
            const LArPointingCluster::Vertex truncatedVertexJ(clusterVertexJ.GetCluster(), clusterVertexJ.IsInnerVertex(), 10);
            const float cosTheta(-truncatedVertexJ.GetDirection().GetDotProduct(vertexDirectionI));

            float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
            LArVertexHelper::GetImpactParameters(truncatedVertexJ.GetPosition(), truncatedVertexJ.GetDirection(), vertexPositionI, rL1, rT1);
            LArVertexHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, truncatedVertexJ.GetPosition(), rL2, rT2);

            if ((cosTheta > 0.985f) && (truncatedVertexJ.GetRms() < 0.5f) && (std::fabs(rL1) < 10.f) && (std::fabs(rL2) < 10.f) &&
                (rT1 < m_spatialResolution) && (rT2 < m_spatialResolution))
            {
                associationType = ClusterAssociation::EMISSION;
                strengthType = ClusterAssociation::STRONG;
            }
        }
    }

    if (strengthType > ClusterAssociation::UNASSOCIATED)
    {
        const ClusterAssociation::VertexType vertexTypeI(clusterVertexI.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
        const ClusterAssociation::VertexType vertexTypeJ(clusterVertexJ.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
        (void) clusterAssociationMatrix[pClusterI].insert(ClusterAssociationMap::value_type(pClusterJ, ClusterAssociation(vertexTypeI, vertexTypeJ, associationType, strengthType, clusterLengthJ)));
        (void) clusterAssociationMatrix[pClusterJ].insert(ClusterAssociationMap::value_type(pClusterI, ClusterAssociation(vertexTypeJ, vertexTypeI, associationType, strengthType, clusterLengthI)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalExtensionAlgorithm::AreClustersAssociated(Cluster *pCluster1, Cluster *pCluster2, const ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    // Search for the pCluster1 <-> pCluster2 association
    ClusterAssociationMatrix::const_iterator iter1A = clusterAssociationMatrix.find(pCluster1);

    if (clusterAssociationMatrix.end() == iter1A)
        return false;

    const ClusterAssociationMap &clusterAssociationMap1(iter1A->second);

    ClusterAssociationMap::const_iterator iter1B = clusterAssociationMap1.find(pCluster2);

    if (clusterAssociationMap1.end() == iter1B)
        return false;

    const ClusterAssociation &clusterAssociationCandidate1(iter1B->second);

    if (clusterAssociationCandidate1.GetStrength() < ClusterAssociation::STANDARD)
        return false;

    // Bail out if there is a stronger association, or if this is a prong
    for (ClusterAssociationMap::const_iterator iter1C = clusterAssociationMap1.begin(), iter1CEnd = clusterAssociationMap1.end(); iter1C != iter1CEnd; ++iter1C)
    {
        const Cluster *pClusterAlternative1(iter1C->first);

        if (pClusterAlternative1 == pCluster2)
            continue;

        const ClusterAssociation &clusterAssociationAlternative1(iter1C->second);

        if (clusterAssociationCandidate1.GetParent() == clusterAssociationAlternative1.GetParent())
        {
            if (clusterAssociationAlternative1.GetAssociation() == ClusterAssociation::NODE)
                return false;

            if ((clusterAssociationAlternative1.GetStrength() > clusterAssociationCandidate1.GetStrength()) ||
                ((clusterAssociationAlternative1.GetStrength() == clusterAssociationCandidate1.GetStrength()) &&
                (clusterAssociationAlternative1.GetFigureOfMerit() > clusterAssociationCandidate1.GetFigureOfMerit())) )
            {
                return false;
            }
        }
    }

    // Search for the pCluster2 <-> pCluster1 association
    ClusterAssociationMatrix::const_iterator iter2A = clusterAssociationMatrix.find(pCluster2);

    if (clusterAssociationMatrix.end() == iter2A)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ClusterAssociationMap &clusterAssociationMap2(iter2A->second);

    ClusterAssociationMap::const_iterator iter2B = clusterAssociationMap2.find(pCluster1);

    if (clusterAssociationMap2.end() == iter2B)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ClusterAssociation &clusterAssociationCandidate2(iter2B->second);

    if (clusterAssociationCandidate2.GetStrength() < ClusterAssociation::STANDARD)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    // Bail out if there is a stronger association, or if this is a prong
    for (ClusterAssociationMap::const_iterator iter2C = clusterAssociationMap2.begin(), iter2CEnd = clusterAssociationMap2.end(); iter2C != iter2CEnd; ++iter2C)
    {
        const Cluster *pClusterAlternative2(iter2C->first);

        if (pClusterAlternative2 == pCluster1)
            continue;

        const ClusterAssociation &clusterAssociationAlternative2(iter2C->second);

        if (clusterAssociationCandidate2.GetParent() == clusterAssociationAlternative2.GetParent())
        {
            if (clusterAssociationAlternative2.GetAssociation() == ClusterAssociation::NODE)
                return false;

            if ((clusterAssociationAlternative2.GetStrength() > clusterAssociationCandidate2.GetStrength()) ||
                ((clusterAssociationAlternative2.GetStrength() == clusterAssociationCandidate2.GetStrength()) &&
                (clusterAssociationAlternative2.GetFigureOfMerit() > clusterAssociationCandidate2.GetFigureOfMerit())) )
            {
                return false;
            }
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_spatialResolution = 1.16f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SpatialResolution", m_spatialResolution));

    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
