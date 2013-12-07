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

        if (LArClusterHelper::GetLengthSquared(pCluster) < 25.f)
            continue;

        if (LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75f)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    LongitudinalAssociationMatrix longitudinalAssociationMatrix;
    this->FillAssociationMatrix(clusterVector, longitudinalAssociationMatrix);

    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster *pClusterI = *iterI;

        for (ClusterVector::const_iterator iterJ = clusterVector.begin(), iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster *pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
                continue;

            clusterAssociationMatrix.SetOverlapResult(pClusterI,pClusterJ,NULL,
            this->AreClustersAssociated(pClusterI,pClusterJ,longitudinalAssociationMatrix));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const ClusterVector &clusterVector, LongitudinalAssociationMatrix &longitudinalAssociationMatrix) const
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

            this->FillAssociationMatrix(clusterI, clusterJ, longitudinalAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ,
    LongitudinalAssociationMatrix &longitudinalAssociationMatrix) const
{
    const bool useInnerI(true), useInnerJ(true);
    const bool useOuterI(false), useOuterJ(false);

    this->FillAssociationMatrix(clusterI, clusterJ, useInnerI, useInnerJ, longitudinalAssociationMatrix);
    this->FillAssociationMatrix(clusterI, clusterJ, useInnerI, useOuterJ, longitudinalAssociationMatrix);   
    this->FillAssociationMatrix(clusterI, clusterJ, useOuterI, useInnerJ, longitudinalAssociationMatrix);
    this->FillAssociationMatrix(clusterI, clusterJ, useOuterI, useOuterJ, longitudinalAssociationMatrix);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ, 
    const bool useInnerI, const bool useInnerJ, LongitudinalAssociationMatrix &longitudinalAssociationMatrix) const
{
    const Cluster *const pClusterI(clusterI.GetCluster());
    const Cluster *const pClusterJ(clusterJ.GetCluster());

    if (pClusterI == pClusterJ)
        return;

    // Target vertices should satisfy a minimum displacement
    const LArPointingCluster::Vertex &targetVertexI(useInnerI ? clusterI.GetInnerVertex() : clusterI.GetOuterVertex());
    const LArPointingCluster::Vertex &targetVertexJ(useInnerJ ? clusterJ.GetInnerVertex() : clusterJ.GetOuterVertex());

    const float distSquared_targetI_to_targetJ((targetVertexI.GetPosition() - targetVertexJ.GetPosition()).GetMagnitudeSquared());

    if( distSquared_targetI_to_targetJ > 400.f )
        return;

    // Target vertices should be the closest pair of vertices
    const LArPointingCluster::Vertex &oppositeVertexI(useInnerI ? clusterI.GetOuterVertex() : clusterI.GetInnerVertex());
    const LArPointingCluster::Vertex &oppositeVertexJ(useInnerJ ? clusterJ.GetOuterVertex() : clusterJ.GetInnerVertex());
    
    const float distSquared_targetI_to_oppositeJ((targetVertexI.GetPosition() - oppositeVertexJ.GetPosition()).GetMagnitudeSquared());
    const float distSquared_oppositeI_to_targetJ((oppositeVertexI.GetPosition() - targetVertexJ.GetPosition()).GetMagnitudeSquared());
    const float distSquared_oppositeI_to_oppositeJ((oppositeVertexI.GetPosition() - oppositeVertexJ.GetPosition()).GetMagnitudeSquared());

    if (distSquared_targetI_to_targetJ > std::min(distSquared_oppositeI_to_oppositeJ, std::min(distSquared_targetI_to_oppositeJ, distSquared_oppositeI_to_targetJ)))  
        return;

    if (distSquared_oppositeI_to_oppositeJ < std::max(distSquared_targetI_to_targetJ, std::max(distSquared_targetI_to_oppositeJ, distSquared_oppositeI_to_targetJ)))  
        return;

    // Check that new layer occupancy would be reasonable
    if (LArClusterHelper::GetLayerOccupancy(pClusterI, pClusterJ) < 0.75)
        return;

    // Check that vertices have a reasonable linear fit
    if (targetVertexI.GetRms() > 1.f || targetVertexJ.GetRms() > 1.f)
        return;
    

    // Association types
    LongitudinalAssociation::AssociationType associationType(LongitudinalAssociation::UNASSOCIATED);
    LongitudinalAssociation::StrengthType strengthType(LongitudinalAssociation::NONE);


    // Requirements for Nodes
    const CartesianVector &vertexPositionI(targetVertexI.GetPosition());
    const CartesianVector &vertexPositionJ(targetVertexJ.GetPosition());
    const CartesianVector &vertexDirectionI(targetVertexI.GetDirection());
    const CartesianVector &vertexDirectionJ(targetVertexJ.GetDirection());

    const float distanceSquared((vertexPositionI - vertexPositionJ).GetMagnitudeSquared());

    if (distanceSquared < 2.f * m_spatialResolution * m_spatialResolution)
    {
        associationType = LongitudinalAssociation::NODE;
        strengthType = LongitudinalAssociation::WEAK;

        if (distanceSquared < m_spatialResolution * m_spatialResolution)
        {
            const float cosTheta(-vertexDirectionI.GetDotProduct(vertexDirectionJ));

            if (cosTheta > 0.906f)
            {
                strengthType = LongitudinalAssociation::STANDARD;

                if (cosTheta > 0.966f)
                    strengthType = LongitudinalAssociation::STRONG;
            }
        }
    }


    // Requirements for Emissions
    const float clusterLengthI((targetVertexI.GetPosition() - oppositeVertexI.GetPosition()).GetMagnitude());
    const float clusterLengthJ((targetVertexJ.GetPosition() - oppositeVertexJ.GetPosition()).GetMagnitude());

    if ((clusterLengthI > 5.f) && (clusterLengthJ > 5.f))
    {
        // Calculate impact parameters
        if (strengthType < LongitudinalAssociation::STRONG)
        {
            const float cosTheta(-vertexDirectionI.GetDotProduct(vertexDirectionJ));
            const float cosThetaI((vertexPositionI - vertexPositionJ).GetUnitVector().GetDotProduct(vertexDirectionI));
            const float cosThetaJ((vertexPositionJ - vertexPositionI).GetUnitVector().GetDotProduct(vertexDirectionJ));

            float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
            LArVertexHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, vertexPositionJ, rL1, rT1);
            LArVertexHelper::GetImpactParameters(vertexPositionJ, vertexDirectionJ, vertexPositionI, rL2, rT2);

            if ((cosTheta > 0.985f) && (std::fabs(cosThetaI) > 0.25f) && (std::fabs(cosThetaJ) > 0.25f) && 
                (targetVertexI.GetRms() < 0.5f && targetVertexJ.GetRms() < 0.5f) &&
                (rL1 > -5.f && rL1 < +15.f && rL1 + clusterLengthJ > +15.f) && 
                (rL2 > -5.f && rL2 < +15.f && rL2 + clusterLengthI > +15.f) && 
                (rT1 < 2.f * m_spatialResolution) && (rT2 < 2.f * m_spatialResolution))
            {
                associationType = LongitudinalAssociation::EMISSION;
                strengthType = LongitudinalAssociation::STRONG;
            }
        }
    }

    if (strengthType > LongitudinalAssociation::NONE)
    {
        const LongitudinalAssociation::VertexType vertexTypeI(targetVertexI.IsInnerVertex() ? LongitudinalAssociation::INNER : LongitudinalAssociation::OUTER);
        const LongitudinalAssociation::VertexType vertexTypeJ(targetVertexJ.IsInnerVertex() ? LongitudinalAssociation::INNER : LongitudinalAssociation::OUTER);
        (void) longitudinalAssociationMatrix[pClusterI].insert(LongitudinalAssociationMap::value_type(pClusterJ, LongitudinalAssociation(vertexTypeI, vertexTypeJ, associationType, strengthType, clusterLengthJ)));
        (void) longitudinalAssociationMatrix[pClusterJ].insert(LongitudinalAssociationMap::value_type(pClusterI, LongitudinalAssociation(vertexTypeJ, vertexTypeI, associationType, strengthType, clusterLengthI)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalExtensionAlgorithm::AreClustersAssociated(Cluster *pCluster1, Cluster *pCluster2, const LongitudinalAssociationMatrix &longitudinalAssociationMatrix) const
{
    // Search for the pCluster1 <-> pCluster2 association
    LongitudinalAssociationMatrix::const_iterator iter1A = longitudinalAssociationMatrix.find(pCluster1);

    if (longitudinalAssociationMatrix.end() == iter1A)
        return false;

    const LongitudinalAssociationMap &longitudinalAssociationMap1(iter1A->second);

    LongitudinalAssociationMap::const_iterator iter1B = longitudinalAssociationMap1.find(pCluster2);

    if (longitudinalAssociationMap1.end() == iter1B)
        return false;

    const LongitudinalAssociation &longitudinalAssociationCandidate1(iter1B->second);

    if (longitudinalAssociationCandidate1.GetStrength() < LongitudinalAssociation::STANDARD)
        return false;

    // Bail out if there is a stronger association, or if this is a prong
    for (LongitudinalAssociationMap::const_iterator iter1C = longitudinalAssociationMap1.begin(), iter1CEnd = longitudinalAssociationMap1.end(); iter1C != iter1CEnd; ++iter1C)
    {
        const Cluster *pClusterAlternative1(iter1C->first);

        if (pClusterAlternative1 == pCluster2)
            continue;

        const LongitudinalAssociation &longitudinalAssociationAlternative1(iter1C->second);

        if (longitudinalAssociationCandidate1.GetParent() == longitudinalAssociationAlternative1.GetParent())
        {
            if (longitudinalAssociationAlternative1.GetAssociation() == LongitudinalAssociation::NODE)
                return false;

            if ((longitudinalAssociationAlternative1.GetStrength() > longitudinalAssociationCandidate1.GetStrength()) ||
                ((longitudinalAssociationAlternative1.GetStrength() == longitudinalAssociationCandidate1.GetStrength()) &&
                (longitudinalAssociationAlternative1.GetFigureOfMerit() > longitudinalAssociationCandidate1.GetFigureOfMerit())) )
            {
                return false;
            }
        }
    }

    // Search for the pCluster2 <-> pCluster1 association
    LongitudinalAssociationMatrix::const_iterator iter2A = longitudinalAssociationMatrix.find(pCluster2);

    if (longitudinalAssociationMatrix.end() == iter2A)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const LongitudinalAssociationMap &longitudinalAssociationMap2(iter2A->second);

    LongitudinalAssociationMap::const_iterator iter2B = longitudinalAssociationMap2.find(pCluster1);

    if (longitudinalAssociationMap2.end() == iter2B)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const LongitudinalAssociation &longitudinalAssociationCandidate2(iter2B->second);

    if (longitudinalAssociationCandidate2.GetStrength() < LongitudinalAssociation::STANDARD)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    // Bail out if there is a stronger association, or if this is a prong
    for (LongitudinalAssociationMap::const_iterator iter2C = longitudinalAssociationMap2.begin(), iter2CEnd = longitudinalAssociationMap2.end(); iter2C != iter2CEnd; ++iter2C)
    {
        const Cluster *pClusterAlternative2(iter2C->first);

        if (pClusterAlternative2 == pCluster1)
            continue;

        const LongitudinalAssociation &longitudinalAssociationAlternative2(iter2C->second);

        if (longitudinalAssociationCandidate2.GetParent() == longitudinalAssociationAlternative2.GetParent())
        {
            if (longitudinalAssociationAlternative2.GetAssociation() == LongitudinalAssociation::NODE)
                return false;

            if ((longitudinalAssociationAlternative2.GetStrength() > longitudinalAssociationCandidate2.GetStrength()) ||
                ((longitudinalAssociationAlternative2.GetStrength() == longitudinalAssociationCandidate2.GetStrength()) &&
                (longitudinalAssociationAlternative2.GetFigureOfMerit() > longitudinalAssociationCandidate2.GetFigureOfMerit())) )
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
