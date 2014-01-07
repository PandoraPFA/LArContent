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

void LongitudinalExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    LArPointingClusterList pointingClusterList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        pointingClusterList.push_back(LArPointingCluster(*iter));
    }

    ClusterAssociationMatrix intermediateAssociationMatrix;

    for (LArPointingClusterList::const_iterator iterI = pointingClusterList.begin(), iterEndI = pointingClusterList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster &clusterI = *iterI;

        for (LArPointingClusterList::const_iterator iterJ = iterI, iterEndJ = pointingClusterList.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster &clusterJ = *iterJ;

            if (clusterI.GetCluster() == clusterJ.GetCluster())
                continue;

            this->FillAssociationMatrix(clusterI, clusterJ, intermediateAssociationMatrix);
        }
    }

    this->FillReducedAssociationMatrix(intermediateAssociationMatrix, clusterAssociationMatrix);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ,
    ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    const bool useInnerI(true), useInnerJ(true);
    const bool useOuterI(false), useOuterJ(false);

    this->FillAssociationMatrix(clusterI, clusterJ, useInnerI, useInnerJ, clusterAssociationMatrix);
    this->FillAssociationMatrix(clusterI, clusterJ, useInnerI, useOuterJ, clusterAssociationMatrix);   
    this->FillAssociationMatrix(clusterI, clusterJ, useOuterI, useInnerJ, clusterAssociationMatrix);
    this->FillAssociationMatrix(clusterI, clusterJ, useOuterI, useOuterJ, clusterAssociationMatrix);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ, 
    const bool useInnerI, const bool useInnerJ, ClusterAssociationMatrix &clusterAssociationMatrix) const
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
    ClusterAssociation::AssociationType associationType(ClusterAssociation::NONE);


    // Requirements for Nodes
    const CartesianVector &vertexPositionI(targetVertexI.GetPosition());
    const CartesianVector &vertexPositionJ(targetVertexJ.GetPosition());
    const CartesianVector &vertexDirectionI(targetVertexI.GetDirection());
    const CartesianVector &vertexDirectionJ(targetVertexJ.GetDirection());

    const float distanceSquared((vertexPositionI - vertexPositionJ).GetMagnitudeSquared());

    if (distanceSquared < 2.f * m_spatialResolution * m_spatialResolution)
    {
        associationType = ClusterAssociation::WEAK;

        if (distanceSquared < m_spatialResolution * m_spatialResolution)
        {
            const float cosTheta(-vertexDirectionI.GetDotProduct(vertexDirectionJ));

            if (cosTheta > 0.906f)
            {
                associationType = ClusterAssociation::STRONG;
            }
        }
    }


    // Requirements for Emissions
    const float clusterLengthI((targetVertexI.GetPosition() - oppositeVertexI.GetPosition()).GetMagnitude());
    const float clusterLengthJ((targetVertexJ.GetPosition() - oppositeVertexJ.GetPosition()).GetMagnitude());

    if ((clusterLengthI > 5.f) && (clusterLengthJ > 5.f))
    {
        // Calculate impact parameters
        if (associationType < ClusterAssociation::STRONG)
        {
            const float cosTheta(-vertexDirectionI.GetDotProduct(vertexDirectionJ));
            const float cosThetaI((vertexPositionI - vertexPositionJ).GetUnitVector().GetDotProduct(vertexDirectionI));
            const float cosThetaJ((vertexPositionJ - vertexPositionI).GetUnitVector().GetDotProduct(vertexDirectionJ));

            float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
            LArVertexHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, vertexPositionJ, rL1, rT1);
            LArVertexHelper::GetImpactParameters(vertexPositionJ, vertexDirectionJ, vertexPositionI, rL2, rT2);

            if ((rL1 > -2.5f && rL1 < +15.f && rL1 < 2.f * clusterLengthJ) && 
                (rL2 > -2.5f && rL2 < +15.f && rL2 < 2.f * clusterLengthI) && 
		(rT1 < 3.f * m_spatialResolution) && (rT2 < 3.f * m_spatialResolution))
	    {
                associationType = ClusterAssociation::WEAK;

                if ((rT1 < 2.f * m_spatialResolution) && (rT2 < 2.f * m_spatialResolution) &&
                    (targetVertexI.GetRms() < 0.5f && targetVertexJ.GetRms() < 0.5f) &&
                    (cosTheta > 0.985f) && (std::fabs(cosThetaI) > 0.25f) && (std::fabs(cosThetaJ) > 0.25f))
                {
                    associationType = ClusterAssociation::STRONG;
		}
	    }
	}
    }

    if (associationType > ClusterAssociation::NONE)
    {
        const ClusterAssociation::VertexType vertexTypeI(targetVertexI.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
        const ClusterAssociation::VertexType vertexTypeJ(targetVertexJ.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
        (void) clusterAssociationMatrix[pClusterI].insert(ClusterAssociationMap::value_type(pClusterJ, ClusterAssociation(vertexTypeI, vertexTypeJ, associationType, clusterLengthJ)));
        (void) clusterAssociationMatrix[pClusterJ].insert(ClusterAssociationMap::value_type(pClusterI, ClusterAssociation(vertexTypeJ, vertexTypeI, associationType, clusterLengthI)));

// ---- BEGIN DISPLAY ----
// if(associationType == ClusterAssociation::STRONG) 
// std::cout << " --- STRONG --- " << std::endl; 
// else std::cout << " --- WEAK --- " << std::endl;
// ClusterList tempListI, tempListJ;
// tempListI.insert((Cluster*)pClusterI);
// tempListJ.insert((Cluster*)pClusterJ);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f); 
// PandoraMonitoringApi::VisualizeClusters(&tempListI, "ClusterI", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempListJ, "ClusterJ", GREEN);
// PandoraMonitoringApi::ViewEvent();
// ---- END DISPLAY ----
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillReducedAssociationMatrix(const ClusterAssociationMatrix &inputAssociationMatrix, ClusterAssociationMatrix &outputAssociationMatrix) const
{
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
                    association23.GetParent() == association21.GetParent())
		{
		    isAssociated = false;
                    break;
		}
	    }

            if (isAssociated)
	    {
                (void) outputAssociationMatrix[pCluster1].insert(ClusterAssociationMap::value_type(pCluster2, association12));
                (void) outputAssociationMatrix[pCluster2].insert(ClusterAssociationMap::value_type(pCluster1, association21));
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const
{
    ClusterAssociationMatrix intermediateAssociationMatrix;

    for (ClusterAssociationMatrix::const_iterator iter1 = clusterAssociationMatrix.begin(), iterEnd1 = clusterAssociationMatrix.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster* pParentCluster(iter1->first);
        const ClusterAssociationMap &clusterAssociationMap(iter1->second);

        Cluster* pBestClusterInner = NULL;
        ClusterAssociation bestAssociationInner(ClusterAssociation::INNER, ClusterAssociation::INNER,
                                                ClusterAssociation::NONE, 0.f);

        Cluster* pBestClusterOuter = NULL;
        ClusterAssociation bestAssociationOuter(ClusterAssociation::OUTER, ClusterAssociation::OUTER,
                                                ClusterAssociation::NONE, 0.f);

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

                    if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
		        pBestClusterInner = (Cluster*)pDaughterCluster;
                    else 
                        pBestClusterInner = NULL;
		}
	    }

            // Outer associations
            if (clusterAssociation.GetParent() == ClusterAssociation::OUTER)
	    {
	        if (clusterAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit())
	        {
		    bestAssociationOuter = clusterAssociation;

                    if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
		        pBestClusterOuter = (Cluster*)pDaughterCluster;
                    else 
                        pBestClusterOuter = NULL;
		}
	    }
	}

        if (pBestClusterInner)
	    (void) intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterInner, bestAssociationInner));
	
        if (pBestClusterOuter)
	    (void) intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterOuter, bestAssociationOuter));
    }


    for (ClusterAssociationMatrix::const_iterator iter3 = intermediateAssociationMatrix.begin(), iterEnd3 = intermediateAssociationMatrix.end(); iter3 != iterEnd3; ++iter3)
    {
        const Cluster* pParentCluster(iter3->first);
        const ClusterAssociationMap &parentAssociationMap(iter3->second);

        for (ClusterAssociationMap::const_iterator iter4 = parentAssociationMap.begin(), iterEnd4 = parentAssociationMap.end(); iter4 != iterEnd4; ++iter4)
	{
	    const Cluster* pDaughterCluster(iter4->first);
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
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f); 
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "ParentCluster", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "DaughterCluster", BLACK);
// PandoraMonitoringApi::ViewEvent();
// ---- END DISPLAY ----
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_spatialResolution = 1.25f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SpatialResolution", m_spatialResolution));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
