/**
 *  @file   LArContent/src/LArClusterAssociation/LongitudinalConsolidationAlgorithm.cc
 * 
 *  @brief  Implementation of the transverse extension algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/LongitudinalConsolidationAlgorithm.h"

#include "LArHelpers/LArVertexHelper.h"
#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{
  
void LongitudinalConsolidationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        clusterVector.push_back(*iter);

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void LongitudinalConsolidationAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{   
    ClusterMergeMap intermediateMergeMap;

    // Loop over long clusters (possible parents)
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        const Cluster* pParentCluster = *iterI;

        if (LArClusterHelper::GetLengthSquared(pParentCluster) < m_longClusterLengthCut * m_longClusterLengthCut)
	    continue;

        LArClusterHelper::TwoDSlidingFitResult slidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(pParentCluster, 20, slidingFitResult);

        // Loop over short clusters (possible daughters)
        for (ClusterVector::const_iterator iterJ = clusterVector.begin(), iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
	    const Cluster* pDaughterCluster = *iterJ;
    
            if (pParentCluster == pDaughterCluster)
	        continue;

            if (LArClusterHelper::GetLengthSquared(pDaughterCluster) > m_shortClusterLengthCut * m_shortClusterLengthCut)
	        continue;

            // Try to fit the short clusters into gaps in the long clusters
            if (this->IsAssociated(slidingFitResult, pDaughterCluster))
	    {

// ClusterList tempList1, tempList2;
// tempList1.insert((Cluster*)pParentCluster);
// tempList2.insert((Cluster*)pDaughterCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "ParentCluster", RED));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "DaughterCluster", BLUE));
// PANDORA_MONITORING_API(ViewEvent());

                intermediateMergeMap[pDaughterCluster].insert((Cluster*)pParentCluster);
	    }
	}
    }
  
    // Merge daughter clusters with their longest parent cluster
    for (ClusterMergeMap::const_iterator iter1 = intermediateMergeMap.begin(), iterEnd1 = intermediateMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster* pDaughterCluster = iter1->first;

        Cluster* pParentCluster = NULL;
        float bestFigureOfMerit(0.f);

        for (ClusterList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster *pCandidateCluster = *iter2;
            const float figureOfMerit(LArClusterHelper::GetLengthSquared(pCandidateCluster));

            if (figureOfMerit > bestFigureOfMerit)
	    {  
                bestFigureOfMerit = figureOfMerit;
	        pParentCluster = (Cluster*)pCandidateCluster;
	    }
	}

        if (pParentCluster)
	{
            clusterMergeMap[pParentCluster].insert((Cluster*)pDaughterCluster);
            clusterMergeMap[pDaughterCluster].insert((Cluster*)pParentCluster);
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalConsolidationAlgorithm::IsAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const Cluster *const pCluster) const
{
    // Project each end of the candidate cluster onto the sliding fit and require
    // a sufficiently small displacement

    bool isFirstAssociated(false), isSecondAssociated(false);
    float firstL(0.f), firstT(0.f), secondT(0.f), secondL(0.f); 

    CartesianVector firstCoordinate(0.f, 0.f, 0.f);
    CartesianVector secondCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster, firstCoordinate, secondCoordinate);

    const float clusterLengthSquared(LArClusterHelper::GetLengthSquared(pCluster));

    try{
        CartesianVector projectedFirstCoordinate(0.f, 0.f, 0.f);
        slidingFitResult.GetGlobalFitProjection(firstCoordinate, projectedFirstCoordinate);
        slidingFitResult.GetLocalPosition(projectedFirstCoordinate, firstL, firstT);

        const float firstDisplacementSquared((projectedFirstCoordinate - firstCoordinate).GetMagnitudeSquared());
        isFirstAssociated = (firstDisplacementSquared < m_maxClusterDisplacementSquared && 
                             firstDisplacementSquared < clusterLengthSquared + 0.5 * m_maxClusterDisplacementSquared);
    }
    catch (StatusCodeException &)
    {
    }

    try{
        CartesianVector projectedSecondCoordinate(0.f, 0.f, 0.f);
        slidingFitResult.GetGlobalFitProjection(secondCoordinate, projectedSecondCoordinate);
        slidingFitResult.GetLocalPosition(projectedSecondCoordinate, secondL, secondT);

        const float secondDisplacementSquared((projectedSecondCoordinate - secondCoordinate).GetMagnitudeSquared());
        isSecondAssociated = (secondDisplacementSquared < m_maxClusterDisplacementSquared &&
                              secondDisplacementSquared < clusterLengthSquared + 0.5 * m_maxClusterDisplacementSquared);
    }
    catch (StatusCodeException &)
    {
    }

    if (!isFirstAssociated || !isSecondAssociated)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_longClusterLengthCut = 5.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LongClusterLengthCut", m_longClusterLengthCut));

    m_shortClusterLengthCut = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShortClusterLengthCut", m_shortClusterLengthCut));

    float m_maxClusterDisplacement = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterDisplacement", m_maxClusterDisplacement));
    m_maxClusterDisplacementSquared = m_maxClusterDisplacement * m_maxClusterDisplacement;

    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
