/**
 *  @file   TrackInEMShowerAlgorithm.cc
 *
 *  @brief  Implementation of the track in em shower algorithm class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/TrackInEMShowerAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackInEMShowerAlgorithm::ClusterAssociation::ClusterAssociation(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster,const CartesianVector &innerMergePoint,
        const CartesianVector &innerMergeDirection, const CartesianVector &outerMergePoint, const CartesianVector &outerMergeDirection) :
    m_pInnerCluster(pInnerCluster),
    m_pOuterCluster(pOuterCluster),
    m_innerMergePoint(innerMergePoint),
    m_innerMergeDirection(innerMergeDirection),
    m_outerMergePoint(outerMergePoint),
    m_outerMergeDirection(outerMergeDirection),
    m_connectingLineDirection(0.f, 0.f, 0.f)
{
    const CartesianVector connectingLineDirection(m_outerMergePoint.GetX() - m_innerMergePoint.GetX(), 0.f, m_outerMergePoint.GetZ() - m_innerMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackInEMShowerAlgorithm::ClusterAssociation::ClusterAssociation() :
    m_pInnerCluster(nullptr),
    m_pOuterCluster(nullptr),
    m_innerMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_innerMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_outerMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_outerMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_connectingLineDirection(CartesianVector(0.f, 0.f, 0.f))
{
}


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
    
TrackInEMShowerAlgorithm::TrackInEMShowerAlgorithm() :
    m_minCaloHits(25), 
    m_minSeparationDistance(27.f),
    m_maxXSeparation(5.f),
    m_maxZSeparation(5.f),
    m_slidingFitWindow(20),
    m_mergePointMinCosAngleDeviation(0.9995),
    m_minClusterLengthSum(75.f),
    m_minDirectionDeviationCosAngle(0.99),
    m_distanceFromLine(0.35),
    m_maxTrackGaps(2),
    m_lineSegmentLength(7.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool TrackInEMShowerAlgorithm::SortByDistanceToLine::operator() (const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs)
{
    const CartesianVector lhsDistanceVector(pLhs->GetPositionVector() - m_startPoint);
    const CartesianVector rhsDistanceVector(pRhs->GetPositionVector() - m_startPoint);

    const float lhsProjectedDistance(lhsDistanceVector.GetDotProduct(m_lineDirection));
    const float rhsProjectedDistance(rhsDistanceVector.GetDotProduct(m_lineDirection));

    return (lhsProjectedDistance < rhsProjectedDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackInEMShowerAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    TwoDSlidingFitResultMap microSlidingFitResultMap;
    TwoDSlidingFitResultMap macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    // Find longest length association and conitnue merging until no further merges can be made
    bool mergeMade(false);
    do
    {
        ClusterAssociation clusterAssociation;
        if(!this->FindBestClusterAssociation(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap, clusterAssociation))
            break;
                        
        CaloHitToParentClusterMap caloHitToParentClusterMap;
        CaloHitVector extrapolatedCaloHitVector;
        this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, extrapolatedCaloHitVector, caloHitToParentClusterMap);

        if (extrapolatedCaloHitVector.empty())
            break;
        
        if (!this->IsTrackContinuous(clusterAssociation, extrapolatedCaloHitVector))
            break;

        //////////////////////////
        ClusterList theInner, theOuter;
        theInner.push_back(clusterAssociation.GetInnerCluster()), theOuter.push_back(clusterAssociation.GetOuterCluster());
        CartesianVector innerMergePoint(clusterAssociation.GetInnerMergePoint()), outerMergePoint(clusterAssociation.GetOuterMergePoint());
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theInner, "INNER", RED);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theOuter, "OUTER", BLUE);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerMergePoint, "INNER MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerMergePoint, "OUTER MERGE POINT", BLUE, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        //////////////////////////
        
        this->RefineTracks(clusterAssociation, extrapolatedCaloHitVector, microSlidingFitResultMap);
        
        this->AddHitsToCluster(clusterAssociation, caloHitToParentClusterMap, extrapolatedCaloHitVector, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

        this->UpdateSlidingFitResultMap(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
        
        mergeMade = true;
    } while (mergeMade);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHits)
            continue;

        clusterVector.push_back(pCluster);
    }
    
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::InitialiseSlidingFitResultMaps(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
            (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, 1000000, slidingFitPitch)));
        }
        catch (StatusCodeException &) {}
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::UpdateSlidingFitResultMap(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResultMap::const_iterator microFitResultIter(microSlidingFitResultMap.find(pCluster));
        
        if(microFitResultIter == microSlidingFitResultMap.end())
        {
            try
            {
                (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
            }
            catch (StatusCodeException &) {}
        }

        const TwoDSlidingFitResultMap::const_iterator macroFitResultIter(macroSlidingFitResultMap.find(pCluster));
        
        if(macroFitResultIter == macroSlidingFitResultMap.end())
        {
            try
            {
                (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, 1000000, slidingFitPitch)));
            }
            catch (StatusCodeException &) {}
        }
    }
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::FindBestClusterAssociation(ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap,
    const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation) const
{
    bool foundAssociation(false);
    float maxLength(0.f);

    // Find the associated cluster pair with the longest length sum
    for (ClusterVector::const_iterator currentIter = clusterVector.begin(); currentIter != clusterVector.end(); ++currentIter)
    {
        const Cluster *const pCurrentCluster(*currentIter);
        
        const TwoDSlidingFitResultMap::const_iterator currentMicroFitIter(microSlidingFitResultMap.find(pCurrentCluster));
        if (currentMicroFitIter == microSlidingFitResultMap.end())
            return false;

        const TwoDSlidingFitResultMap::const_iterator currentMacroFitIter(macroSlidingFitResultMap.find(pCurrentCluster));
        if (currentMacroFitIter == macroSlidingFitResultMap.end())
            return false;
   
        for (ClusterVector::const_iterator testIter = std::next(currentIter); testIter != clusterVector.end(); ++testIter)
        {
            const Cluster *const pTestCluster = *testIter;

            const float lengthSum(LArClusterHelper::GetLength(pCurrentCluster) + LArClusterHelper::GetLength(pTestCluster)); 
            if ((lengthSum < maxLength) || (lengthSum < m_minClusterLengthSum))
                continue;

            const TwoDSlidingFitResultMap::const_iterator testMicroFitIter(microSlidingFitResultMap.find(pTestCluster));
            if (testMicroFitIter == microSlidingFitResultMap.end())
                continue;

            const TwoDSlidingFitResultMap::const_iterator testMacroFitIter(macroSlidingFitResultMap.find(pTestCluster));
            if (testMacroFitIter == macroSlidingFitResultMap.end())
                continue;
                        
            const bool isCurrentInner(LArClusterHelper::SortByPosition(pCurrentCluster, pTestCluster));

            //PandoraMonitoringApi::ViewEvent(this->GetPandora());

            
            CartesianVector currentMergePoint(0.f, 0.f, 0.f), testMergePoint(0.f, 0.f, 0.f), currentMergeDirection(0.f, 0.f, 0.f), testMergeDirection(0.f, 0.f, 0.f);
            if (!this->GetClusterMergingCoordinates(currentMicroFitIter->second, currentMacroFitIter->second, testMacroFitIter->second, currentMergePoint, currentMergeDirection, isCurrentInner) ||
                !this->GetClusterMergingCoordinates(testMicroFitIter->second, testMacroFitIter->second, currentMacroFitIter->second, testMergePoint, testMergeDirection, !isCurrentInner))
            {
                continue;
            }


            if (((currentMacroFitIter->second.GetGlobalMinLayerPosition().GetZ() > testMacroFitIter->second.GetGlobalMinLayerPosition().GetZ()) && (currentMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ() < testMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ())) || ((testMacroFitIter->second.GetGlobalMinLayerPosition().GetZ() > currentMacroFitIter->second.GetGlobalMinLayerPosition().GetZ()) && (testMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ() < currentMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ())))
           {
               std::cout << "CLUSTER INSIDE ANOTHER" << std::endl;
               ClusterList theCurrent, theTest;
                theCurrent.push_back(pCurrentCluster), theTest.push_back(pTestCluster);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCurrent, "CURRENT", BLUE);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theTest, "TEST", BLACK);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentMergePoint, "CURRENT POINT", BLUE, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testMergePoint, "TEST POINT", BLACK, 2);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                //continue;
           }


            

                


                

            if ((isCurrentInner && !AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection)) ||
                (!isCurrentInner && !AreClustersAssociated(testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection)))
            {
                continue;
            }

            foundAssociation = true;
            maxLength = lengthSum;

            if (isCurrentInner)
            {
                clusterAssociation = ClusterAssociation(pCurrentCluster, pTestCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection);
            }
            else
            {
                clusterAssociation = ClusterAssociation(pTestCluster, pCurrentCluster, testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection);
            }
        }
    }

    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult,
    const TwoDSlidingFitResult &associatedMacroFitResult, CartesianVector &currentMergePosition, CartesianVector &currentMergeDirection, const bool isInner) const
{
    CartesianVector currentAverageDirection(0.f, 0.f, 0.f);
    CartesianVector associatedAverageDirection(0.f, 0.f, 0.f);
    currentMacroFitResult.GetGlobalDirection(currentMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), currentAverageDirection);
    associatedMacroFitResult.GetGlobalDirection(associatedMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), associatedAverageDirection);    

    const LayerFitResultMap& currentMicroLayerFitResultMap(currentMicroFitResult.GetLayerFitResultMap());
    const unsigned int startLayer(isInner ? currentMicroFitResult.GetMaxLayer() : currentMicroFitResult.GetMinLayer());
    const unsigned int endLayer(isInner ? currentMicroFitResult.GetMinLayer() : currentMicroFitResult.GetMaxLayer());
    const unsigned int loopTerminationLayer(endLayer + (isInner ? -1 : 1));
    unsigned int goodPositionCount(0);
    unsigned int gradientStabilityHitWindow(std::ceil(currentMicroFitResult.GetCluster()->GetNCaloHits() * 0.1));

    //std::cout << "STABILITY WINDOW: " << gradientStabilityHitWindow << std::endl;
    
    for (unsigned int i = startLayer; i != loopTerminationLayer; i += isInner ? -1 : 1)
    {
        const auto microIter(currentMicroLayerFitResultMap.find(i));

        if (microIter == currentMicroLayerFitResultMap.end())
            continue;
        
        CartesianVector microDirection(0.f, 0.f, 0.f);
        currentMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
    
        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(associatedAverageDirection));
        //std::cout << cosDirectionOpeningAngle << std::endl;
        if (cosDirectionOpeningAngle > m_mergePointMinCosAngleDeviation)
        {
            if (goodPositionCount == 0)
            {
                // so that direction vectors face one another
                currentMergeDirection = currentAverageDirection * (isInner ? 1.f : -1.f);
                currentMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), currentMergePosition);
            }
            
            ++goodPositionCount;
        }
        else
        {
            goodPositionCount = 0;
        }

        if (goodPositionCount > gradientStabilityHitWindow)
            break;
        
        // abort merging process if cannot find a stable region with gradient close to associated average 
        if (i == endLayer)
            return false;                                         
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::AreClustersAssociated(const CartesianVector &innerPoint, const CartesianVector &innerDirection, const CartesianVector &outerPoint,
    const CartesianVector &outerDirection) const
{
    if (outerPoint.GetZ() < innerPoint.GetZ())
    {
        std::cout << "Z" << std::endl;
        return false;
    }

    // check that clusters are reasonably far away
    const float separationDistance(std::sqrt(innerPoint.GetDistanceSquared(outerPoint)));
    if (separationDistance < m_minSeparationDistance)
        return false;
    
    // check that opening angle is not too large
    // ATTN: cluster directions are pointing towards one another
    if (innerDirection.GetCosOpeningAngle(outerDirection * (-1.0)) < m_minDirectionDeviationCosAngle)
    {
        std::cout << "angle" << std::endl;
        return false;
    }
    
    // check that fit allows you to get from one merge point to other merge point
    const CartesianVector extrapolatedInnerPoint(innerPoint + (innerDirection * separationDistance));
    const CartesianVector extrapolatedOuterPoint(outerPoint + (outerDirection * separationDistance));

    if ((extrapolatedInnerPoint.GetX() > outerPoint.GetX() + m_maxXSeparation) || (extrapolatedInnerPoint.GetX() < outerPoint.GetX() - m_maxXSeparation))
        return false;

    if ((extrapolatedInnerPoint.GetZ() > outerPoint.GetZ() + m_maxZSeparation) || (extrapolatedInnerPoint.GetZ() < outerPoint.GetZ() - m_maxZSeparation))
        return false;

    if ((extrapolatedOuterPoint.GetX() > innerPoint.GetX() + m_maxXSeparation) || (extrapolatedOuterPoint.GetX() < innerPoint.GetX() - m_maxXSeparation))
        return false;

    if ((extrapolatedOuterPoint.GetZ() > innerPoint.GetZ() + m_maxZSeparation) || (extrapolatedOuterPoint.GetZ() < innerPoint.GetZ() - m_maxZSeparation))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const ClusterList *const pClusterList,
    CaloHitVector &extrapolatedCaloHitVector, CaloHitToParentClusterMap &caloHitToParentClusterMap) const
{
    const CartesianVector &innerPoint(clusterAssociation.GetInnerMergePoint()), &outerPoint(clusterAssociation.GetOuterMergePoint());
    const float minX(std::min(innerPoint.GetX(), outerPoint.GetX())), maxX(std::max(innerPoint.GetX(), outerPoint.GetX()));

    // ATTN: consider hits from merging clusters as merging point may not be cluster end
    for (const Cluster *const pCluster : *pClusterList) 
    {
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second) 
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < innerPoint.GetZ()) || (hitPosition.GetZ() > outerPoint.GetZ()))
                    continue;
                
                const float distanceFromLine(clusterAssociation.GetConnectingLineDirection().GetCrossProduct(hitPosition - innerPoint).GetMagnitude());
                
                if (distanceFromLine > m_distanceFromLine)
                    continue;
                
                extrapolatedCaloHitVector.push_back(pCaloHit);
                caloHitToParentClusterMap[pCaloHit] = pCluster;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsTrackContinuous(const ClusterAssociation &clusterAssociation, CaloHitVector &extrapolatedCaloHitVector) const
{
    const CartesianVector &innerMergePoint(clusterAssociation.GetInnerMergePoint()), &outerMergePoint(clusterAssociation.GetOuterMergePoint());
    const CartesianVector &trackDirection(clusterAssociation.GetConnectingLineDirection());
    const CartesianVector trackStep(trackDirection * m_lineSegmentLength);

    // to handle final segment - merge track remainder with preceding segment and if track remainder is more than half trackstep split into two
    const float trackLength((outerMergePoint - innerMergePoint).GetMagnitude());
    const unsigned int fullSegments(floor(trackLength / m_lineSegmentLength));
    const float lengthOfTrackRemainder(trackLength - (fullSegments * m_lineSegmentLength));

    // sort hits by projected distance from the innerMergePoint
    std::sort(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), SortByDistanceToLine(innerMergePoint, trackDirection));

    unsigned int segmentsWithoutHits(0);
    CaloHitVector::const_iterator caloHitIter(extrapolatedCaloHitVector.begin());

    for (unsigned int i = 0; i < (fullSegments + 1); ++i)
    {
        // if have run out of hits
        if (caloHitIter == extrapolatedCaloHitVector.end())
        {
            ++segmentsWithoutHits;
            if (segmentsWithoutHits > m_maxTrackGaps)
                return false;
            continue;
        }
        
        CartesianVector lowerBoundary(innerMergePoint + (trackStep * static_cast<float>(i)));
        CartesianVector upperBoundary(innerMergePoint + (trackStep * static_cast<float>(i + 1.f)));

        if (i >= fullSegments - 1)
        {
            if (lengthOfTrackRemainder > m_lineSegmentLength * 0.5f)
            {
                lowerBoundary = lowerBoundary - (trackStep * 0.5f * (i - fullSegments + 1.f)) + (trackDirection * 0.5f * (i - fullSegments + 1.f) * lengthOfTrackRemainder);
                upperBoundary = upperBoundary - (trackStep * 0.5f * (i - fullSegments + 2.f)) + (trackDirection * 0.5f * (i - fullSegments + 2.f) * lengthOfTrackRemainder);
            }
            else
            {
                upperBoundary = outerMergePoint;
            }
        }

        unsigned int hitsInSegment(0);
        while (this->IsInLineSegment(lowerBoundary, upperBoundary, (*caloHitIter)->GetPositionVector()))
        {
            ++hitsInSegment;
            ++caloHitIter;

            if (caloHitIter == extrapolatedCaloHitVector.end())
                break;
        }

        segmentsWithoutHits = hitsInSegment ? 0 : segmentsWithoutHits + 1;

        if (segmentsWithoutHits > m_maxTrackGaps)
            return false;

        // leave loop early if final merged segment is not split into two
        if (i == (fullSegments - 1) && !(lengthOfTrackRemainder > m_lineSegmentLength * 0.5))
            return true;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
{
    const float gradient = (-1.0)*(upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ() + gradient*upperBoundary.GetX())/gradient);
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ() + gradient*lowerBoundary.GetX())/gradient);

    const CartesianVector upper(xPointOnUpperLine, 0.f, point.GetZ());
    const CartesianVector lower(xPointOnLowerLine, 0.f, point.GetZ());

    if ((point.GetX() > xPointOnUpperLine) && (point.GetX() > xPointOnLowerLine))
        return false;

    if ((point.GetX() < xPointOnUpperLine) && (point.GetX() < xPointOnLowerLine))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::RefineTracks(const ClusterAssociation &clusterAssociation, const CaloHitVector &extrapolatedCaloHitVector,
    const TwoDSlidingFitResultMap &microFitResultMap) const
{
    this->RefineTrack(clusterAssociation.GetInnerCluster(), clusterAssociation.GetInnerMergePoint(), extrapolatedCaloHitVector, microFitResultMap, true);
    this->RefineTrack(clusterAssociation.GetOuterCluster(), clusterAssociation.GetOuterMergePoint(), extrapolatedCaloHitVector, microFitResultMap, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::RefineTrack(const Cluster *const pCluster, const CartesianVector &splitPosition, const CaloHitVector &extrapolatedCaloHitVector,
    const TwoDSlidingFitResultMap &microFitResultMap, const bool isInner) const
{
    const TwoDSlidingFitResult microFitResult(microFitResultMap.at(pCluster));
    float rL(0.f), rT(0.f);
        
    microFitResult.GetLocalPosition(splitPosition, rL, rT);

    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const CartesianVector hitPosition(pCaloHit->GetPositionVector());
            float thisL(0.f), thisT(0.f);

            microFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            if ((thisL > rL && isInner) || ((thisL < rL && !isInner)))
            {
                if (std::find(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), pCaloHit) != extrapolatedCaloHitVector.end())
                    continue;

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::AddHitsToCluster(const ClusterAssociation &clusterAssociation, const CaloHitToParentClusterMap &caloHitToParentClusterMap,
    const CaloHitVector &extrapolatedCaloHitVector, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    const Cluster *const pClusterToEnlarge(clusterAssociation.GetInnerCluster());
    const Cluster *const pClusterToDelete(clusterAssociation.GetOuterCluster());
    
    // remove from sliding fit map clusters that will be deleted or whose constituents will change
    std::vector<TwoDSlidingFitResultMap*> slidingFitResultVector({&microSlidingFitResultMap, &macroSlidingFitResultMap});
    this->RemoveClusterFromSlidingFitResultMaps(pClusterToEnlarge, slidingFitResultVector);
    this->RemoveClusterFromSlidingFitResultMaps(pClusterToDelete, slidingFitResultVector);
    
    for (const CaloHit *const pCaloHit : extrapolatedCaloHitVector)
    {
        CaloHitToParentClusterMap::const_iterator caloHitParentIter(caloHitToParentClusterMap.find(pCaloHit));

        if (caloHitParentIter == caloHitToParentClusterMap.end())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        if ((caloHitParentIter->second == pClusterToEnlarge) || (caloHitParentIter->second == pClusterToDelete))
            continue;

        const StatusCode statusCode(PandoraContentApi::RemoveFromCluster(*this, caloHitParentIter->second, pCaloHit));
            
        if (statusCode == STATUS_CODE_SUCCESS)
        {
            RemoveClusterFromSlidingFitResultMaps(caloHitParentIter->second, slidingFitResultVector);
            PandoraContentApi::AddToCluster(*this, pClusterToEnlarge, pCaloHit);
        }
        else if (statusCode == STATUS_CODE_NOT_ALLOWED)
        {
            RemoveClusterFromClusterVector(caloHitParentIter->second, clusterVector);
            PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, caloHitParentIter->second);
        }
        else
        {
            throw statusCode;
        }
    }

    this->RemoveClusterFromClusterVector(pClusterToDelete, clusterVector);
    PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pClusterToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::RemoveClusterFromSlidingFitResultMaps(const Cluster *const pCluster, std::vector<TwoDSlidingFitResultMap*> &slidingFitResultMapVector) const
{
    for (TwoDSlidingFitResultMap *const pSlidingFitResultMap : slidingFitResultMapVector)
    {
        const TwoDSlidingFitResultMap::const_iterator fitToDelete(pSlidingFitResultMap->find(pCluster));
        if (fitToDelete != pSlidingFitResultMap->end())
            pSlidingFitResultMap->erase(fitToDelete);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::RemoveClusterFromClusterVector(const Cluster *const pCluster, ClusterVector &clusterVector) const
{
    ClusterVector::const_iterator clusterToDelete(std::find(clusterVector.begin(), clusterVector.end(), pCluster));
    if (clusterToDelete != clusterVector.end())
        clusterVector.erase(clusterToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode TrackInEMShowerAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZSeparation", m_maxZSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparationDistance", m_minSeparationDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergePointMinCosAngleDeviation", m_mergePointMinCosAngleDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLengthSum", m_minClusterLengthSum));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDirectionDeviationCosAngle", m_minDirectionDeviationCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFromLine", m_distanceFromLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrackGaps", m_maxTrackGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LineSegmentLength", m_lineSegmentLength));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content


    
