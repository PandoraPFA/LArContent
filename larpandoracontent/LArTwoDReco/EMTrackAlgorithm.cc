/**
 *  @file   EMTrackAlgorithm.cc
 *
 *  @brief  Implementation of the em track algorithm class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/EMTrackAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{


EMTrackAlgorithm::EMTrackAlgorithm() :
    m_caloHitListName(),
    m_caloHitToParentClusterMap(),
    m_minCaloHits(25),
    m_minSeparationDistance(27),
    m_maxXSeparation(10),
    m_maxZSeparation(10),
    m_slidingFitWindow(20)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EMTrackAlgorithm::ClusterAssociation::ClusterAssociation(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster, const CartesianVector &innerMergePoint, const CartesianVector &innerMergeDirection, const CartesianVector &outerMergePoint, const CartesianVector &outerMergeDirection, const CartesianVector &connectingLineDirection) :
    m_pInnerCluster(pInnerCluster),
    m_pOuterCluster(pOuterCluster),
    m_innerMergePoint(innerMergePoint),
    m_innerMergeDirection(innerMergeDirection),
    m_outerMergePoint(outerMergePoint),
    m_outerMergeDirection(outerMergeDirection),
    m_connectingLineDirection(connectingLineDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EMTrackAlgorithm::ClusterAssociation::ClusterAssociation() :
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
    
bool EMTrackAlgorithm::SortByDistanceToLine::operator() (const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs)
{
    const CartesianVector lhsDistanceVector(pLhs->GetPositionVector() - m_referencePoint);
    const CartesianVector rhsDistanceVector(pRhs->GetPositionVector() - m_referencePoint);

    const float lhsProjectedDistance(lhsDistanceVector.GetDotProduct(m_referenceDirection));
    const float rhsProjectedDistance(rhsDistanceVector.GetDotProduct(m_referenceDirection));

    return (lhsProjectedDistance < rhsProjectedDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EMTrackAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    TwoDSlidingFitResultMap microSlidingFitResultMap;
    TwoDSlidingFitResultMap macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    // loop through clusters, finding associations and merge until no further merges can be made
    bool mergeMade(true);
    while (mergeMade)
    {
        mergeMade = false;
        for (const Cluster *const innerCluster : clusterVector)
        {
            ClusterAssociation clusterAssociation;
            
            if(!this->FindBestClusterAssociation(innerCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap, clusterAssociation))
                continue;
                        
            // ATTN map will gain hanging pointers during the merging process, these are not accessed
            CaloHitToParentClusterMap caloHitToParentClusterMap;
            CaloHitVector extrapolatedCaloHitVector;
            this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, extrapolatedCaloHitVector, caloHitToParentClusterMap);

            if (extrapolatedCaloHitVector.empty())
                continue;
        
            if (!this->IsTrackContinuous(clusterAssociation, extrapolatedCaloHitVector))
                continue;

            this->RefineTracks(clusterAssociation, extrapolatedCaloHitVector, microSlidingFitResultMap);
            
            this->AddHitsToCluster(clusterAssociation, caloHitToParentClusterMap, extrapolatedCaloHitVector, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

            this->UpdateSlidingFitResultMap(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

            mergeMade = true;
            break;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector)
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHits)
            continue;

        clusterVector.push_back(pCluster);
    }
    
    // sort hits by Z, then X and then Y ISOBEL THIS DESCRIPTION ISN'T CORRECT!
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::InitialiseSlidingFitResultMaps(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap)
{
    microSlidingFitResultMap.clear();
    macroSlidingFitResultMap.clear();

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

void EMTrackAlgorithm::UpdateSlidingFitResultMap(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResultMap::const_iterator fitResultIter(microSlidingFitResultMap.find(pCluster));
        //ISOBEL - DO I GUARD AGAINST NOT BEING IN BOTH MAPS??
        
        if(fitResultIter == microSlidingFitResultMap.end())
        {
            try
            {
                (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
                (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, 1000000, slidingFitPitch)));
            }
            catch (StatusCodeException &) {}
        }
    }
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::FindBestClusterAssociation(const Cluster *const pCurrentCluster, const ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap, const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation)
{
    const ClusterVector::const_iterator currentIter(std::find(clusterVector.begin(), clusterVector.end(), pCurrentCluster));

    const TwoDSlidingFitResultMap::const_iterator currentMicroFitIter(microSlidingFitResultMap.find(pCurrentCluster));
    if (currentMicroFitIter == microSlidingFitResultMap.end())
        return false;

    const TwoDSlidingFitResultMap::const_iterator currentMacroFitIter(macroSlidingFitResultMap.find(pCurrentCluster));
    if (currentMacroFitIter == macroSlidingFitResultMap.end())
        return false;

    bool foundAssociation(false);
    unsigned int maxHits(0);    
    for (ClusterVector::const_iterator testIter = std::next(currentIter); testIter != clusterVector.end(); ++testIter)
    {
        const Cluster *const pTestCluster = *testIter;
        
        if (pTestCluster->GetNCaloHits() < maxHits)
            continue;

        const TwoDSlidingFitResultMap::const_iterator testMicroFitIter(microSlidingFitResultMap.find(pTestCluster));
        if (testMicroFitIter == microSlidingFitResultMap.end())
            continue;

        const TwoDSlidingFitResultMap::const_iterator testMacroFitIter(macroSlidingFitResultMap.find(pTestCluster));
        if (testMacroFitIter == macroSlidingFitResultMap.end())
            continue;
        
        CartesianVector currentMergePoint(0.f, 0.f, 0.f), testMergePoint(0.f, 0.f, 0.f), currentMergeDirection(0.f, 0.f, 0.f), testMergeDirection(0.f, 0.f, 0.f);
        if (!this->GetClusterMergingCoordinates(currentMicroFitIter->second, currentMacroFitIter->second, testMacroFitIter->second, currentMergePoint, currentMergeDirection, true)
            || !this->GetClusterMergingCoordinates(testMicroFitIter->second, testMacroFitIter->second, currentMacroFitIter->second, testMergePoint, testMergeDirection, false))
        {
            continue;
        }

        if (!AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection))
            continue;

        foundAssociation = true;
        maxHits = pTestCluster->GetNCaloHits();

        //ClusterList theCurrent, theTest;
        //theCurrent.push_back(pCurrentCluster);
        //theTest.push_back(pTestCluster);
        CartesianVector connectingLineDirection(testMergePoint.GetX() - currentMergePoint.GetX(), 0.f, testMergePoint.GetZ() - currentMergePoint.GetZ());
        if (LArClusterHelper::SortByPosition(pCurrentCluster, pTestCluster))
        {
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCurrent, "INNER", BLUE, 2);
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theTest, "OUTER", RED, 2);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            clusterAssociation = ClusterAssociation(pCurrentCluster, pTestCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection, connectingLineDirection.GetUnitVector());
        }
        else
        {
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCurrent, "OUTER", BLUE, 2);
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theTest, "INNER", RED, 2);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            clusterAssociation = ClusterAssociation(pTestCluster, pCurrentCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection, connectingLineDirection.GetUnitVector());
        }
        //theCurrent.clear();
        //theTest.clear();
        
    }
    
    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult,
    const TwoDSlidingFitResult &associatedMacroFitResult, CartesianVector &currentMergePosition, CartesianVector &currentMergeDirection, const bool isInner)
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
    
    for (unsigned int i = startLayer; i != loopTerminationLayer; i += isInner ? -1 : 1)
    {
        const auto microIter(currentMicroLayerFitResultMap.find(i));

        if (microIter == currentMicroLayerFitResultMap.end())
        {
            if (i == endLayer)
                return false;
                    
            continue;
        }
        
        CartesianVector microDirection(0.f, 0.f, 0.f);
        currentMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
    
        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(associatedAverageDirection));
        if (cosDirectionOpeningAngle > 0.9995)
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
        
        // Abort merging process if cannot find a stable region with gradient close to associated average 
        if (i == endLayer)
            return false;                                         
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::AreClustersAssociated(const CartesianVector &currentPoint, const CartesianVector &currentDirection, const CartesianVector &testPoint, const CartesianVector &testDirection)
{
    if (testPoint.GetZ() < currentPoint.GetZ())
        return false;

    // check that clusters are reasonably far away
    const float separationDistance(std::sqrt(currentPoint.GetDistanceSquared(testPoint)));
    if (separationDistance < m_minSeparationDistance)
        return false;
    
    // check that opening angle is not too large
    // ATTN - HAVE TURNED THE DIRECTION AROUND IN PREVIOUS FUNCTION
    if (currentDirection.GetCosOpeningAngle(testDirection * (-1.0)) < 0.99)
        return false;
    
    // check that fit allows you to get from one merge point to other merge point
    const CartesianVector extrapolatedCurrentPoint(currentPoint + (currentDirection * separationDistance));
    const CartesianVector extrapolatedTestPoint(testPoint + (testDirection * separationDistance));
    
    if ((extrapolatedCurrentPoint.GetX() > testPoint.GetX() + m_maxXSeparation) || (extrapolatedCurrentPoint.GetX() < testPoint.GetX() - m_maxXSeparation))
        return false;

    if ((extrapolatedCurrentPoint.GetZ() > testPoint.GetZ() + m_maxZSeparation) || (extrapolatedCurrentPoint.GetZ() < testPoint.GetZ() - m_maxZSeparation))
        return false;

    if ((extrapolatedTestPoint.GetX() > currentPoint.GetX() + m_maxXSeparation) || (extrapolatedTestPoint.GetX() < currentPoint.GetX() - m_maxXSeparation))
        return false;

    if ((extrapolatedTestPoint.GetZ() > currentPoint.GetZ() + m_maxZSeparation) || (extrapolatedTestPoint.GetZ() < currentPoint.GetZ() - m_maxZSeparation))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const ClusterList *const pClusterList, CaloHitVector &extrapolatedCaloHitVector, CaloHitToParentClusterMap &caloHitToParentClusterMap)
{
    const CartesianVector &innerPoint(clusterAssociation.GetInnerMergePoint());
    const CartesianVector &outerPoint(clusterAssociation.GetOuterMergePoint());
    const CartesianVector &connectingLineDirection(clusterAssociation.GetConnectingLineDirection());

    const float m_distanceFromLine(0.35);
    const float minX(std::min(innerPoint.GetX(), outerPoint.GetX()));
    const float maxX(std::max(innerPoint.GetX(), outerPoint.GetX()));
    const float minZ(std::min(innerPoint.GetZ(), outerPoint.GetZ()));
    const float maxZ(std::max(innerPoint.GetZ(), outerPoint.GetZ()));

    // ATTN: Consider hits from merging clusters as merging point may not be cluster end
    for (const Cluster *const pCluster : *pClusterList) 
    {
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second) 
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < minZ) || (hitPosition.GetZ() > maxZ))
                    continue;
                
                const float distanceFromLine(connectingLineDirection.GetCrossProduct(hitPosition - innerPoint).GetMagnitude());
                if (distanceFromLine > m_distanceFromLine)
                    continue;
                
                extrapolatedCaloHitVector.push_back(pCaloHit);
                caloHitToParentClusterMap[pCaloHit] = pCluster;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::IsTrackContinuous(const ClusterAssociation &clusterAssociation, CaloHitVector &extrapolatedCaloHitVector)
{
    const float m_maxTrackGaps(2.f);
    const float m_lineSegmentLength(7.f);
    
    const CartesianVector &innerMergePoint(clusterAssociation.GetInnerMergePoint()), &outerMergePoint(clusterAssociation.GetOuterMergePoint());
    const CartesianVector &trackDirection(clusterAssociation.GetConnectingLineDirection());
    const CartesianVector trackStep(trackDirection * m_lineSegmentLength);

    // To handle the final segment
    // Merge remainder segment into preceding segment and split into two if remainder segment length exceeded half of m_lineSegmentLength
    const float trackLength((outerMergePoint - innerMergePoint).GetMagnitude());
    const unsigned int fullSegments(floor(trackLength / m_lineSegmentLength));
    const float lengthOfTrackRemainder(trackLength - (fullSegments * m_lineSegmentLength));

    // sort hits by projected distance from the innerMergePoint
    std::sort(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), SortByDistanceToLine(innerMergePoint, trackDirection));

    CaloHitVector::const_iterator caloHitIter(extrapolatedCaloHitVector.begin());
    unsigned int hitsInSegment(0), segmentsWithoutHits(0);
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

        while (this->IsInLineSegment(lowerBoundary, upperBoundary, (*caloHitIter)->GetPositionVector()))
        {
            ++hitsInSegment;
            ++caloHitIter;

            if (caloHitIter == extrapolatedCaloHitVector.end())
                break;
        }

        // if number of hits in segment exceeds threshold then reset segmentsWithoutHits
        if (!hitsInSegment)
        {
            ++segmentsWithoutHits;
        }
        else
        {
            segmentsWithoutHits = 0;
        }

        if (segmentsWithoutHits > m_maxTrackGaps)
            return false;

        // case in which final two segments are merged need to leave the loop early
        if (i == (fullSegments - 1) && !(lengthOfTrackRemainder > m_lineSegmentLength * 0.5))
            return true;

        hitsInSegment = 0;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point)
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

void EMTrackAlgorithm::RefineTracks(const ClusterAssociation &clusterAssociation, const CaloHitVector &extrapolatedCaloHitVector,
    const TwoDSlidingFitResultMap &microFitResultMap)
{
    this->RefineTrack(clusterAssociation.GetInnerCluster(), clusterAssociation.GetInnerMergePoint(), extrapolatedCaloHitVector, microFitResultMap, true);
    this->RefineTrack(clusterAssociation.GetOuterCluster(), clusterAssociation.GetOuterMergePoint(), extrapolatedCaloHitVector, microFitResultMap, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::RefineTrack(const Cluster *const pCluster, const CartesianVector &splitPosition, const CaloHitVector &extrapolatedCaloHitVector,
    const TwoDSlidingFitResultMap &microFitResultMap, const bool isInner)
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

void EMTrackAlgorithm::AddHitsToCluster(const ClusterAssociation &clusterAssociation, const CaloHitToParentClusterMap &caloHitToParentClusterMap, const CaloHitVector &extrapolatedCaloHitVector, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap)
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

void EMTrackAlgorithm::RemoveClusterFromSlidingFitResultMaps(const Cluster *const pCluster, std::vector<TwoDSlidingFitResultMap*> &slidingFitResultMapVector)
{
    for (TwoDSlidingFitResultMap *const pSlidingFitResultMap : slidingFitResultMapVector)
    {
        const TwoDSlidingFitResultMap::const_iterator fitToDelete(pSlidingFitResultMap->find(pCluster));
        if (fitToDelete != pSlidingFitResultMap->end())
            pSlidingFitResultMap->erase(fitToDelete);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::RemoveClusterFromClusterVector(const Cluster *const pCluster, ClusterVector &clusterVector)
{
    ClusterVector::const_iterator clusterToDelete(std::find(clusterVector.begin(), clusterVector.end(), pCluster));
    if (clusterToDelete != clusterVector.end())
        clusterVector.erase(clusterToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode EMTrackAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListName", m_caloHitListName));

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

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content


    
