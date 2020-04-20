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
    m_slidingFitWindow(20),
    m_limitZ(false),
    m_useOtherCluster(false),
    m_abortIfNoPosition(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EMTrackAlgorithm::ClusterAssociation::ClusterAssociation(const Cluster *const pAssociatedCluster, const CartesianVector &innerMergePoint, const CartesianVector &innerMergeDirection, const CartesianVector &outerMergePoint, const CartesianVector &outerMergeDirection, const CartesianVector &connectingLineDirection) :
    m_pAssociatedCluster(pAssociatedCluster),
    m_innerMergePoint(innerMergePoint),
    m_innerMergeDirection(innerMergeDirection),
    m_outerMergePoint(outerMergePoint),
    m_outerMergeDirection(outerMergeDirection),
    m_connectingLineDirection(connectingLineDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EMTrackAlgorithm::ClusterAssociation::ClusterAssociation() :
    m_pAssociatedCluster(nullptr),
    m_innerMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_innerMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_outerMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_outerMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_connectingLineDirection(CartesianVector(0.f, 0.f, 0.f))
{
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
    this->InitialiseSlidingFitResultMap(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
    
    //PandoraMonitoringApi::Create(this->GetPandora());
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    /*
    ClusterList clusterList(clusterVector.begin(), clusterVector.end());
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &clusterList, "Selected Clusters", BLACK);
    
    for (const Cluster *const pCluster : clusterVector)
    {
        const CartesianVector position(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POINT", RED, 2);
        PandoraMonitoringApi::Pause(this->GetPandora());
    }
    
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    
    bool mergeMade(true);

    while (mergeMade)
    {
        mergeMade = false;
        for (const Cluster *const innerCluster : clusterVector)
        {

            ClusterAssociation clusterAssociation;
            
            // if no association is found, skip to the next cluster
            if(!this->FindBestClusterAssociation(innerCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap, clusterAssociation))
                continue;
            
            ClusterList currentCluster;
            currentCluster.push_back(innerCluster);

            ClusterList associatedCluster;
            associatedCluster.push_back(clusterAssociation.GetAssociatedCluster());

            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "CURRENT", BLACK);
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &associatedCluster, "ASSOCIATED", BLUE);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            
            
            // ATTN map will gain hanging pointers during the merging process, these are not accessed
            CaloHitToParentClusterMap caloHitToParentClusterMap;
            CaloHitVector extrapolatedCaloHitVector;
            this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, extrapolatedCaloHitVector, caloHitToParentClusterMap);

            if (extrapolatedCaloHitVector.empty())
                continue;
        
            if (!this->IsTrackContinuous(clusterAssociation, extrapolatedCaloHitVector))
                continue;

            //std::cout << "HERE" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());


            this->RemoveClusteringErrors(innerCluster, clusterAssociation, extrapolatedCaloHitVector, microSlidingFitResultMap);
            
            this->AddHitsToCluster(innerCluster, clusterAssociation.GetAssociatedCluster(), caloHitToParentClusterMap, extrapolatedCaloHitVector, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

            /////////////////////////////
            //MONITORING PURPOSES
            //ClusterList theCluster;
            //theCluster.push_back(innerCluster);
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "AFTER", GREEN);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            /////////////////////////////

            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), pClusterList, "ALL", BLUE);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());


            this->UpdateSlidingFitResultMap(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

            mergeMade = true;
            break;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

//Call this REFINE TRACK
void EMTrackAlgorithm::RemoveClusteringErrors(const Cluster *const innerCluster, const ClusterAssociation &clusterAssociation, const CaloHitVector &extrapolatedCaloHitVector, const TwoDSlidingFitResultMap &microFitResultMap)
{
    ClusterList mergingClusters({innerCluster, clusterAssociation.GetAssociatedCluster()});

    for (const Cluster *const pCluster : mergingClusters)
    {
        const bool isInner(pCluster == innerCluster);
        const TwoDSlidingFitResult microFitResult(microFitResultMap.at(pCluster));
        const CartesianVector &splitPosition(isInner ? clusterAssociation.GetInnerMergePoint() : clusterAssociation.GetOuterMergePoint());
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
                    {
                        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "CORRECT TO MERGE HIT", SPRING, 2);
                        continue;
                    }
                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "ERROR HIT", RED, 2);
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));
                }
            }
        }
    }

    
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector)
{

    //std::cout << "MIN HITS: " << m_minCaloHits << std::endl;
    
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHits)
            continue;

        clusterVector.push_back(pCluster);
    }
    
    // sort hits by Z, then X and then Y
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByPosition);
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

void EMTrackAlgorithm::RemoveClusterFromClusterVector(const Cluster *const pCluster, ClusterVector &clusterVector)
{
    ClusterVector::const_iterator clusterToDelete(std::find(clusterVector.begin(), clusterVector.end(), pCluster));
    if (clusterToDelete != clusterVector.end())
        clusterVector.erase(clusterToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::RemoveClusterFromSlidingFitResultMap(const Cluster *const pCluster, TwoDSlidingFitResultMap &slidingFitResultMap)
{
    const TwoDSlidingFitResultMap::const_iterator fitToDelete(slidingFitResultMap.find(pCluster));
    if (fitToDelete != slidingFitResultMap.end())
        slidingFitResultMap.erase(fitToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::AddHitsToCluster(const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete, const CaloHitToParentClusterMap &caloHitToParentClusterMap, const CaloHitVector &extrapolatedCaloHitVector, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap)
{
    // remove from sliding fit map clusters that will be deleted or whose constituents will change
    this->RemoveClusterFromSlidingFitResultMap(pClusterToEnlarge, microSlidingFitResultMap);
    this->RemoveClusterFromSlidingFitResultMap(pClusterToEnlarge, macroSlidingFitResultMap);
    this->RemoveClusterFromSlidingFitResultMap(pClusterToDelete, microSlidingFitResultMap);
    this->RemoveClusterFromSlidingFitResultMap(pClusterToDelete, macroSlidingFitResultMap);
    
    for (const CaloHit *const pCaloHit : extrapolatedCaloHitVector)
    {
        CaloHitToParentClusterMap::const_iterator caloHitParentIter(caloHitToParentClusterMap.find(pCaloHit));

        if (caloHitParentIter == caloHitToParentClusterMap.end())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        if (caloHitParentIter->second == pClusterToEnlarge || caloHitParentIter->second == pClusterToDelete)
            continue;

        const StatusCode statusCode(PandoraContentApi::RemoveFromCluster(*this, caloHitParentIter->second, pCaloHit));
            
        if (statusCode == STATUS_CODE_SUCCESS)
        {
            RemoveClusterFromSlidingFitResultMap(caloHitParentIter->second, microSlidingFitResultMap);
            RemoveClusterFromSlidingFitResultMap(caloHitParentIter->second, macroSlidingFitResultMap);
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

//-------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::IsTrackContinuous(const ClusterAssociation &clusterAssociation, CaloHitVector &extrapolatedCaloHitVector)
{
    const float m_maxTrackGaps(2);
    const float m_lineSegmentLength(7);
    
    const CartesianVector &innerMergePoint(clusterAssociation.GetInnerMergePoint()), &outerMergePoint(clusterAssociation.GetOuterMergePoint());
    const CartesianVector &trackDirection(clusterAssociation.GetConnectingLineDirection());
    const CartesianVector trackStep(trackDirection * m_lineSegmentLength);

    // sort hits by project distance from the innerMergePoint
    std::sort(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), SortByDistanceToLine(innerMergePoint, trackDirection));

    unsigned int hitsInSegment(0), segmentsWithoutHits(0);

    CaloHitVector::const_iterator caloHitIter(extrapolatedCaloHitVector.begin());

    float trackLength((outerMergePoint - innerMergePoint).GetMagnitude());
    unsigned int fullSegments(floor(trackLength / m_lineSegmentLength));
    float lengthOfTrackRemainder(trackLength - fullSegments * m_lineSegmentLength);

    //std::cout << "TRACK LENGTH: " << trackLength << std::endl;
    //std::cout << "LINE SEGMENT LENGTH: " << m_lineSegmentLength << std::endl;
    //std::cout << "FULL SEGMENTS: " << fullSegments << std::endl;
    //std::cout << "LENGTH OF TRACK REMAINDER: " << lengthOfTrackRemainder << std::endl;

    //std::cout << "TRACK STEP: " << trackStep << std::endl;
    
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
        
        CartesianVector lowerBoundary(innerMergePoint + trackStep * static_cast<float>(i));
        CartesianVector upperBoundary(innerMergePoint + trackStep * static_cast<float>(i + 1.0));

        // merge remainder segment into preceding segment, split into two if remainder segment length is over half of the m_lineSegmentLength
        if (i == fullSegments - 1)
        {
            if (lengthOfTrackRemainder > m_lineSegmentLength * 0.5)
            {
                upperBoundary = lowerBoundary + (trackStep * 0.5) + (trackDirection * lengthOfTrackRemainder * 0.5);
            }
            else
            {
                upperBoundary = outerMergePoint;
            }
        }

        if (i == fullSegments)
        {
            lowerBoundary = outerMergePoint - (trackStep * 0.5) - (trackDirection * lengthOfTrackRemainder * 0.5);
            upperBoundary = outerMergePoint;
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

        //std::cout << "HITS IN SEGMENT: " << hitsInSegment << std::endl;
        //std::cout << "SEGMENTS WITHOUT HITS: " << segmentsWithoutHits << std::endl;

        if (segmentsWithoutHits > m_maxTrackGaps)
            return false;

        // case in which final two segments are merged, therefore need to leave the loop early
        if (i == (fullSegments - 1) && !(lengthOfTrackRemainder > m_lineSegmentLength * 0.5))
            return true;

        hitsInSegment = 0;
    }

    return true;
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

bool EMTrackAlgorithm::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point)
{

    const float gradient = (-1.0)*(upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ() + gradient*upperBoundary.GetX())/gradient);
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ() + gradient*lowerBoundary.GetX())/gradient);

    const CartesianVector upper(xPointOnUpperLine, 0.f, point.GetZ());
    const CartesianVector lower(xPointOnLowerLine, 0.f, point.GetZ());

    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "POINT", GREEN, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &lowerBoundary, "LOWER BOUNDARY", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upperBoundary, "UPPER BOUNDARY", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upper, "UPPER", VIOLET, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &lower, "LOWER", VIOLET, 2);
    */
    if (point.GetX() > xPointOnUpperLine && point.GetX() > xPointOnLowerLine)
    {
        //std::cout << "failed higher than x" << std::endl;
        return false;
    }

    if (point.GetX() < xPointOnUpperLine && point.GetX() < xPointOnLowerLine)
    {
        //std::cout << "failer lower than x" << std::endl;
        return false;
    }

    //std::cout << "passed" << std::endl;

    //PandoraMonitoringApi::Pause(this->GetPandora());
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const ClusterList *const pClusterList, CaloHitVector &extrapolatedCaloHitVector, CaloHitToParentClusterMap &caloHitToParentClusterMap)
{

    const CartesianVector &innerPoint(clusterAssociation.GetInnerMergePoint());
    const CartesianVector &outerPoint(clusterAssociation.GetOuterMergePoint());
    const CartesianVector &connectingLineDirection(clusterAssociation.GetConnectingLineDirection());
    
    /////////////////////////////
    //MONITORING PURPOSES
    /*
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &innerPoint, &outerPoint, "CONNECTING LINE", BLUE, 2, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerPoint, "INNER", BLUE, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerPoint, "OUTER", BLUE, 2);
    
    ClusterList innerClusterList;
    innerClusterList.push_back(innerCluster);
    ClusterList outerClusterList;
    outerClusterList.push_back(clusterAssociation.GetAssociatedCluster());

    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &innerClusterList, "INNER CLUSTER", BLUE);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &outerClusterList, "OUTER CLUSTER", BLUE);
    */
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    /////////////////////////////

    float m_distanceFromLine(0.35);
    float minX(std::min(innerPoint.GetX(), outerPoint.GetX()));
    float maxX(std::max(innerPoint.GetX(), outerPoint.GetX()));
    float minZ(std::min(innerPoint.GetZ(), outerPoint.GetZ()));
    float maxZ(std::max(innerPoint.GetZ(), outerPoint.GetZ()));
    /*
    CartesianVector TR(maxX, 0.f, maxZ);
    CartesianVector TL(minX, 0.f, maxZ);
    CartesianVector BR(maxX, 0.f, minZ);
    CartesianVector BL(minX, 0.f, minZ);
    
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &TR, &TL, "LINE", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &TR, &BR, "LINE", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &TL, &BL, "LINE", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &BR, &BL, "LINE", RED, 2, 2);
    */
    for (const Cluster *const pCluster : *pClusterList) 
    {
        // do consider hits from 'merging clusters', as hits could be after merging point
        //if (pCluster == innerCluster || pCluster == clusterAssociation.GetAssociatedCluster())
        //continue;
        
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second) 
            {
                const CartesianVector hitPosition(pCaloHit->GetPositionVector());

                if (hitPosition.GetX() < minX || hitPosition.GetX() > maxX)
                    continue;

                
                if (hitPosition.GetZ() < minZ || hitPosition.GetZ() > maxZ)
                    continue;
                
                //float distanceFromLine(std::fabs((gradient*(hitPosition.GetX()) - hitPosition.GetZ() + zIntercept)/(std::sqrt(1 + std::pow(gradient, 2)))));
                float distanceFromLine(connectingLineDirection.GetCrossProduct(hitPosition - innerPoint).GetMagnitude());
     
                if (distanceFromLine > m_distanceFromLine)
                {
                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POTENTIAL HIT", RED, 2);
                    continue;
                }

                //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POTENTIAL HIT", GREEN, 2);
                
                extrapolatedCaloHitVector.push_back(pCaloHit);
                caloHitToParentClusterMap[pCaloHit] = pCluster;
            }
        }
    }
    
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}


//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::InitialiseSlidingFitResultMap(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap)
{
    macroSlidingFitResultMap.clear();
    microSlidingFitResultMap.clear();
    
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

    /*
    for (const Cluster *const pCluster : clusterVector)
    {
        const auto microFitIter(microSlidingFitResultMap.find(pCluster));
        if (microFitIter == microSlidingFitResultMap.end())
            continue;

        const auto macroFitIter(macroSlidingFitResultMap.find(pCluster));
        if (macroFitIter == macroSlidingFitResultMap.end())
            continue;

        ClusterList theCluster;
        theCluster.push_back(pCluster);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "Cluster", BLACK);

        const CartesianVector& microClusterDirection(microFitIter->second.GetGlobalMinLayerDirection());
        const CartesianVector& microClusterPosition(microFitIter->second.GetGlobalMinLayerPosition());

        const CartesianVector& macroClusterDirection(macroFitIter->second.GetGlobalMinLayerDirection());
        const CartesianVector& macroClusterPosition(macroFitIter->second.GetGlobalMinLayerPosition());

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &microClusterPosition, "MICRO POSITION", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &macroClusterPosition, "MACRO POSITION", BLACK, 2);
        
        const CartesianVector microPoint1(microClusterPosition + microClusterDirection*40.0);
        const CartesianVector microPoint2(microClusterPosition - microClusterDirection*40.0);

        const CartesianVector macroPoint1(macroClusterPosition + macroClusterDirection*40.0);
        const CartesianVector macroPoint2(macroClusterPosition - macroClusterDirection*40.0);

        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &microPoint1, &microPoint2, "MICRO DIRECTION", VIOLET, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &macroPoint1, &macroPoint2, "MACRO DIRECTION", BLUE, 2, 2);
        
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
        
    }
    */
   
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::FindBestClusterAssociation(const Cluster *const pCurrentCluster, const ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap, const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation)
{
    const ClusterVector::const_iterator currentIter(std::find(clusterVector.begin(), clusterVector.end(), pCurrentCluster));

    const TwoDSlidingFitResultMap::const_iterator currentMicroFitIter(microSlidingFitResultMap.find(pCurrentCluster));
    if (currentMicroFitIter == microSlidingFitResultMap.end())
        return false;

    // this should not happen
    const TwoDSlidingFitResultMap::const_iterator currentMacroFitIter(macroSlidingFitResultMap.find(pCurrentCluster));
    if (currentMacroFitIter == macroSlidingFitResultMap.end())
    {
        std::cout << "ISOBEL: BEST CLUSTER ASSOCIATION, CLUSTER IN MICRO AND NOT IN MACRO" << std::endl;
        return false;
    }

    bool foundAssociation(false);
    unsigned int maxHits(0);    
    for (ClusterVector::const_iterator testIter = currentIter; testIter != clusterVector.end(); ++testIter)
    {
        const Cluster *const pTestCluster = *testIter;
        
        if (pTestCluster == pCurrentCluster)
            continue;

        if (pTestCluster->GetNCaloHits() < maxHits)
            continue;

        const TwoDSlidingFitResultMap::const_iterator testMicroFitIter(microSlidingFitResultMap.find(pTestCluster));
        if (testMicroFitIter == microSlidingFitResultMap.end())
            continue;

        const TwoDSlidingFitResultMap::const_iterator testMacroFitIter(macroSlidingFitResultMap.find(pTestCluster));
        if (testMacroFitIter == macroSlidingFitResultMap.end())
        {
            std::cout << "ISOBEL: BEST CLUSTER ASSOCIATION, CLUSTER IN MICRO AND NOT IN MACRO" << std::endl;
            continue;
        }
        
        CartesianVector currentMergePoint(0.f, 0.f, 0.f), testMergePoint(0.f, 0.f, 0.f), currentMergeDirection(0.f, 0.f, 0.f), testMergeDirection(0.f, 0.f, 0.f);
        if (!this->GetClusterMergingCoordinates(currentMicroFitIter->second, currentMacroFitIter->second, currentMergePoint, currentMergeDirection, testMicroFitIter->second, testMacroFitIter->second, testMergePoint, testMergeDirection))
            continue;

        if (!AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection))
            continue;

        foundAssociation = true;
        maxHits = pTestCluster->GetNCaloHits();
        
        CartesianVector connectingLineDirection(testMergePoint.GetX() - currentMergePoint.GetX(), 0.f, testMergePoint.GetZ() - currentMergePoint.GetZ());
        clusterAssociation = ClusterAssociation(pTestCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection, connectingLineDirection.GetUnitVector());
    }
    
    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::AreClustersAssociated(const CartesianVector &currentPoint, const CartesianVector &currentDirection, const CartesianVector &testPoint, const CartesianVector &testDirection)
{
    if (m_limitZ)
    {
        if (testPoint.GetZ() < currentPoint.GetZ())
        {
            /*
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "Z", RED, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "Z", RED, 2);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            */
            return false;
        }
    }
    
    // check that clusters are reasonably far away
    if (std::sqrt(currentPoint.GetDistanceSquared(testPoint)) < m_minSeparationDistance)
    {
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "TOO CLOSE", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TOO CLOSE", RED, 2);
        std::cout << "DISTANCE: " << std::sqrt(currentPoint.GetDistanceSquared(testPoint)) << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        return false;
     
    }
    
    // check that opening angle is not too large
    // ATTN - HAVE TURNED THE DIRECTION AROUND IN PREVIOUS FUNCTION
    // this may be useless (if going to have it drop out at an earlier stage)
    if (currentDirection.GetCosOpeningAngle(testDirection * (-1.0)) < 0.99)
    {
        /*
        std::string reason("ANGLE" + std::to_string(currentDirection.GetCosOpeningAngle(testDirection * (-1.0))));
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, reason, RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, reason, RED, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        return false;
    }
    
    // check that fit allows you to get from one merge point to other merge point
    float separationDistance(std::sqrt(currentPoint.GetDistanceSquared(testPoint)));
    CartesianVector extrapolatedCurrentPoint(currentPoint + currentDirection*separationDistance);
    CartesianVector extrapolatedTestPoint(testPoint + testDirection*separationDistance);
    
    /*
    CartesianVector testTR(testPoint.GetX() + m_maxXSeparation, 0, testPoint.GetZ() + m_maxZSeparation);
    CartesianVector testTL(testPoint.GetX() - m_maxXSeparation, 0, testPoint.GetZ() + m_maxZSeparation);
    CartesianVector testBR(testPoint.GetX() + m_maxXSeparation, 0, testPoint.GetZ() - m_maxZSeparation);
    CartesianVector testBL(testPoint.GetX() - m_maxXSeparation, 0, testPoint.GetZ() - m_maxZSeparation);

    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTR, &testTL, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTR, &testBR, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTL, &testBL, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testBR, &testBL, "TEST BOX", RED, 2, 2);
    */
    
    if (extrapolatedCurrentPoint.GetX() > testPoint.GetX() + m_maxXSeparation || extrapolatedCurrentPoint.GetX() < testPoint.GetX() - m_maxXSeparation)
    {
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        return false;
    }

    if (extrapolatedCurrentPoint.GetZ() > testPoint.GetZ() + m_maxZSeparation || extrapolatedCurrentPoint.GetZ() < testPoint.GetZ() - m_maxZSeparation)
    {
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        return false;
    }

    /*
        CartesianVector currentTR(currentPoint.GetX() + m_maxXSeparation, 0, currentPoint.GetZ() + m_maxZSeparation);
        CartesianVector currentTL(currentPoint.GetX() - m_maxXSeparation, 0, currentPoint.GetZ() + m_maxZSeparation);
        CartesianVector currentBR(currentPoint.GetX() + m_maxXSeparation, 0, currentPoint.GetZ() - m_maxZSeparation);
        CartesianVector currentBL(currentPoint.GetX() - m_maxXSeparation, 0, currentPoint.GetZ() - m_maxZSeparation);
        
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTR, &currentTL, "CURRENT BOX", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTR, &currentBR, "CURRENT BOX", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTL, &currentBL, "CURRENT BOX", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentBR, &currentBL, "CURRENT BOX", RED, 2, 2);
    */
    
    
    if (extrapolatedTestPoint.GetX() > currentPoint.GetX() + m_maxXSeparation || extrapolatedTestPoint.GetX() < currentPoint.GetX() - m_maxXSeparation)
    {
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        return false;
    }


    if (extrapolatedTestPoint.GetZ() > currentPoint.GetZ() + m_maxZSeparation || extrapolatedTestPoint.GetZ() < currentPoint.GetZ() - m_maxZSeparation)
    {
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        return false;
    }

    /*
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTR, &testTL, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTR, &testBR, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTL, &testBL, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testBR, &testBL, "TEST BOX", RED, 2, 2);

    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTR, &currentTL, "CURRENT BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTR, &currentBR, "CURRENT BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTL, &currentBL, "CURRENT BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentBR, &currentBL, "CURRENT BOX", RED, 2, 2);
    
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", GREEN, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", VIOLET, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", GREEN, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", VIOLET, 2);
    */
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "PASS", GREEN, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "PASS", GREEN, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult, CartesianVector &currentMergePosition, CartesianVector &currentMergeDirection, const TwoDSlidingFitResult &testMicroFitResult, const TwoDSlidingFitResult &testMacroFitResult, CartesianVector &testMergePosition, CartesianVector &testMergeDirection)
{
    //////////
    //CartesianVector origin(0.f,0.f,0.f);
    //////////
    
    //CURRENT CLUSTER
    const LayerFitResultMap& currentMicroLayerFitResultMap(currentMicroFitResult.GetLayerFitResultMap());
    const LayerFitResultMap& currentMacroLayerFitResultMap(currentMacroFitResult.GetLayerFitResultMap());
    const LayerFitResultMap& testMicroLayerFitResultMap(testMicroFitResult.GetLayerFitResultMap());
    const LayerFitResultMap& testMacroLayerFitResultMap(testMacroFitResult.GetLayerFitResultMap());

    CartesianVector currentAverageDirection(0.f, 0.f, 0.f);
    CartesianVector testAverageDirection(0.f, 0.f, 0.f);

    currentMacroFitResult.GetGlobalDirection(currentMacroLayerFitResultMap.begin()->second.GetGradient(), currentAverageDirection);
    testMacroFitResult.GetGlobalDirection(testMacroLayerFitResultMap.begin()->second.GetGradient(), testAverageDirection);    

    // potentially, could be different 
    int maxLayer(std::min(currentMicroFitResult.GetMaxLayer(), currentMacroFitResult.GetMaxLayer()));
    int minLayer(std::max(currentMicroFitResult.GetMinLayer(), currentMacroFitResult.GetMinLayer()));

    unsigned int goodPositionCount(0);
    unsigned int currentStabilityHitWindow(std::ceil(currentMicroFitResult.GetCluster()->GetNCaloHits() * 0.1));
    for (int i = maxLayer; i >= minLayer; --i)
    {
        const auto microIter(currentMicroLayerFitResultMap.find(i));

        // sometimes there isn't a fit result for the layer
        if (microIter == currentMicroLayerFitResultMap.end())
            continue;
        
        CartesianVector microDirection(0.f, 0.f, 0.f);
        currentMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);

        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(m_useOtherCluster ? testAverageDirection : currentAverageDirection));
       
        //////////
        //std::cout << "COS OPENING ANGLE: " << cosDirectionOpeningAngle << std::endl;
        //currentMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), currentMergePosition);
        //////////
        
        if (cosDirectionOpeningAngle > 0.9995)
        {
            if (goodPositionCount == 0)
            {
                currentMergeDirection = currentAverageDirection;
                currentMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), currentMergePosition);
            }
            ++goodPositionCount;
        }
        else
        {
            goodPositionCount = 0;
        }

        if (goodPositionCount > currentStabilityHitWindow)
        {
            break;
        }
            
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentMergePosition, "POSITION", VIOLET, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &origin, "ORIGIN", BLUE, 2);
        //PandoraMonitoringApi::Pause(this->GetPandora());
        
        // if cannot find a point that is near the average 
        if (i == minLayer)
        {
            if (m_abortIfNoPosition)
            {
                return false;
            }
            
            currentMergePosition = currentMicroFitResult.GetGlobalMaxLayerPosition();
            currentMergeDirection = currentMicroFitResult.GetGlobalMaxLayerDirection();
            
        }
    }

    //TEST CLUSTER
    maxLayer = std::min(testMicroFitResult.GetMaxLayer(), testMacroFitResult.GetMaxLayer());
    minLayer = std::max(testMicroFitResult.GetMinLayer(), testMacroFitResult.GetMinLayer());

    goodPositionCount = 0;
    unsigned int testStabilityHitWindow(std::ceil(testMicroFitResult.GetCluster()->GetNCaloHits() * 0.1));
    //std::cout << "TEST STABILITY WINDOW: " << testStabilityHitWindow << std::endl;
    for (int i = minLayer; i <= maxLayer; ++i)
    {
        const auto microIter(testMicroLayerFitResultMap.find(i));

        if (microIter == testMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        testMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
        
        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(m_useOtherCluster ? currentAverageDirection : testAverageDirection));
        
        ///////////////
        //std::cout << "COS OPENING ANGLE: " << microDirection.GetCosOpeningAngle(m_useOtherCluster ? currentAverageDirection : testAverageDirection) << std::endl;
        //testMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), testMergePosition);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testMergePosition, "POSITION", VIOLET, 2);
        ///////////////
        
        if (cosDirectionOpeningAngle > 0.9995)
        {
            if (goodPositionCount == 0)
            {
                // so that direction vectors face one another
                testMergeDirection = testAverageDirection * (-1.0);
                testMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), testMergePosition);
            }
            ++goodPositionCount;
        }
        else
        {
            goodPositionCount = 0;
        }

        //std::cout << "POSITION COUNT: " << goodPositionCount  << std::endl;
        
        if (goodPositionCount > testStabilityHitWindow)
        {
            break;
        }

        
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &origin, "ORIGIN", BLUE, 2);
        //PandoraMonitoringApi::Pause(this->GetPandora());

        // if cannot find a point that is near the average 
        if (i == maxLayer)
        {
            if (m_abortIfNoPosition)
            {
                return false;
            }
            
	        //std::cout << "ISOBEL: COULDN'T FIND TEST AVERAGE POINT" << std::endl;
            testMergePosition = testMicroFitResult.GetGlobalMinLayerPosition();
            testMergeDirection = testMicroFitResult.GetGlobalMinLayerDirection() * (-1.0);
        }
    }
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentMergePosition, "CURRENT MERGE POSITION", BLUE, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testMergePosition, "TEST MERGE POSITION", BLUE, 2);
    PandoraMonitoringApi::Pause(this->GetPandora());
    */
    return true;

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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LimitZ", m_limitZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseOtherCluster", m_useOtherCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AbortIfNoPosition", m_abortIfNoPosition));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content


    
