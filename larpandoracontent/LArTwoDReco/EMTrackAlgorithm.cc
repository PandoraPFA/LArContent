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

using namespace pandora;

namespace lar_content
{


EMTrackAlgorithm::EMTrackAlgorithm() :
    m_caloHitListName(),
    m_caloHitToParentClusterMap(),
    m_minCaloHits(25),
    m_maxXSeparation(10),
    m_maxZSeparation(10),
    m_useUnclusteredHits(),
    m_useAxisDirection()
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

    TwoDSlidingFitResultMap slidingFitResultMap;
    this->InitialiseSlidingFitResultMap(clusterVector, slidingFitResultMap);
    
    //PandoraMonitoringApi::Create(this->GetPandora());
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    //ClusterList clusterList(clusterVector.begin(), clusterVector.end());
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &clusterList, "Selected Clusters", BLACK);
    /*
    for (const Cluster *const pCluster : clusterVector)
    {
        const CartesianVector position(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POINT", RED, 2);
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
            if(!this->FindBestClusterAssociation(innerCluster, clusterVector, slidingFitResultMap, clusterAssociation))
                continue;

            //ClusterList currentCluster;
            //currentCluster.push_back(innerCluster);

            //ClusterList associatedCluster;
            //associatedCluster.push_back(clusterAssociation.GetAssociatedCluster());

            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "CURRENT", GREEN);
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &associatedCluster, "ASSOCIATED", BLUE);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
            // Get a list of associated calo hits.
            CaloHitToParentClusterMap caloHitToParentClusterMap;
            CaloHitVector extrapolatedCaloHitVector;
            this->GetExtrapolatedCaloHits(innerCluster, clusterAssociation, pCaloHitList, pClusterList, extrapolatedCaloHitVector, caloHitToParentClusterMap);

            if (extrapolatedCaloHitVector.empty())
                continue;
        
            if (!this->IsTrackContinuous(clusterAssociation, extrapolatedCaloHitVector))
                continue;

            //std::cout << "HERE" << std::endl;
            
            this->AddHitsToCluster(innerCluster, clusterAssociation.GetAssociatedCluster(), clusterVector, slidingFitResultMap, caloHitToParentClusterMap, extrapolatedCaloHitVector);

            /////////////////////////////
            //MONITORING PURPOSES
            //ClusterList theCluster;
            //theCluster.push_back(innerCluster);
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "AFTER", GREEN);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            /////////////////////////////

            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), pClusterList, "ALL", BLUE);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());


            this->UpdateSlidingFitResultMap(clusterVector, slidingFitResultMap);

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
    
    // sort hits by Z, then X and then Y
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void EMTrackAlgorithm::UpdateSlidingFitResultMap(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const float slidingFitWindow(20);
    
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResultMap::const_iterator fitResultIter(slidingFitResultMap.find(pCluster));

        if(fitResultIter == slidingFitResultMap.end())
        {
            try {(void) slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, slidingFitWindow, slidingFitPitch)));}
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

void EMTrackAlgorithm::AddHitsToCluster(const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete, ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap, const CaloHitToParentClusterMap &caloHitToParentClusterMap, const CaloHitVector &extrapolatedCaloHitVector)
{

    // remove from sliding fit map clusters that will be deleted or whose constituents will change
    this->RemoveClusterFromSlidingFitResultMap(pClusterToEnlarge, slidingFitResultMap);
    this->RemoveClusterFromSlidingFitResultMap(pClusterToDelete, slidingFitResultMap);
    
    for (const CaloHit *const pCaloHit : extrapolatedCaloHitVector)
    {
        CaloHitToParentClusterMap::const_iterator caloHitParentIter(caloHitToParentClusterMap.find(pCaloHit));

        if (caloHitParentIter == caloHitToParentClusterMap.end())
        {
            std::cout << "CALO HIT HAS NO PARENT" << std::endl;
            PandoraContentApi::AddToCluster(*this, pClusterToEnlarge, pCaloHit);
        }
        else
        {
            const StatusCode statusCode(PandoraContentApi::RemoveFromCluster(*this, caloHitParentIter->second, pCaloHit));
            
            if (statusCode == STATUS_CODE_SUCCESS)
            {
                RemoveClusterFromSlidingFitResultMap(caloHitParentIter->second, slidingFitResultMap);
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

        //std::cout << "ITERATION" << i << std::endl;
        
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

        // merge remainder segment into preceding segment, split into two if remainder segment lenght is over half of the m_lineSegmentLength
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

        //std::cout << "LOWER BOUNDARY: " << lowerBoundary << std::endl;
        //std::cout << "UPPER BOUNDARY: " << upperBoundary << std::endl;
        
        while (this->IsInLineSegment(lowerBoundary, upperBoundary, (*caloHitIter)->GetPositionVector()))
        {
            ++hitsInSegment;
            ++caloHitIter;

            if (caloHitIter == extrapolatedCaloHitVector.end())
                break;
        }

        if (!hitsInSegment)
            ++segmentsWithoutHits;

        //std::cout << "HITS IN SEGMENT: " << hitsInSegment << std::endl;
        //std::cout << "SEGMENTS WITHOUT HITS: " << segmentsWithoutHits << std::endl;

        if (segmentsWithoutHits > m_maxTrackGaps)
            return false;

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

    
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "POINT", GREEN, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &lowerBoundary, "LOWER BOUNDARY", BLACK, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upperBoundary, "UPPER BOUNDARY", BLACK, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upper, "UPPER", VIOLET, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &lower, "LOWER", VIOLET, 2);
            
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

void EMTrackAlgorithm::GetExtrapolatedCaloHits(const Cluster *const innerCluster, const ClusterAssociation &clusterAssociation, const CaloHitList *const pCaloHitList, const ClusterList *const pClusterList, CaloHitVector &extrapolatedCaloHitVector, CaloHitToParentClusterMap &caloHitToParentClusterMap)
{

    const CartesianVector &innerPoint(clusterAssociation.GetInnerMergePoint());
    const CartesianVector &outerPoint(clusterAssociation.GetOuterMergePoint());
    const CartesianVector &connectingLineDirection(clusterAssociation.GetConnectingLineDirection());
    

    /////////////////////////////
    //MONITORING PURPOSES

    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &innerPoint, &outerPoint, "CONNECTING LINE", BLUE, 2, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerPoint, "INNER", BLUE, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerPoint, "OUTER", BLUE, 2);

    //ClusterList innerClusterList;
    //innerClusterList.push_back(innerCluster);
    //ClusterList outerClusterList;
    //outerClusterList.push_back(clusterAssociation.GetAssociatedCluster());

    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &innerClusterList, "INNER CLUSTER", BLACK);
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &outerClusterList, "OUTER CLUSTER", BLACK);

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    /////////////////////////////

    float m_distanceFromLine(0.35);
    float minX(std::min(innerPoint.GetX(), outerPoint.GetX()));
    float maxX(std::max(innerPoint.GetX(), outerPoint.GetX()));
    float minZ(std::min(innerPoint.GetZ(), outerPoint.GetZ()));
    float maxZ(std::max(innerPoint.GetZ(), outerPoint.GetZ()));

    
    for (const Cluster *const pCluster : *pClusterList) 
    {
        // do not consider hits from 'merging clusters'
        if (pCluster == innerCluster || pCluster == clusterAssociation.GetAssociatedCluster())
            continue;
        
        
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

                /////////////////////////////
                //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POTENTIAL HIT", GREEN, 2);
                //std::cout << "DISTANCE FROM LINE: " << distanceFromLine << std::endl;
                //std::cout << "NEW DISTANCE FROM LINE: " << distanceFromLine << std::endl;
                //PandoraMonitoringApi::Pause(this->GetPandora());
                /////////////////////////////
                
                extrapolatedCaloHitVector.push_back(pCaloHit);
                caloHitToParentClusterMap[pCaloHit] = pCluster;
            }
        }
    }

    if (m_useUnclusteredHits)
    {
        // now collect any unclustered hits
        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

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

            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POTENTIAL HIT", VIOLET, 2);
            extrapolatedCaloHitVector.push_back(pCaloHit);                  
        }
    }
    
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}


//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::InitialiseSlidingFitResultMap(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap)
{
    slidingFitResultMap.clear();
    
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const float slidingFitWindow(20);
    
    for (const Cluster *const pCluster : clusterVector)
    {
        try {(void) slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, slidingFitWindow, slidingFitPitch)));}
        catch (StatusCodeException &) {}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::FindBestClusterAssociation(const Cluster *const pCurrentCluster, const ClusterVector &clusterVector, const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterAssociation &clusterAssociation)
{
    const ClusterVector::const_iterator currentIter(std::find(clusterVector.begin(), clusterVector.end(), pCurrentCluster));

    const TwoDSlidingFitResultMap::const_iterator currentFitIter(slidingFitResultMap.find(pCurrentCluster));
    if (currentFitIter == slidingFitResultMap.end())
        return false;

    bool foundAssociation(false);
    unsigned int maxHits(0);

    
    for (ClusterVector::const_iterator testIter = currentIter; testIter != clusterVector.end(); ++testIter)
    {
        const Cluster *const pTestCluster = *testIter;
        
        if (pTestCluster == pCurrentCluster)
            continue;

        if (pTestCluster->GetNCaloHits() < maxHits)
            continue;

        const TwoDSlidingFitResultMap::const_iterator testFitIter(slidingFitResultMap.find(pTestCluster));
        if (testFitIter == slidingFitResultMap.end())
            continue;
        
        CartesianVector currentMergePoint(0.f, 0.f, 0.f), testMergePoint(0.f, 0.f, 0.f), currentMergeDirection(0.f, 0.f, 0.f), testMergeDirection(0.f, 0.f, 0.f);
        this->GetClusterMergingCoordinates(currentFitIter->second, currentMergePoint, currentMergeDirection, testFitIter->second, testMergePoint, testMergeDirection);

        
        if (!AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection))
            continue;

        foundAssociation = true;
        maxHits = pTestCluster->GetNCaloHits();
        
        CartesianVector connectingLineDirection(testMergePoint.GetX() - currentMergePoint.GetX(), 0.f, testMergePoint.GetZ() - currentMergePoint.GetZ());
        clusterAssociation = ClusterAssociation(pTestCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection, connectingLineDirection.GetUnitVector());
    }
    

    return foundAssociation;
    
    /*
    for (ClusterVector::const_iterator testIter = currentIter; testIter != clusterVector.end(); ++testIter)
    {
        const Cluster *const pTestCluster = *testIter;
        
        if (pTestCluster == pCurrentCluster)
            continue;

        if (pTestCluster->GetNCaloHits() < maxHits)
            continue;

        const TwoDSlidingFitResultMap::const_iterator testFitIter(slidingFitResultMap.find(pTestCluster));
        if (testFitIter == slidingFitResultMap.end())
            continue;
        
        CartesianVector currentMergePoint(0.f, 0.f, 0.f), testMergePoint(0.f, 0.f, 0.f), currentMergeDirection(0.f, 0.f, 0.f), testMergeDirection(0.f, 0.f, 0.f);
        this->GetClusterMergingCoordinates(currentFitIter->second, currentMergePoint, currentMergeDirection, testFitIter->second, testMergePoint, testMergeDirection);

        if (!AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection))
            continue;

        foundAssociation = true;
        maxHits = pTestCluster->GetNCaloHits();
        
        CartesianVector connectingLineDirection(testMergePoint.GetX() - currentMergePoint.GetX(), 0.f, testMergePoint.GetZ() - currentMergePoint.GetZ());
        clusterAssociation = ClusterAssociation(pTestCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection, connectingLineDirection.GetUnitVector());
        break;
    }
    

    return foundAssociation;
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::AreClustersAssociated(const CartesianVector &currentPoint, const CartesianVector &currentDirection, const CartesianVector &testPoint, const CartesianVector &testDirection)
{

    // check that clusters are reasonably far away
    if (std::sqrt(currentPoint.GetDistanceSquared(testPoint)) < 30)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "TOO FAR AWAY", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TOO FAR AWAY", BLUE, 2);
        //std::cout << "DISTANCE: " << std::sqrt(currentPoint.GetDistanceSquared(testPoint)) << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    
    // check that opening angle is not too large
    if (currentDirection.GetCosOpeningAngle(testDirection*(-1.0)) < 0.99)
    {
        //std::string reason("ANGLE" + std::to_string(currentDirection.GetCosOpeningAngle(testDirection*(-1.0))));
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "ANGLE", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "ANGLE", BLUE, 2);
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
    */
    
    /*
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTR, &testTL, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTR, &testBR, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTL, &testBL, "TEST BOX", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testBR, &testBL, "TEST BOX", RED, 2, 2);
    */
    
    if (extrapolatedCurrentPoint.GetX() > testPoint.GetX() + m_maxXSeparation || extrapolatedCurrentPoint.GetX() < testPoint.GetX() - m_maxXSeparation)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    if (extrapolatedCurrentPoint.GetZ() > testPoint.GetZ() + m_maxZSeparation || extrapolatedCurrentPoint.GetZ() < testPoint.GetZ() - m_maxZSeparation)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

        /*
        CartesianVector currentTR(currentPoint.GetX() + m_maxXSeparation, 0, currentPoint.GetZ() + m_maxZSeparation);
        CartesianVector currentTL(currentPoint.GetX() - m_maxXSeparation, 0, currentPoint.GetZ() + m_maxZSeparation);
        CartesianVector currentBR(currentPoint.GetX() + m_maxXSeparation, 0, currentPoint.GetZ() - m_maxZSeparation);
        CartesianVector currentBL(currentPoint.GetX() - m_maxXSeparation, 0, currentPoint.GetZ() - m_maxZSeparation);
        */
        
        /*
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTR, &currentTL, "CURRENT BOX", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTR, &currentBR, "CURRENT BOX", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTL, &currentBL, "CURRENT BOX", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentBR, &currentBL, "CURRENT BOX", RED, 2, 2);
        */

    
    if (extrapolatedTestPoint.GetX() > currentPoint.GetX() + m_maxXSeparation || extrapolatedTestPoint.GetX() < currentPoint.GetX() - m_maxXSeparation)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }


    if (extrapolatedTestPoint.GetZ() > currentPoint.GetZ() + m_maxZSeparation || extrapolatedTestPoint.GetZ() < currentPoint.GetZ() - m_maxZSeparation)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
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
    
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

    
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "PASS", GREEN, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "PASS", GREEN, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    void EMTrackAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &cluster1FitResult, CartesianVector &cluster1MergePoint, CartesianVector &cluster1Direction, const TwoDSlidingFitResult &cluster2FitResult, CartesianVector &cluster2MergePoint, CartesianVector &cluster2Direction)
{

    CartesianVector minLayerPosition1(cluster1FitResult.GetGlobalMinLayerPosition()), maxLayerPosition1(cluster1FitResult.GetGlobalMaxLayerPosition());
    CartesianVector minLayerPosition2(cluster2FitResult.GetGlobalMinLayerPosition()), maxLayerPosition2(cluster2FitResult.GetGlobalMaxLayerPosition());

    CartesianVector minLayerDirection1(0.f, 0.f, 0.f), maxLayerDirection1(0.f, 0.f, 0.f), minLayerDirection2(0.f, 0.f, 0.f), maxLayerDirection2(0.f, 0.f, 0.f);


    if (m_useAxisDirection)
    {
        minLayerDirection1 = cluster1FitResult.GetAxisDirection();
        maxLayerDirection1 = cluster1FitResult.GetAxisDirection();
        minLayerDirection2 = cluster2FitResult.GetAxisDirection();
        maxLayerDirection2 = cluster2FitResult.GetAxisDirection();
    }
    else
    {
        minLayerDirection1 = cluster1FitResult.GetGlobalMinLayerDirection();
        maxLayerDirection1 = cluster1FitResult.GetGlobalMaxLayerDirection();
        minLayerDirection2 = cluster2FitResult.GetGlobalMinLayerDirection();
        maxLayerDirection2 = cluster2FitResult.GetGlobalMaxLayerDirection();
    }
    
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &minLayerPosition1, "minLayerPosition1", VIOLET, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &maxLayerPosition1, "maxLayerPosition1", VIOLET, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &minLayerPosition2, "minLayerPosition2", VIOLET, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &maxLayerPosition2, "maxLayerPosition2", VIOLET, 2);
    
    PandoraMonitoringApi::Pause(this->GetPandora());
    */
    
    float separationDistance(std::numeric_limits<float>::max());

     float ii(minLayerPosition1.GetDistanceSquared(minLayerPosition2));
     if (ii < separationDistance)
     {
         cluster1MergePoint = minLayerPosition1;
         cluster2MergePoint = minLayerPosition2;
         cluster1Direction = minLayerDirection1;
         cluster2Direction = minLayerDirection2;
         separationDistance = ii;
     }

     float io(minLayerPosition1.GetDistanceSquared(maxLayerPosition2));
     if (io < separationDistance)
     {
         cluster1MergePoint = minLayerPosition1;
         cluster2MergePoint = maxLayerPosition2;
         cluster1Direction = minLayerDirection1;
         cluster2Direction = maxLayerDirection2;
         separationDistance = io;
     }

     float oi(maxLayerPosition1.GetDistanceSquared(minLayerPosition2));
     if (oi < separationDistance)
     {
         cluster1MergePoint = maxLayerPosition1;
         cluster2MergePoint = minLayerPosition2;
         cluster1Direction = maxLayerDirection1;
         cluster2Direction = minLayerDirection2;
         separationDistance = oi;
     }

     float oo(maxLayerPosition1.GetDistanceSquared(maxLayerPosition2));
     if (oo < separationDistance)
     {
         cluster1MergePoint = maxLayerPosition1;
         cluster2MergePoint = maxLayerPosition2;
         cluster1Direction = maxLayerDirection1;
         cluster2Direction = maxLayerDirection2;
         separationDistance = oo;
     }                 
     
     //make sure directions are pointing at each other
     CartesianVector cluster1To2(cluster2MergePoint - cluster1MergePoint);
     CartesianVector cluster2To1(cluster1MergePoint - cluster2MergePoint);
     if (cluster1To2.GetCosOpeningAngle(cluster1Direction*(-1)) > cluster1To2.GetCosOpeningAngle(cluster1Direction))
         cluster1Direction *= -1.0;

     if (cluster2To1.GetCosOpeningAngle(cluster2Direction*(-1)) > cluster2To1.GetCosOpeningAngle(cluster2Direction))
         cluster2Direction *= -1.0;
     
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode EMTrackAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "UseAxisDirection", m_useAxisDirection));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnclusteredHits", m_useUnclusteredHits));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZSeparation", m_maxZSeparation));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content

/*
void EMTrackAlgorithm::ConnectByLine(const CartesianVector &innerCoordinate, const CartesianVector &outerCoordinate, float &gradient, float &zIntercept)
{
    gradient = (outerCoordinate.GetZ() - innerCoordinate.GetZ())/(outerCoordinate.GetX() - innerCoordinate.GetX());
    zIntercept = outerCoordinate.GetZ() - gradient*outerCoordinate.GetX();
}
*/


    
