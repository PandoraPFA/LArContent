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
    m_minCaloHits(30),
    m_maxXSeparation(10),
    m_maxZSeparation(10)
{
}

EMTrackAlgorithm::ClusterAssociation::ClusterAssociation(const Cluster *const pAssociatedCluster, const CartesianVector &innerMergePoint, const CartesianVector &innerMergeDirection, const CartesianVector &outerMergePoint, const CartesianVector &outerMergeDirection) :
    m_pAssociatedCluster(pAssociatedCluster),
    m_innerMergePoint(innerMergePoint),
    m_innerMergeDirection(innerMergeDirection),
    m_outerMergePoint(outerMergePoint),
    m_outerMergeDirection(outerMergeDirection)
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
    SelectCleanClusters(pClusterList, clusterVector);

    TwoDSlidingFitResultMap slidingFitResultMap;
    InitialiseSlidingFitResultMap(clusterVector, slidingFitResultMap);
    
    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    for (const Cluster *const innerCluster : clusterVector)
    {
        // ISOBEL - NEED TO IMPLEMENT AN EMPTY CONSTRUCTOR
        ClusterAssociation clusterAssociation(nullptr, CartesianVector(0,0,0), CartesianVector(0,0,0), CartesianVector(0,0,0), CartesianVector(0,0,0));

        // if no association is found, skip to the next cluster
        if(!this->FindBestClusterAssociation(innerCluster, clusterVector, slidingFitResultMap, clusterAssociation))
            continue;

        ClusterList currentCluster;
        currentCluster.push_back(innerCluster);

        ClusterList associatedCluster;
        associatedCluster.push_back(clusterAssociation.GetAssociatedCluster());

        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "CURRENT", GREEN);
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &associatedCluster, "ASSOCIATED", BLUE);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
        // Get a list of associated calo hits.
        CaloHitToParentClusterMap caloHitToParentClusterMap;
        CaloHitVector extrapolatedCaloHitVector;
        this->GetExtrapolatedCaloHits(innerCluster, clusterAssociation, pCaloHitList, pClusterList, extrapolatedCaloHitVector, caloHitToParentClusterMap);
        
        if (!this->IsTrackContinuous(clusterAssociation, extrapolatedCaloHitVector))
            continue;
    
        this->AddHitsToCluster(innerCluster, clusterAssociation.GetAssociatedCluster(), clusterVector, slidingFitResultMap, caloHitToParentClusterMap, extrapolatedCaloHitVector);


        // Need to add the last clusters together, delete the inner cluster fit result from map and delete the associated cluster from the cluster vector
        // then need to redo fits and then loop again.
        // maybe combine functions below in a snazzy way?
       

        /////////////////////////////
        //MONITORING PURPOSES
        ClusterList theCluster;
        theCluster.push_back(innerCluster);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "AFTER", GREEN);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        /////////////////////////////

        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), pClusterList, "ALL", BLUE);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
        // NEED TO UPDATE THE SLIDING FITS
    }

    
    
    return STATUS_CODE_SUCCESS;
}

    
//------------------------------------------------------------------------------------------------------------------------------------------
    /*
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
    */



    
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

    this->RemoveClusterFromSlidingFitResultMap(pClusterToEnlarge, slidingFitResultMap);
    this->RemoveClusterFromSlidingFitResultMap(pClusterToDelete, slidingFitResultMap);
    
    for (const CaloHit *const pCaloHit : extrapolatedCaloHitVector)
    {
        CaloHitToParentClusterMap::const_iterator caloHitParentIter(caloHitToParentClusterMap.find(pCaloHit));

        if (caloHitParentIter == caloHitToParentClusterMap.end())
        {
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
    // this is very important
    std::sort(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), LArClusterHelper::SortHitsByPosition);

    const float m_maxTrackGaps(2);
    const float m_lineSegmentLength(7);
    
    const CartesianVector &innerMergePoint(clusterAssociation.GetInnerMergePoint()), &outerMergePoint(clusterAssociation.GetOuterMergePoint());
    const CartesianVector trackDirection(CartesianVector(outerMergePoint.GetX() - innerMergePoint.GetX(), 0.f, outerMergePoint.GetZ() - innerMergePoint.GetZ()).GetUnitVector());

    std::cout << trackDirection << std::endl;
    const float zStep(std::fabs(m_lineSegmentLength ));
    unsigned int hitsInSegment(0), segmentsWithoutHits(0);

    

    CaloHitVector::const_iterator caloHitIter(extrapolatedCaloHitVector.begin());
    for (unsigned int i = 0; i < floor(std::sqrt(innerMergePoint.GetDistanceSquared(outerMergePoint)) / m_lineSegmentLength); ++i)
    {
        float zLimit(innerMergePoint.GetZ() + zStep * (static_cast<float>(i)));
        float xStep = std::sqrt(std::pow(m_lineSegmentLength, 2) - std::pow(zStep, 2));
        CartesianVector pointOnLine(innerMergePoint.GetX() - (xStep * (static_cast<float>(i))), 0, zLimit);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &pointOnLine, "BOUNDARY", BLACK, 2);
        
        while ((*caloHitIter)->GetPositionVector().GetZ() < zLimit)
        {
            ++hitsInSegment;
            ++caloHitIter;
        }

        if (!hitsInSegment)
            ++segmentsWithoutHits;

        std::cout << "HITS IN SEGMENT: " << hitsInSegment << std::endl;
        std::cout << "SEGMENTS WITHOUT HITS: " << segmentsWithoutHits << std::endl;

        PandoraMonitoringApi::Pause(this->GetPandora());

        if (segmentsWithoutHits > m_maxTrackGaps)
            return false;

        hitsInSegment = 0;
    }

    //THEN DO THE LAST ONE BUT BACKWARDS
    
    return true;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::GetExtrapolatedCaloHits(const Cluster *const innerCluster, const ClusterAssociation &clusterAssociation, const CaloHitList *const pCaloHitList, const ClusterList *const pClusterList, CaloHitVector &extrapolatedCaloHitVector, CaloHitToParentClusterMap &caloHitToParentClusterMap)
{

    CartesianVector innerPoint(clusterAssociation.GetInnerMergePoint());
    CartesianVector outerPoint(clusterAssociation.GetOuterMergePoint());
    
    float gradient(0), zIntercept(0);
    ConnectByLine(innerPoint, outerPoint, gradient, zIntercept);

    /////////////////////////////
    //MONITORING PURPOSES
    
    CartesianVector point1(innerPoint.GetX(), 0, gradient*(innerPoint.GetX()) + zIntercept);
    CartesianVector point2(outerPoint.GetX(), 0, gradient*(outerPoint.GetX()) + zIntercept);

    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &point1, &point2, "CONNECTING LINE", BLUE, 2, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerPoint, "INNER", BLUE, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerPoint, "OUTER", BLUE, 2);

    ClusterList innerClusterList;
    innerClusterList.push_back(innerCluster);
    ClusterList outerClusterList;
    outerClusterList.push_back(clusterAssociation.GetAssociatedCluster());

    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &innerClusterList, "INNER CLUSTER", BLACK);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &outerClusterList, "OUTER CLUSTER", BLACK);

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
        for (OrderedCaloHitList::const_iterator iterJ = orderedCaloHitList.begin(); iterJ != orderedCaloHitList.end(); ++iterJ){
            for(CaloHitList::const_iterator hitIter = iterJ->second->begin(); hitIter != iterJ->second->end(); ++hitIter)
            {
                const CaloHit *const pCaloHit = *hitIter;
                const CartesianVector hitPosition(pCaloHit->GetPositionVector());


                if (hitPosition.GetX() < minX || hitPosition.GetX() > maxX)
                    continue;

                if (hitPosition.GetZ() < minZ || hitPosition.GetZ() > maxZ)
                    continue;
                
                float distanceFromLine(std::fabs((gradient*(hitPosition.GetX()) - hitPosition.GetZ() + zIntercept)/(std::sqrt(1 + std::pow(gradient, 2)))));
                
                if (distanceFromLine > m_distanceFromLine)
                {
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POTENTIAL HIT", RED, 2);
                    continue;
                }


                /////////////////////////////
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POTENTIAL HIT", GREEN, 2);
                //std::cout << "DISTANCE FROM LINE" << distanceFromLine << std::endl;
                //PandoraMonitoringApi::Pause(this->GetPandora());
                /////////////////////////////
                
                extrapolatedCaloHitVector.push_back(pCaloHit);
                caloHitToParentClusterMap[pCaloHit] = pCluster;

            }
        }
    }

    
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
                
        float distanceFromLine(std::fabs((gradient*(hitPosition.GetX()) - hitPosition.GetZ() + zIntercept)/(std::sqrt(1 + std::pow(gradient, 2)))));

        if (distanceFromLine > m_distanceFromLine)
        {
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POTENTIAL HIT", RED, 2);
            continue;
        }

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POTENTIAL HIT", GREEN, 2);
        extrapolatedCaloHitVector.push_back(pCaloHit);                  
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::ConnectByLine(const CartesianVector &innerCoordinate, const CartesianVector &outerCoordinate, float &gradient, float &zIntercept)
{
    gradient = (outerCoordinate.GetZ() - innerCoordinate.GetZ())/(outerCoordinate.GetX() - innerCoordinate.GetX());
    zIntercept = outerCoordinate.GetZ() - gradient*outerCoordinate.GetX();
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

    unsigned int maxHits(0);
    bool foundAssociation(false);
    
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

        CartesianVector currentPoint(0,0,0), testPoint(0,0,0);
        CartesianVector currentDirection(0,0,0), testDirection(0,0,0);
        GetClusterMergingCoordinates(currentFitIter->second, currentPoint, currentDirection, testFitIter->second, testPoint, testDirection);

        if (!AreClustersAssociated(currentPoint, currentDirection, testPoint, testDirection))
            continue;

        foundAssociation = true;
        maxHits = pTestCluster->GetNCaloHits();
        clusterAssociation = ClusterAssociation(pTestCluster, currentPoint, currentDirection, testPoint, testDirection);
    }

    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::AreClustersAssociated(const CartesianVector &currentPoint, const CartesianVector &currentDirection, const CartesianVector &testPoint, const CartesianVector &testDirection)
{

    // check that clusters are reasonably far away
    if (std::sqrt(currentPoint.GetDistanceSquared(testPoint) < 50))
    {
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "DISTANCE: " + std::to_string(std::sqrt(currentPoint.GetDistanceSquared(testPoint))), RED);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    
    // check that opening angle is not too large
    if (currentDirection.GetCosOpeningAngle(testDirection*(-1.0)) < 0.97)
    {
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "ANGLE: " + std::to_string(std::fabs(currentDirection.GetCosOpeningAngle(testDirection))), RED);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    
    // check that fit allows you to get from one merge point to other merge point
    float separationDistance(std::sqrt(currentPoint.GetDistanceSquared(testPoint)));
    CartesianVector extrapolatedCurrentPoint(currentPoint + currentDirection*separationDistance);
    CartesianVector extrapolatedTestPoint(testPoint + testDirection*separationDistance);

    //CartesianVector testTR(testPoint.GetX() + m_maxXSeparation, 0, testPoint.GetZ() + m_maxZSeparation);
    //CartesianVector testTL(testPoint.GetX() - m_maxXSeparation, 0, testPoint.GetZ() + m_maxZSeparation);
    //CartesianVector testBR(testPoint.GetX() + m_maxXSeparation, 0, testPoint.GetZ() - m_maxZSeparation);
    //CartesianVector testBL(testPoint.GetX() - m_maxXSeparation, 0, testPoint.GetZ() - m_maxZSeparation);

    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTR, &testTL, "TEST BOX", RED, 2, 2);
    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTR, &testBR, "TEST BOX", RED, 2, 2);
    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testTL, &testBL, "TEST BOX", RED, 2, 2);
    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &testBR, &testBL, "TEST BOX", RED, 2, 2);

    
    if (extrapolatedCurrentPoint.GetX() > testPoint.GetX() + m_maxXSeparation || extrapolatedCurrentPoint.GetX() < testPoint.GetX() - m_maxXSeparation)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "EXTRAPOLATED CURRENT X TOO FAR", RED);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

        if (extrapolatedCurrentPoint.GetZ() > testPoint.GetZ() + m_maxZSeparation || extrapolatedCurrentPoint.GetZ() < testPoint.GetZ() - m_maxZSeparation)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "EXTRAPOLATED CURRENT Z TOO FAR", RED);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    //CartesianVector currentTR(currentPoint.GetX() + m_maxXSeparation, 0, currentPoint.GetZ() + m_maxZSeparation);
    //CartesianVector currentTL(currentPoint.GetX() - m_maxXSeparation, 0, currentPoint.GetZ() + m_maxZSeparation);
    //CartesianVector currentBR(currentPoint.GetX() + m_maxXSeparation, 0, currentPoint.GetZ() - m_maxZSeparation);
    //CartesianVector currentBL(currentPoint.GetX() - m_maxXSeparation, 0, currentPoint.GetZ() - m_maxZSeparation);

    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTR, &currentTL, "CURRENT BOX", RED, 2, 2);
    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTR, &currentBR, "CURRENT BOX", RED, 2, 2);
    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentTL, &currentBL, "CURRENT BOX", RED, 2, 2);
    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &currentBR, &currentBL, "CURRENT BOX", RED, 2, 2);


    
    if (extrapolatedTestPoint.GetX() > currentPoint.GetX() + m_maxXSeparation || extrapolatedTestPoint.GetX() < currentPoint.GetX() - m_maxXSeparation)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "EXTRAPOLATED TEST X TOO FAR", RED);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }


    if (extrapolatedTestPoint.GetZ() > currentPoint.GetZ() + m_maxZSeparation || extrapolatedTestPoint.GetZ() < currentPoint.GetZ() - m_maxZSeparation)
    {
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "EXTRAPOLATED TEST Z TOO FAR", RED);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "PASSED!", GREEN);
    
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    void EMTrackAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &cluster1FitResult, CartesianVector &cluster1MergePoint, CartesianVector &cluster1Direction, const TwoDSlidingFitResult &cluster2FitResult, CartesianVector &cluster2MergePoint, CartesianVector &cluster2Direction)
{

    CartesianVector minLayerPosition1(cluster1FitResult.GetGlobalMinLayerPosition()), maxLayerPosition1(cluster1FitResult.GetGlobalMaxLayerPosition());
    CartesianVector minLayerPosition2(cluster2FitResult.GetGlobalMinLayerPosition()), maxLayerPosition2(cluster2FitResult.GetGlobalMaxLayerPosition());

    CartesianVector minLayerDirection1(cluster1FitResult.GetGlobalMinLayerDirection()), maxLayerDirection1(cluster1FitResult.GetGlobalMaxLayerDirection());
    CartesianVector minLayerDirection2(cluster2FitResult.GetGlobalMinLayerDirection()), maxLayerDirection2(cluster2FitResult.GetGlobalMaxLayerDirection());

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
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZSeparation", m_maxZSeparation));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content