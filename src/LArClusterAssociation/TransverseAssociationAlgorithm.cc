/**
 *  @file   LArContent/src/LArClusterAssociation/TransverseAssociationAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster association algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/TransverseAssociationAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

void TransverseAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // clear vector
    clusterVector.clear();

    // loop over input cluster list
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        clusterVector.push_back(*iter);

    // sort clusters
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &inputClusters, ClusterAssociationMap &clusterAssociationMap) const
{  
    // Cluster Vectors
    ClusterVector transverseClusters, longitudinalClusters, seedClusters, nonSeedClusters;

     // Separate transverse and longitudinal clusters
    // (transverse clusters are either short or not longitudinal)
    this->GetTransverseClusters(inputClusters, transverseClusters, longitudinalClusters);

    // Separate transverse clusters into seed and non-seed types
    // (seed clusters are not already associated with a longitudinal cluster)
    this->GetSeedClusters(transverseClusters, longitudinalClusters, seedClusters, nonSeedClusters);

    // Use seed clusters to create transverse cluster objects
    // (these are protoclusters with a direction and inner/outer vertices).
    LArTransverseClusterMap transverseClusterMap;
    this->FillTransverseClusterMap(seedClusters, transverseClusters, transverseClusterMap);

// ---- BEGIN DISPLAY ----
// ClusterList tempList1, tempList2, tempList3, tempList4;
// tempList1.insert(transverseClusters.begin(), transverseClusters.end());
// tempList2.insert(longitudinalClusters.begin(), longitudinalClusters.end());
// tempList3.insert(seedClusters.begin(), seedClusters.end());
// tempList4.insert(nonSeedClusters.begin(), nonSeedClusters.end());
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f); 
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "TransverseClusters", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "LongitudinalClusters", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempList3, "SeedClusters", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempList4, "NonSeedClusters", YELLOW);

// for (LArTransverseClusterMap::const_iterator iterMarker = transverseClusterMap.begin(), iterEndMarker = transverseClusterMap.end(); iterMarker != iterEndMarker; ++iterMarker)
// {
// LArTransverseCluster transverseCluster = iterMarker->second;

// float minX = transverseCluster.GetInnerVertex().GetX();
// float maxX = transverseCluster.GetOuterVertex().GetX();
// float minZ = transverseCluster.GetInnerVertex().GetZ();
// float maxZ = transverseCluster.GetOuterVertex().GetZ();
// unsigned numMarkers = 20;

// for (unsigned int n=0; n<numMarkers; ++n){
// float posX = minX + (((float)n+0.5)/((float)numMarkers))*(maxX-minX);
// float posZ = minZ + (((float)n+0.5)/((float)numMarkers))*(maxZ-minZ);
// CartesianVector marker(posX,0.f,posZ);
// PANDORA_MONITORING_API(AddMarkerToVisualization(&marker, "hit", RED, 1.));
// }  
// }

//PandoraMonitoringApi::ViewEvent();
// ---- END DISPLAY ----

    // Use transverse cluster objects to define forward/backward associations
    this->FillClusterAssociationMap(transverseClusterMap, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
    float currentMinX(0.f), currentMaxX(0.f), testMinX(0.f), testMaxX(0.f);

    this->GetExtremalCoordinatesX(pCurrentCluster, currentMinX, currentMaxX);
    this->GetExtremalCoordinatesX(pTestCluster, testMinX, testMaxX);

    if (isForward)
        return (testMaxX > currentMaxX);

    return (testMinX < currentMinX); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetExtremalCoordinatesX(const Cluster *const pCluster, float &minX, float &maxX) const
{
    minX = std::min(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()).GetX(),
        pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()).GetX());
    maxX = std::max(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()).GetX(),
        pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()).GetX());

    const OrderedCaloHitList& orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iterI = orderedCaloHitList.begin(), iterEndI = orderedCaloHitList.end(); iterI != iterEndI; ++iterI)
    {
        for (CaloHitList::const_iterator iterJ = iterI->second->begin(), iterEndJ = iterI->second->end(); iterJ != iterEndJ; ++iterJ)
        {
            const CaloHit* pCaloHit = *iterJ;

            if (pCaloHit->GetPositionVector().GetX() < minX)
                minX = pCaloHit->GetPositionVector().GetX();

            if (pCaloHit->GetPositionVector().GetX() > maxX)
                maxX = pCaloHit->GetPositionVector().GetX();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetTransverseClusters(const ClusterVector &inputVector, ClusterVector &transverseVector, ClusterVector &longitudinalVector) const
{
    // clear vectors  
    transverseVector.clear();
    longitudinalVector.clear();

    // loop over input cluster list
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;  

        // separate short and long clusters
        if (pCluster->GetNCaloHits() <= m_clusterCaloHits)
        {
            transverseVector.push_back(pCluster);
        }

        // separate transverse and longitudinal clusters
        else
        {
            const CartesianVector longitudinalDirection(0.f, 0.f, 1.f);
            const ClusterHelper::ClusterFitResult &clusterFitResult(pCluster->GetFitToAllHitsResult());

            if (clusterFitResult.IsFitSuccessful())
            {
                const CartesianVector fitDirection(clusterFitResult.GetDirection());

                // calculate dot product with respect to forward direction
                if (std::fabs(fitDirection.GetDotProduct(longitudinalDirection)) < m_clusterCosAngle)
                {
                    transverseVector.push_back(pCluster);
                }
                else
                {
                    longitudinalVector.push_back(pCluster);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetSeedClusters(const ClusterVector &shortVector, const ClusterVector &longVector, ClusterVector &seedVector,
    ClusterVector &nonSeedVector) const
{
    seedVector.clear();
    nonSeedVector.clear();

    // loop over short cluster list
    for (ClusterVector::const_iterator iterShort = shortVector.begin(), iterEndShort = shortVector.end(); iterShort != iterEndShort; ++iterShort)
    {
        Cluster *pClusterShort = *iterShort;
        bool isSeedCluster(true);

        if (pClusterShort->GetNCaloHits() <= m_clusterCaloHits)
        {
            CartesianVector position((pClusterShort->GetCentroid(pClusterShort->GetInnerPseudoLayer())
                + pClusterShort->GetCentroid(pClusterShort->GetOuterPseudoLayer())) * 0.5f);

            // loop over long cluster list
            for (ClusterVector::const_iterator iterLong = longVector.begin(), iterEndLong = longVector.end(); iterLong != iterEndLong; ++iterLong)
            {
                Cluster *pClusterLong = *iterLong;

                if (LArClusterHelper::GetClosestDistance(position, pClusterLong) < m_clusterWindow)
                {
                    isSeedCluster = false;
                    break;
                }
            }
        }

        if (isSeedCluster)
        {
            seedVector.push_back(pClusterShort);
        }
        else
        {
            nonSeedVector.push_back(pClusterShort);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetTransverseAssociatedClusters(const Cluster *const inputCluster, const ClusterVector &inputClusterVector,
    ClusterVector &outputClusterVector) const
{
    outputClusterVector.clear();

    // input cluster
    Cluster *pCluster = (Cluster*)inputCluster;  // John won't approve of this...

    //long clusters 
    if (pCluster->GetNCaloHits() > m_clusterCaloHits)
    {
        outputClusterVector.push_back(pCluster);
        return;
    }

    // short clusters
    for (ClusterVector::const_iterator iter = inputClusterVector.begin(), iterEnd = inputClusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pNearbyCluster = *iter;

        if ((pCluster->GetNCaloHits() > m_clusterCaloHits) || (pNearbyCluster->GetNCaloHits() > m_clusterCaloHits))
            continue;

        if (pCluster == pNearbyCluster)
            continue;

        if (this->IsTransverseAssociated(pCluster, pNearbyCluster))
            outputClusterVector.push_back(pNearbyCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillTransverseClusterMap(const ClusterVector &seedClusters, const ClusterVector &transverseClusters,
    LArTransverseClusterMap &transverseClusterMap) const
{
    ClusterVector associatedClusters;

    for (ClusterVector::const_iterator iter = seedClusters.begin(), iterEnd = seedClusters.end(); iter != iterEnd; ++iter)
    {
        const Cluster *pCluster = *iter;
        associatedClusters.clear();

        this->GetTransverseAssociatedClusters(pCluster, transverseClusters, associatedClusters);

        if (associatedClusters.empty())
            continue;

        (void) transverseClusterMap.insert(LArTransverseClusterMap::value_type(pCluster, LArTransverseCluster(pCluster,associatedClusters)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillClusterAssociationMap(const LArTransverseClusterMap &transverseClusterMap, ClusterAssociationMap &clusterAssociationMap) const
{
    // Compile list of all possible forward and backward associations
    LArClusterMergeMap forwardMergeMap, backwardMergeMap;

    for (LArTransverseClusterMap::const_iterator iter = transverseClusterMap.begin(), iterEnd = transverseClusterMap.end(); iter != iterEnd; ++iter)
    {
        const LArTransverseCluster& transCluster = iter->second;
        Cluster *pCentralCluster = (Cluster*)(transCluster.GetCluster());

        for (unsigned int n = 0; n < transCluster.GetNumClusters(); ++n)
        {
            Cluster *pCluster = (Cluster*)(transCluster.GetCluster(n));

            if (pCentralCluster == pCluster)
                continue;

            if (this->IsExtremalCluster(true, pCentralCluster, pCluster) && !this->IsExtremalCluster(false, pCentralCluster, pCluster))
            {
                if (forwardMergeMap[pCentralCluster].find(pCluster) == forwardMergeMap[pCentralCluster].end())
                    forwardMergeMap[pCentralCluster].insert(pCluster);

                if (backwardMergeMap[pCluster].find(pCentralCluster) == backwardMergeMap[pCluster].end())
                    backwardMergeMap[pCluster].insert(pCentralCluster);
            }

            if (this->IsExtremalCluster(false, pCentralCluster, pCluster) && !this->IsExtremalCluster(true, pCentralCluster, pCluster))
            {
                if (forwardMergeMap[pCluster].find(pCentralCluster) == forwardMergeMap[pCluster].end())
                    forwardMergeMap[pCluster].insert(pCentralCluster);

                if (backwardMergeMap[pCentralCluster].find(pCluster) == backwardMergeMap[pCentralCluster].end())
                    backwardMergeMap[pCentralCluster].insert(pCluster);
            }
        }
    }

    for (LArTransverseClusterMap::const_iterator iter1 = transverseClusterMap.begin(), iterEnd1 = transverseClusterMap.end(); iter1 != iterEnd1; ++iter1)
    {
        const LArTransverseCluster &innerCluster = iter1->second;
        Cluster *pInnerCluster = (Cluster*)(innerCluster.GetCluster());

        for (LArTransverseClusterMap::const_iterator iter2 = transverseClusterMap.begin(), iterEnd2 = transverseClusterMap.end(); iter2 != iterEnd2; ++iter2)
        {
            const LArTransverseCluster& outerCluster = iter2->second;
            Cluster* pOuterCluster = (Cluster*)(outerCluster.GetCluster());

            if (pInnerCluster == pOuterCluster)
                continue;

            if (this->IsExtremalCluster(true, pInnerCluster, pOuterCluster) && this->IsExtremalCluster(false, pOuterCluster, pInnerCluster))
            {
                if (this->IsTransverseAssociated(innerCluster, outerCluster))
                {
                    if (forwardMergeMap[pInnerCluster].find(pOuterCluster) == forwardMergeMap[pInnerCluster].end())
                        forwardMergeMap[pInnerCluster].insert(pOuterCluster);

                    if (backwardMergeMap[pOuterCluster].find(pInnerCluster) == forwardMergeMap[pOuterCluster].end())
                        backwardMergeMap[pOuterCluster].insert(pInnerCluster);
                }
            }
        }
    }

    // Select neighbouring forward associations
    for (LArClusterMergeMap::const_iterator iter1 = forwardMergeMap.begin(), iterEnd1 = forwardMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pInnerCluster = iter1->first;

        for (ClusterList::iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster *pOuterCluster = *iter2;

            if (pOuterCluster == pInnerCluster)
                throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            bool isNeighbouringCluster(true);

            for (ClusterList::iterator iter3 = iter1->second.begin(), iterEnd3 = iter1->second.end(); iter3 != iterEnd3; ++iter3)
            {
                Cluster *pMiddleCluster = *iter3;

                if (pMiddleCluster == pOuterCluster)
                    continue;

                if (forwardMergeMap[pMiddleCluster].find(pOuterCluster) != forwardMergeMap[pMiddleCluster].end())
                {
                    isNeighbouringCluster = false;
                    break;
                }
            }

            if (isNeighbouringCluster)
            {
                if (clusterAssociationMap[pInnerCluster].m_forwardAssociations.find(pOuterCluster) == clusterAssociationMap[pInnerCluster].m_forwardAssociations.end())
                    clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);

                if (clusterAssociationMap[pOuterCluster].m_backwardAssociations.find(pInnerCluster) == clusterAssociationMap[pOuterCluster].m_backwardAssociations.end())
                    clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
            }
        }
    }

    // Select neighbouring backward associations
    for (LArClusterMergeMap::const_iterator iter1 = backwardMergeMap.begin(), iterEnd1 = backwardMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pOuterCluster = iter1->first;

        for (ClusterList::iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster *pInnerCluster = *iter2;

            if (pInnerCluster == pOuterCluster)
                throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            bool isNeighbouringCluster(true);

            for (ClusterList::iterator iter3 = iter1->second.begin(), iterEnd3 = iter1->second.end(); iter3 != iterEnd3; ++iter3)
            {
                Cluster *pMiddleCluster = *iter3;

                if (pMiddleCluster == pInnerCluster)
                    continue;

                if (backwardMergeMap[pMiddleCluster].find(pInnerCluster) != backwardMergeMap[pMiddleCluster].end())
                {
                    isNeighbouringCluster = false;
                    break;
                }
            }

            if (isNeighbouringCluster)
            {
                if (clusterAssociationMap[pInnerCluster].m_forwardAssociations.find(pOuterCluster) == clusterAssociationMap[pInnerCluster].m_forwardAssociations.end())
                    clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);

                if (clusterAssociationMap[pOuterCluster].m_backwardAssociations.find(pInnerCluster) == clusterAssociationMap[pOuterCluster].m_backwardAssociations.end())
                     clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
            }
        }
    }

    // Prune single associations
    bool carryOn(true);

    while (carryOn)
    {
        carryOn = false;

        for (ClusterAssociationMap::iterator iter1 = clusterAssociationMap.begin(), iterEnd1 = clusterAssociationMap.end(); iter1 != iterEnd1; ++iter1)
        {
            Cluster *pCluster = iter1->first;

            if (pCluster->GetNCaloHits() > m_clusterCaloHits)
                continue;

            if (!iter1->second.m_backwardAssociations.empty() && !iter1->second.m_forwardAssociations.empty())
                continue;

            // check forward associations
            if (iter1->second.m_forwardAssociations.size() == 1)
            {
                Cluster *pForwardCluster = *(iter1->second.m_forwardAssociations.begin());

                if (pForwardCluster->GetNCaloHits() < m_clusterCaloHits)
                {
                    ClusterAssociationMap::iterator iter2 = clusterAssociationMap.find(pForwardCluster);

                    if (iter2 == clusterAssociationMap.end())
                        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

                    if (iter2->second.m_forwardAssociations.empty())
                    {
                        iter1->second.m_forwardAssociations.erase(pForwardCluster);

                        if (iter2->second.m_backwardAssociations.find(pCluster) != iter2->second.m_backwardAssociations.end())
                            iter2->second.m_backwardAssociations.erase(pCluster);

                        carryOn = true;
                    }
                }
            }

            // check backward associations
            if (iter1->second.m_backwardAssociations.size() == 1)
            {
                Cluster *pBackwardCluster = *(iter1->second.m_backwardAssociations.begin());

                if (pBackwardCluster->GetNCaloHits() < m_clusterCaloHits)
                {
                    ClusterAssociationMap::iterator iter2 = clusterAssociationMap.find(pBackwardCluster);

                    if (iter2 == clusterAssociationMap.end())
                        throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

                    if (iter2->second.m_backwardAssociations.empty())
                    {
                        iter1->second.m_backwardAssociations.erase(pBackwardCluster);

                        if (iter2->second.m_forwardAssociations.find(pCluster) != iter2->second.m_forwardAssociations.end())
                            iter2->second.m_forwardAssociations.erase(pCluster);

                        carryOn = true;
                    }
                }
            }
        }
    }

// ---- BEGIN DISPLAY ----
// for (ClusterAssociationMap::iterator iter = clusterAssociationMap.begin(), iterEnd = clusterAssociationMap.end(); iter != iterEnd; ++iter)
// {
// Cluster* pCentralCluster  = iter->first;
// ClusterList forwardList  = iter->second.m_forwardAssociations;
// ClusterList backwardList = iter->second.m_backwardAssociations;  
// ClusterList tempList; tempList.insert(pCentralCluster);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f); 
// PandoraMonitoringApi::VisualizeClusters(&tempList, "CentralClusters", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&forwardList, "Forward", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&backwardList, "Backward", RED);
// PandoraMonitoringApi::ViewEvent();
// }
// ---- END DISPLAY ----
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const LArTransverseCluster &transCluster1, const LArTransverseCluster &transCluster2) const
{
    // relative angle
    const CartesianVector &direction1(transCluster1.GetDirection());
    const CartesianVector &direction2(transCluster2.GetDirection());

    if (direction1.GetDotProduct(direction2) < m_minCosRelativeAngle)
        return false;

    // relative position (InnerX transCluster1 -> OuterX transCluster2)
    const Cluster *pCluster1(transCluster1.GetCluster());
    const Cluster *pCluster2(transCluster2.GetCluster());
 
    if ((transCluster1.GetDirection().GetDotProduct(transCluster2.GetInnerVertex()-transCluster1.GetInnerVertex()) > 0) &&
        (transCluster1.GetDirection().GetDotProduct(transCluster2.GetOuterVertex()-transCluster1.GetOuterVertex()) > 0) &&
        (transCluster2.GetDirection().GetDotProduct(transCluster1.GetInnerVertex()-transCluster2.GetInnerVertex()) < 0) &&
        (transCluster2.GetDirection().GetDotProduct(transCluster1.GetOuterVertex()-transCluster2.GetOuterVertex()) < 0))
    {
        return (this->IsTransverseAssociated(transCluster1, transCluster2.GetInnerVertex()) &&
            this->IsTransverseAssociated(transCluster2, transCluster1.GetOuterVertex()));
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const LArTransverseCluster &transCluster, const pandora::CartesianVector& testVertex) const
{
    const CartesianVector &innerVertex(transCluster.GetInnerVertex());
    const CartesianVector &outerVertex(transCluster.GetOuterVertex());
    const CartesianVector &direction(transCluster.GetDirection());

    if (std::fabs(direction.GetCrossProduct(testVertex - innerVertex).GetMagnitudeSquared()) > m_maxTransverseSeparation * m_maxTransverseSeparation)
    {
        return false;
    }

    if ((direction.GetDotProduct(testVertex - innerVertex) < -m_maxLongitudinalSeparation) ||
        (direction.GetDotProduct(testVertex - outerVertex) > +m_maxLongitudinalSeparation))
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    // compare X and Z centroids
    const float aveX1(0.5f * (pCluster1->GetCentroid(pCluster1->GetInnerPseudoLayer()).GetX() + pCluster1->GetCentroid(pCluster1->GetOuterPseudoLayer()).GetX()));
    const float aveZ1(0.5f * (pCluster1->GetCentroid(pCluster1->GetInnerPseudoLayer()).GetZ() + pCluster1->GetCentroid(pCluster1->GetOuterPseudoLayer()).GetZ()));
    const float aveX2(0.5f * (pCluster2->GetCentroid(pCluster2->GetInnerPseudoLayer()).GetX() + pCluster2->GetCentroid(pCluster2->GetOuterPseudoLayer()).GetX()));
    const float aveZ2(0.5f * (pCluster2->GetCentroid(pCluster2->GetInnerPseudoLayer()).GetZ() + pCluster2->GetCentroid(pCluster2->GetOuterPseudoLayer()).GetZ()));

    if ((std::fabs(aveX2 - aveX1) < m_clusterWindow) && (std::fabs(aveZ2 - aveZ1) < m_clusterWindow) &&
        (std::fabs(aveZ2 - aveZ1) < std::fabs(aveX2 - aveX1) * std::fabs(m_clusterTanAngle)))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TransverseAssociationAlgorithm::LArTransverseCluster::LArTransverseCluster(const Cluster *pSeedCluster, const ClusterVector &associatedClusters) :
    m_innerVertex(0.f,0.f,0.f),
    m_outerVertex(0.f,0.f,0.f),
    m_direction(0.f,0.f,1.f),
    m_rms(std::numeric_limits<float>::max())
{
    float Swzz(0.f);
    float Swxx(0.f);
    float Swzx(0.f);
    float Swz(0.f);
    float Swx(0.f);
    float Sw(0.f);

    float minX(std::numeric_limits<float>::max());
    float maxX(-std::numeric_limits<float>::max());

    float aveX(0.f); 
    float aveZ(0.f);

    Cluster *pCluster = (Cluster*)pSeedCluster;  // John won't approve of this...
    m_clusterVector.push_back(pCluster);

    for (ClusterVector::const_iterator iter = associatedClusters.begin(), iterEnd = associatedClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (pCluster == pSeedCluster)
            continue;

        m_clusterVector.push_back(*iter);
    }

    for (ClusterVector::const_iterator iterI = m_clusterVector.begin(), iterEndI = m_clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        const Cluster *pInputCluster = *iterI;
        const OrderedCaloHitList &orderedCaloHitList(pInputCluster->GetOrderedCaloHitList());

        for (OrderedCaloHitList::const_iterator iterJ = orderedCaloHitList.begin(), iterEndJ = orderedCaloHitList.end(); iterJ != iterEndJ; ++iterJ)
        {
            const CaloHitList *pCaloHitList = iterJ->second;

            for (CaloHitList::const_iterator iterK = pCaloHitList->begin(), iterEndK = pCaloHitList->end(); iterK != iterEndK; ++iterK)
            {
                const CaloHit *pCaloHit = *iterK;

                if (pCaloHit->GetPositionVector().GetX() < minX)
                    minX = pCaloHit->GetPositionVector().GetX();

                if (pCaloHit->GetPositionVector().GetX() > maxX)
                    maxX = pCaloHit->GetPositionVector().GetX();

                Swzz += pCaloHit->GetPositionVector().GetZ() * pCaloHit->GetPositionVector().GetZ();
                Swxx += pCaloHit->GetPositionVector().GetX() * pCaloHit->GetPositionVector().GetX();
                Swzx += pCaloHit->GetPositionVector().GetZ() * pCaloHit->GetPositionVector().GetX();
                Swz  += pCaloHit->GetPositionVector().GetZ();
                Swx  += pCaloHit->GetPositionVector().GetX();
                Sw   += 1.f;
            }
        }
    }

    if (Sw > 0.f)
    {
        aveX = Swx / Sw; 
        aveZ = Swz / Sw;

        if (Sw * Swxx - Swx * Swx > 0.f)
        {
            float m((Sw * Swzx - Swx * Swz) / (Sw * Swxx - Swx * Swx));
            float c((Swz - m * Swx) / Sw);

            float px(1.f / std::sqrt(1.f + m * m));
            float pz(m / std::sqrt(1.f + m * m));

            m_innerVertex.SetValues(minX, 0.f, aveZ + m * (minX - aveX));
            m_outerVertex.SetValues(maxX, 0.f, aveZ + m * (maxX - aveX));
            m_direction.SetValues(px, 0.f, pz);

            m_rms = std::sqrt(((Swzz + m * m * Swxx + c * c * Sw) - 2.f * (m * Swzx + c * Swz - m * c * Swx)) / Sw);
        }
        else
        {
            m_innerVertex.SetValues(aveX, 0.f, aveZ);
            m_outerVertex.SetValues(aveX, 0.f, aveZ);
            m_direction.SetValues(1.f, 0.f, 0.f);
            //m_rms = std::numeric_limits<float>::max();
        }
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
   m_clusterWindow = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterWindow", m_clusterWindow));

    m_clusterCaloHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterCaloHits", m_clusterCaloHits));

    m_clusterAngle = 45.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterAngle", m_clusterAngle));

    m_clusterCosAngle = std::cos(m_clusterAngle * M_PI / 180.f);
    m_clusterTanAngle = std::tan(m_clusterAngle * M_PI / 180.f);

    m_minCosRelativeAngle = 0.866f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    m_maxTransverseSeparation = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseSeparation", m_maxTransverseSeparation));

    m_maxLongitudinalSeparation = 7.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalSeparation", m_maxLongitudinalSeparation));

    m_minTransverseLength = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTransverseLength", m_minTransverseLength));

    m_minTransverseLayers = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTransverseLayers", m_minTransverseLayers));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
