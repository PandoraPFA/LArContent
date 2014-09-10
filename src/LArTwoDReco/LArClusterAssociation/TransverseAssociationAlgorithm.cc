/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the transverse association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void TransverseAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    clusterVector.clear();
    clusterVector.insert(clusterVector.begin(), pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &allClusters, ClusterAssociationMap &clusterAssociationMap) const
{
    TransverseClusterList transverseClusterList;

    try
    {
        // Step 1: Sort the input clusters into sub-samples
        //         (a) shortClusters:  below first length cut
        //         (b) mediumClusters: between first and second length cuts (separated into transverse and longitudinal)
        //         (c) longClusters:   above second length cut
        ClusterVector shortClusters, transverseMediumClusters, longitudinalMediumClusters, longClusters;
        this->SortInputClusters(allClusters, shortClusters, transverseMediumClusters, longitudinalMediumClusters, longClusters);

        ClusterVector transverseClusters;
        transverseClusters.insert(transverseClusters.end(), shortClusters.begin(), shortClusters.end());
        transverseClusters.insert(transverseClusters.end(), transverseMediumClusters.begin(), transverseMediumClusters.end());

        ClusterVector establishedClusters;
        establishedClusters.insert(establishedClusters.end(), transverseMediumClusters.begin(), transverseMediumClusters.end());
        establishedClusters.insert(establishedClusters.end(), longitudinalMediumClusters.begin(), longitudinalMediumClusters.end());
        establishedClusters.insert(establishedClusters.end(), longClusters.begin(), longClusters.end());

        // Step 2: Form loose transverse associations between short clusters,
        //         without hopping over any established clusters
        ClusterAssociationMap firstAssociationMap;
        this->FillReducedAssociationMap(shortClusters, establishedClusters, firstAssociationMap);

        // Step 3: Form transverse cluster objects. Basically, try to assign a direction to each
        //         of the clusters in the 'transverseClusters' list. For the short clusters in
        //         this list, the direction is obtained from a straight line fit to its associated
        //         clusters as selected in the previous step.
        this->FillTransverseClusterList(transverseClusters, firstAssociationMap, transverseClusterList);

// ---- BEGIN DISPLAY ----
// ClusterList tempList1, tempList2, tempList3;
// tempList1.insert(transverseClusters.begin(), transverseClusters.end());
// for (unsigned int nCluster = 0; nCluster<transverseClusterList.size(); ++nCluster)
// {
// LArTransverseCluster* transverseCluster = transverseClusterList.at(nCluster);
// Cluster* pCluster = transverseCluster->GetSeedCluster();
// tempList2.insert(pCluster);
// if( tempList1.count(pCluster)) tempList1.erase(pCluster);
// }
// tempList3.insert(longitudinalMediumClusters.begin(), longitudinalMediumClusters.end());
// tempList3.insert(longClusters.begin(), longClusters.end());
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "ShortClusters", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "TransverseClusters", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempList3, "LongClusters", BLUE);
// PandoraMonitoringApi::ViewEvent();
// ---- END DISPLAY ----

    // Step 4: Form loose transverse associations between transverse clusters
    //         (First, associate medium clusters, without hopping over long clusters
    //          Next, associate all transverse clusters, without hopping over any clusters)
    ClusterAssociationMap secondAssociationMap;
    this->FillReducedAssociationMap(transverseMediumClusters, longClusters, secondAssociationMap);
    this->FillReducedAssociationMap(transverseClusters, allClusters, secondAssociationMap);

    // Step 5: Form associations between transverse cluster objects
    //         (These transverse associations must already exist as loose associations
    //          between transverse clusters as identified in the previous step).
    ClusterAssociationMap transverseAssociationMap;
    this->FillTransverseAssociationMap(transverseClusterList, secondAssociationMap, transverseAssociationMap);

// ---- BEGIN DISPLAY ----
// ClusterList tempList, tempListForward, tempListBackward;
// for (ClusterAssociationMap::iterator iterI = transverseAssociationMap.begin(), iterEndI = transverseAssociationMap.end(); iterI != iterEndI; ++iterI)
// {
// Cluster *pCluster = iterI->first;
// for (ClusterList::iterator iterJ = iterI->second.m_forwardAssociations.begin(), iterEndJ = iterI->second.m_forwardAssociations.end(); iterJ != iterEndJ; ++iterJ)
// tempListForward.insert(*iterJ);
// for (ClusterList::iterator iterJ = iterI->second.m_backwardAssociations.begin(), iterEndJ = iterI->second.m_backwardAssociations.end(); iterJ != iterEndJ; ++iterJ)
// tempListBackward.insert(*iterJ);
// }
// for (ClusterAssociationMap::iterator iterI = transverseAssociationMap.begin(), iterEndI = transverseAssociationMap.end(); iterI != iterEndI; ++iterI)
// {
// Cluster *pCluster = iterI->first;
// if(tempListForward.count(pCluster) > 0 && tempListBackward.count(pCluster) > 0)
// {
// tempList.insert(pCluster);
// tempListForward.erase(pCluster);
// tempListBackward.erase(pCluster);
// }
// }
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "Associations", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempListForward, "ForwardOnly", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempListBackward, "BackwardOnly", GREEN);
// PandoraMonitoringApi::ViewEvent();
// ---- END DISPLAY ----

        // Step 6: Finalise the forward/backward transverse associations by symmetrising the
        //         transverse association map and removing any double-counting
        this->FinalizeClusterAssociationMap(transverseAssociationMap, clusterAssociationMap);
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "TransverseAssociationAlgorithm: exception " << statusCodeException.ToString() << std::endl;
    }

    for (TransverseClusterList::const_iterator iter = transverseClusterList.begin(), iterEnd = transverseClusterList.end(); iter != iterEnd; ++iter)
    {
        delete *iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::SortInputClusters(const ClusterVector &inputVector, ClusterVector &shortVector, ClusterVector &transverseMediumVector, ClusterVector &longitudinalMediumVector, ClusterVector &longVector) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        const float clusterLengthT(this->GetTransverseSpan(pCluster));
        const float clusterLengthL(this->GetLongitudinalSpan(pCluster));
        const float clusterLengthSquared(clusterLengthT * clusterLengthT + clusterLengthL * clusterLengthL);

        if (clusterLengthSquared < m_firstLengthCut * m_firstLengthCut)
        {
            shortVector.push_back(pCluster);
        }
        else if (clusterLengthSquared < m_secondLengthCut * m_secondLengthCut)
        {
            if (clusterLengthL < clusterLengthT * std::fabs(m_clusterTanAngle))
                transverseMediumVector.push_back(pCluster);
            else
                longitudinalMediumVector.push_back(pCluster);
        }
        else
        {
            longVector.push_back(pCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillAssociationMap(const ClusterVector &firstVector, const ClusterVector &secondVector, ClusterAssociationMap &firstAssociationMap, ClusterAssociationMap &secondAssociationMap) const
{
    for (ClusterVector::const_iterator iterI = firstVector.begin(), iterEndI = firstVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster *pClusterI = *iterI;

        for (ClusterVector::const_iterator iterJ = secondVector.begin(), iterEndJ = secondVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster *pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
            continue;

            if (this->IsAssociated(true, pClusterI, pClusterJ))
            {
                firstAssociationMap[pClusterI].m_forwardAssociations.insert(pClusterJ);
                secondAssociationMap[pClusterJ].m_backwardAssociations.insert(pClusterI);
            }

            if (this->IsAssociated(false, pClusterI, pClusterJ))
            {
                firstAssociationMap[pClusterI].m_backwardAssociations.insert(pClusterJ);
                secondAssociationMap[pClusterJ].m_forwardAssociations.insert(pClusterI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillReducedAssociationMap(const ClusterVector &firstVector, const ClusterVector &secondVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // To build a 'reduced' association map, form associations between clusters in the first cluster vector,
    // but prevent these associations from hopping over any clusters in the second cluster vector.
    // i.e. A->B from the first vector is forbidden if A->C->B exists with C from the second vector

    ClusterAssociationMap firstAssociationMap, firstAssociationMapSwapped;
    ClusterAssociationMap secondAssociationMap, secondAssociationMapSwapped;

    this->FillAssociationMap(firstVector, firstVector, firstAssociationMap, firstAssociationMapSwapped);
    this->FillAssociationMap(firstVector, secondVector, secondAssociationMap, secondAssociationMapSwapped);
    this->FillReducedAssociationMap(firstAssociationMap, secondAssociationMap, secondAssociationMapSwapped, clusterAssociationMap);

// ---- BEGIN DISPLAY ----
// for (ClusterAssociationMap::iterator iter = clusterAssociationMap.begin(), iterEnd = clusterAssociationMap.end(); iter != iterEnd; ++iter)
// {
// Cluster *pCluster = iter->first;
// ClusterList tempList;
// tempList.insert(pCluster);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "CentralCluster", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&iter->second.m_forwardAssociations,"ForwardAssociations",RED);
// PandoraMonitoringApi::VisualizeClusters(&iter->second.m_backwardAssociations,"BackwardAssociations",GREEN);
// PandoraMonitoringApi::ViewEvent();
// }
// ---- END DISPLAY ----
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillTransverseClusterList(const ClusterVector &inputVector, const ClusterAssociationMap &inputAssociationMap, TransverseClusterList &transverseClusterList) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;
        ClusterVector associatedClusters;

        this->GetAssociatedClusters(pCluster, inputAssociationMap, associatedClusters);

        if (this->GetTransverseSpan(pCluster, associatedClusters) < m_transverseClusterMinLength)
            continue;

// ---- BEGIN DISPLAY ----
// ClusterList tempList1, tempList2;
// tempList1.insert(pCluster);
// for (unsigned int nCluster = 0; nCluster<associatedClusters.size(); ++nCluster)
// {
// Cluster* pTempCluster = associatedClusters.at(nCluster);
// tempList2.insert(pTempCluster);
// }
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "CentralCluster", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "AssociatedClusters", GREEN);
// PandoraMonitoringApi::ViewEvent();
// ---- END DISPLAY ----
        transverseClusterList.push_back(new LArTransverseCluster(pCluster, associatedClusters));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillTransverseAssociationMap(const TransverseClusterList &transverseClusterList, const ClusterAssociationMap& transverseAssociationMap, ClusterAssociationMap &clusterAssociationMap) const
{
    for (TransverseClusterList::const_iterator iter1 = transverseClusterList.begin(), iterEnd1 = transverseClusterList.end(); iter1 != iterEnd1; ++iter1)
    {
        LArTransverseCluster *pInnerTransverseCluster = *iter1;
        Cluster *pInnerCluster(pInnerTransverseCluster->GetSeedCluster());

        ClusterAssociationMap::const_iterator iterInner = transverseAssociationMap.find(pInnerCluster);
        if (transverseAssociationMap.end() == iterInner)
            continue;

        for (TransverseClusterList::const_iterator iter2 = transverseClusterList.begin(), iterEnd2 = transverseClusterList.end(); iter2 != iterEnd2; ++iter2)
        {
            LArTransverseCluster *pOuterTransverseCluster = *iter2;
            Cluster *pOuterCluster(pOuterTransverseCluster->GetSeedCluster());

            ClusterAssociationMap::const_iterator iterOuter = transverseAssociationMap.find(pOuterCluster);
            if(transverseAssociationMap.end() == iterOuter)
                continue;

            if (pInnerCluster == pOuterCluster)
                continue;

            if (iterInner->second.m_forwardAssociations.count(pOuterCluster) == 0 ||
                iterOuter->second.m_backwardAssociations.count(pInnerCluster) == 0)
                continue;

            if (!this->IsExtremalCluster(true, pInnerCluster, pOuterCluster) ||
                !this->IsExtremalCluster(false, pOuterCluster, pInnerCluster))
                continue;

            if (!this->IsTransverseAssociated(pInnerTransverseCluster, pOuterTransverseCluster))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetAssociatedClusters(Cluster *const pClusterI, const ClusterAssociationMap& associationMap, ClusterVector &associatedVector) const
{
    ClusterAssociationMap::const_iterator iterI = associationMap.find(pClusterI);
    if (associationMap.end() == iterI)
        return;

    for (ClusterList::iterator iterJ = iterI->second.m_forwardAssociations.begin(), iterEndJ = iterI->second.m_forwardAssociations.end(); iterJ != iterEndJ; ++iterJ)
    {
        Cluster *pClusterJ = *iterJ;

        if (this->IsTransverseAssociated(pClusterI, pClusterJ))
            associatedVector.push_back(pClusterJ);
    }

    for (ClusterList::iterator iterJ = iterI->second.m_backwardAssociations.begin(), iterEndJ = iterI->second.m_backwardAssociations.end(); iterJ != iterEndJ; ++iterJ)
    {
        Cluster *pClusterJ = *iterJ;

        if (this->IsTransverseAssociated(pClusterJ, pClusterI))
            associatedVector.push_back(pClusterJ);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsAssociated(const bool isForward, const Cluster *const pFirstCluster, const Cluster *const pSecondCluster) const
{
    CartesianVector firstInner(0.f,0.f,0.f), firstOuter(0.f,0.f,0.f);
    CartesianVector secondInner(0.f,0.f,0.f), secondOuter(0.f,0.f,0.f);
    this->GetExtremalCoordinatesX(pFirstCluster, firstInner, firstOuter);
    this->GetExtremalCoordinatesX(pSecondCluster, secondInner, secondOuter);

    const CartesianVector firstCoordinate(isForward ? firstOuter : firstInner);
    const CartesianVector secondCoordinate(isForward ? secondOuter : secondInner);

    if ((firstCoordinate.GetZ() > std::max(secondInner.GetZ(),secondOuter.GetZ()) + m_maxLongitudinalOverlap) ||
        (firstCoordinate.GetZ() < std::min(secondInner.GetZ(),secondOuter.GetZ()) - m_maxLongitudinalOverlap))
        return false;

    if ((isForward && secondCoordinate.GetX() < firstCoordinate.GetX()) ||
        (!isForward && secondCoordinate.GetX() > firstCoordinate.GetX()))
        return false;

    const CartesianVector firstProjection(LArClusterHelper::GetClosestPosition(firstCoordinate, pSecondCluster));

    if((isForward && firstProjection.GetX() < firstCoordinate.GetX() - m_maxTransverseOverlap) ||
        (!isForward && firstProjection.GetX() > firstCoordinate.GetX() + m_maxTransverseOverlap))
        return false;

    if((isForward && firstProjection.GetX() > firstCoordinate.GetX() + m_clusterWindow) ||
        (!isForward && firstProjection.GetX() < firstCoordinate.GetX() - m_clusterWindow))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
{
    CartesianVector innerInner(0.f,0.f,0.f), innerOuter(0.f,0.f,0.f);
    CartesianVector outerInner(0.f,0.f,0.f), outerOuter(0.f,0.f,0.f);
    this->GetExtremalCoordinatesX(pInnerCluster, innerInner, innerOuter);
    this->GetExtremalCoordinatesX(pOuterCluster, outerInner, outerOuter);

    const CartesianVector innerCentroid((innerInner + innerOuter) * 0.5);
    const CartesianVector outerCentroid((outerInner + outerOuter) * 0.5);

    if ((std::fabs(innerCentroid.GetZ() - outerInner.GetZ()) > std::fabs(innerCentroid.GetX() - outerInner.GetX()) * std::fabs(m_clusterTanAngle)) &&
        (std::fabs(innerCentroid.GetZ() - outerOuter.GetZ()) > std::fabs(innerCentroid.GetX() - outerOuter.GetX()) * std::fabs(m_clusterTanAngle)))
        return false;

    if ((std::fabs(outerCentroid.GetZ() - innerInner.GetZ()) > std::fabs(outerCentroid.GetX() - innerInner.GetX()) * std::fabs(m_clusterTanAngle)) &&
        (std::fabs(outerCentroid.GetZ() - innerOuter.GetZ()) > std::fabs(outerCentroid.GetX() - innerOuter.GetX()) * std::fabs(m_clusterTanAngle)))
        return false;

    const CartesianVector innerProjection(LArClusterHelper::GetClosestPosition(outerInner, pInnerCluster));
    const CartesianVector outerProjection(LArClusterHelper::GetClosestPosition(innerOuter, pOuterCluster));

    if (innerOuter.GetX() > innerProjection.GetX() + m_maxTransverseOverlap ||
    outerInner.GetX() < outerProjection.GetX() - m_maxTransverseOverlap)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const LArTransverseCluster *const pInnerTransverseCluster,
    const LArTransverseCluster *const pOuterTransverseCluster) const
{
    if (pInnerTransverseCluster->GetDirection().GetDotProduct(pOuterTransverseCluster->GetDirection()) < m_transverseClusterMinCosTheta)
        return false;

    if (!this->IsTransverseAssociated(pInnerTransverseCluster, pOuterTransverseCluster->GetInnerVertex()))
        return false;

    if (!this->IsTransverseAssociated(pOuterTransverseCluster, pInnerTransverseCluster->GetOuterVertex()))
        return false;

    if (!this->IsTransverseAssociated(pInnerTransverseCluster->GetSeedCluster(),pOuterTransverseCluster->GetSeedCluster()))
        return false;

    if (this->IsOverlapping(pInnerTransverseCluster->GetSeedCluster(),pOuterTransverseCluster->GetSeedCluster()))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const LArTransverseCluster *const pTransverseCluster, const CartesianVector& testVertex) const
{
    const CartesianVector &innerVertex(pTransverseCluster->GetInnerVertex());
    const CartesianVector &outerVertex(pTransverseCluster->GetOuterVertex());
    const CartesianVector &direction(pTransverseCluster->GetDirection());

    if (std::fabs(direction.GetCrossProduct(testVertex - innerVertex).GetMagnitudeSquared()) > m_transverseClusterMaxDisplacement * m_transverseClusterMaxDisplacement)
        return false;

    if ((direction.GetDotProduct(testVertex - innerVertex) < -m_clusterWindow) || (direction.GetDotProduct(testVertex - outerVertex) > +m_clusterWindow))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsOverlapping(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
{
    CartesianVector innerInner(0.f,0.f,0.f), innerOuter(0.f,0.f,0.f);
    CartesianVector outerInner(0.f,0.f,0.f), outerOuter(0.f,0.f,0.f);
    this->GetExtremalCoordinatesX(pInnerCluster, innerInner, innerOuter);
    this->GetExtremalCoordinatesX(pOuterCluster, outerInner, outerOuter);

    const CartesianVector innerProjection(LArClusterHelper::GetClosestPosition(outerInner, pInnerCluster));
    const CartesianVector outerProjection(LArClusterHelper::GetClosestPosition(innerOuter, pOuterCluster));

    const float innerOverlapSquared((innerProjection - innerOuter).GetMagnitudeSquared());
    const float outerOverlapSquared((outerProjection - outerInner).GetMagnitudeSquared());

    return (std::max(innerOverlapSquared,outerOverlapSquared) > m_maxProjectedOverlap * m_maxProjectedOverlap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TransverseAssociationAlgorithm::GetTransverseSpan(const Cluster *const pCluster) const
{
    float minX(+std::numeric_limits<float>::max());
    float maxX(-std::numeric_limits<float>::max());

    this->GetExtremalCoordinatesX(pCluster, minX, maxX);

    return (maxX - minX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TransverseAssociationAlgorithm::GetLongitudinalSpan(const Cluster *const pCluster) const
{
    float minZ(+std::numeric_limits<float>::max());
    float maxZ(-std::numeric_limits<float>::max());

    this->GetExtremalCoordinatesZ(pCluster, minZ, maxZ);

    return (maxZ - minZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TransverseAssociationAlgorithm::GetTransverseSpan(const Cluster *const pCentralCluster, const ClusterVector &associatedClusters) const
{
    float overallMinX(+std::numeric_limits<float>::max());
    float overallMaxX(-std::numeric_limits<float>::max());

    this->GetExtremalCoordinatesX(pCentralCluster, overallMinX, overallMaxX);

    float localMinX(+std::numeric_limits<float>::max());
    float localMaxX(-std::numeric_limits<float>::max());

    for (ClusterVector::const_iterator iter = associatedClusters.begin(), iterEnd = associatedClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pAssociatedCluster = *iter;

        this->GetExtremalCoordinatesX(pAssociatedCluster, localMinX, localMaxX);

        if (localMinX < overallMinX)
            overallMinX = localMinX;

        if (localMaxX > overallMaxX)
            overallMaxX = localMaxX;
    }

    return (overallMaxX - overallMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{
    float currentMinX(0.f), currentMaxX(0.f);
    this->GetExtremalCoordinatesX(pCurrentCluster, currentMinX, currentMaxX);

    float testMinX(0.f), testMaxX(0.f);
    this->GetExtremalCoordinatesX(pTestCluster, testMinX, testMaxX);

    if (isForward)
    {
        return (testMaxX > currentMaxX);
    }
    else
    {
        return (testMinX < currentMinX);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetExtremalCoordinatesX(const Cluster *const pCluster, float &minX, float &maxX) const
{
    return this->GetExtremalCoordinatesXZ(pCluster, true, minX, maxX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetExtremalCoordinatesZ(const Cluster *const pCluster, float &minZ, float &maxZ) const
{
    return this->GetExtremalCoordinatesXZ(pCluster, false, minZ, maxZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetExtremalCoordinatesXZ(const Cluster *const pCluster, const bool useX, float &minXZ, float &maxXZ) const
{
    minXZ = +std::numeric_limits<float>::max();
    maxXZ = -std::numeric_limits<float>::max();

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const float caloHitXZ(useX ? (*hitIter)->GetPositionVector().GetX() : (*hitIter)->GetPositionVector().GetZ());

            if (caloHitXZ < minXZ)
                minXZ = caloHitXZ;

            if (caloHitXZ > maxXZ)
                maxXZ = caloHitXZ;
        }
    }

    if (maxXZ < minXZ)
    throw pandora::StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetExtremalCoordinatesX(const Cluster *const pCluster, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate) const
{
    CartesianVector firstCoordinate(0.f,0.f,0.f), secondCoordinate(0.f,0.f,0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster, firstCoordinate, secondCoordinate);

    innerCoordinate = (firstCoordinate.GetX() < secondCoordinate.GetX() ? firstCoordinate : secondCoordinate);
    outerCoordinate = (firstCoordinate.GetX() > secondCoordinate.GetX() ? firstCoordinate : secondCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillReducedAssociationMap(const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const
{
    return this->FillReducedAssociationMap(inputAssociationMap, inputAssociationMap, inputAssociationMap, outputAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillReducedAssociationMap(const ClusterAssociationMap &firstAssociationMap, const ClusterAssociationMap &secondAssociationMap, const ClusterAssociationMap &secondAssociationMapSwapped, ClusterAssociationMap &clusterAssociationMap) const
{
    // Remove associations A->B from the first association map
    // if A->C exists in the second map and C->B exists in the reversed second map

    // Method can also be accessed through FillReducedAssociationMap(input,output) method,
    // which will remove association A->B from the input map if an association A->C and C->B
    // already exists in the map.

    for (ClusterAssociationMap::const_iterator iterFirst = firstAssociationMap.begin(), iterEndFirst = firstAssociationMap.end(); iterFirst != iterEndFirst; ++iterFirst)
    {
        Cluster *pCluster = iterFirst->first;

        ClusterAssociationMap::const_iterator iterSecond = secondAssociationMap.find(pCluster);

        // Remove double-counting in forward associations
        for (ClusterList::iterator iterOuter = iterFirst->second.m_forwardAssociations.begin(), iterEndOuter = iterFirst->second.m_forwardAssociations.end(); iterOuter != iterEndOuter; ++iterOuter)
        {
            Cluster *pOuterCluster = *iterOuter;

            bool isNeighbouringCluster(true);

            if (secondAssociationMap.end() != iterSecond)
            {
                for (ClusterList::iterator iterMiddle = iterSecond->second.m_forwardAssociations.begin(), iterEndMiddle = iterSecond->second.m_forwardAssociations.end(); iterMiddle != iterEndMiddle; ++iterMiddle)
                {
                     Cluster *pMiddleCluster = *iterMiddle;

                     ClusterAssociationMap::const_iterator iterSecondCheck = secondAssociationMapSwapped.find(pMiddleCluster);
                     if (secondAssociationMapSwapped.end() == iterSecondCheck)
                     continue;

                     if (iterSecondCheck->second.m_forwardAssociations.count(pOuterCluster) > 0)
                     {
                     isNeighbouringCluster = false;
                     break;
                     }
                }
            }

            if (isNeighbouringCluster)
            clusterAssociationMap[pCluster].m_forwardAssociations.insert(pOuterCluster);
        }


        // Remove double-counting in backward associations
        for (ClusterList::iterator iterInner = iterFirst->second.m_backwardAssociations.begin(), iterEndInner = iterFirst->second.m_backwardAssociations.end(); iterInner != iterEndInner; ++iterInner)
        {
            Cluster *pInnerCluster = *iterInner;

            bool isNeighbouringCluster(true);

            if (secondAssociationMap.end() != iterSecond)
            {
                for (ClusterList::iterator iterMiddle = iterSecond->second.m_backwardAssociations.begin(), iterEndMiddle = iterSecond->second.m_backwardAssociations.end(); iterMiddle != iterEndMiddle; ++iterMiddle)
                {
                     Cluster *pMiddleCluster = *iterMiddle;

                     ClusterAssociationMap::const_iterator iterSecondCheck = secondAssociationMapSwapped.find(pMiddleCluster);
                     if (secondAssociationMapSwapped.end() == iterSecondCheck)
                     continue;

                     if (iterSecondCheck->second.m_backwardAssociations.count(pInnerCluster) > 0)
                     {
                     isNeighbouringCluster = false;
                     break;
                     }
                 }
            }

            if (isNeighbouringCluster)
                clusterAssociationMap[pCluster].m_backwardAssociations.insert(pInnerCluster);
        }

    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillSymmetricAssociationMap(const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const
{
    // Generate a symmetrised association map, so that both A--Fwd-->B and B--Bwd-->A both exist.
    // If A is associated to B through both a backward and forward association (very bad!),
    // try to rationalise this through majority voting, otherwise remove the association.

    for (ClusterAssociationMap::const_iterator iter = inputAssociationMap.begin(), iterEnd = inputAssociationMap.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = iter->first;

        // Symmetrise forward associations
        for (ClusterList::iterator iterForward = iter->second.m_forwardAssociations.begin(), iterEndForward = iter->second.m_forwardAssociations.end(); iterForward != iterEndForward; ++iterForward)
        {
            Cluster *pForwardCluster = *iterForward;

            int nCounter(+1);

            if (iter->second.m_backwardAssociations.count(pForwardCluster))
                --nCounter;

            ClusterAssociationMap::const_iterator iterCheck = inputAssociationMap.find(pForwardCluster);
            if (inputAssociationMap.end() != iterCheck)
            {
                if (iterCheck->second.m_forwardAssociations.count(pCluster))
                    --nCounter;

                if (iterCheck->second.m_backwardAssociations.count(pCluster))
                    ++nCounter;
            }

            if (nCounter > 0)
            {
                if(!(outputAssociationMap[pCluster].m_backwardAssociations.count(pForwardCluster) == 0 &&
                     outputAssociationMap[pForwardCluster].m_forwardAssociations.count(pCluster) == 0))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                outputAssociationMap[pCluster].m_forwardAssociations.insert(pForwardCluster);
                outputAssociationMap[pForwardCluster].m_backwardAssociations.insert(pCluster);
            }
        }

        // Symmetrise backward associations
        for (ClusterList::iterator iterBackward = iter->second.m_backwardAssociations.begin(), iterEndBackward = iter->second.m_backwardAssociations.end(); iterBackward != iterEndBackward; ++iterBackward)
        {
            Cluster *pBackwardCluster = *iterBackward;

            int nCounter(-1);

            if (iter->second.m_forwardAssociations.count(pBackwardCluster))
            ++nCounter;

            ClusterAssociationMap::const_iterator iterCheck = inputAssociationMap.find(pBackwardCluster);
            if (inputAssociationMap.end() != iterCheck)
            {
                if (iterCheck->second.m_backwardAssociations.count(pCluster))
                    ++nCounter;

                if (iterCheck->second.m_forwardAssociations.count(pCluster))
                    --nCounter;
            }

            if (nCounter < 0)
            {
                if(!(outputAssociationMap[pCluster].m_forwardAssociations.count(pBackwardCluster) == 0 &&
                     outputAssociationMap[pBackwardCluster].m_backwardAssociations.count(pCluster) == 0))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                outputAssociationMap[pCluster].m_backwardAssociations.insert(pBackwardCluster);
                outputAssociationMap[pBackwardCluster].m_forwardAssociations.insert(pCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FinalizeClusterAssociationMap(const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const
{
    ClusterAssociationMap intermediateAssociationMap;
    this->FillSymmetricAssociationMap(inputAssociationMap, intermediateAssociationMap);
    this->FillReducedAssociationMap(intermediateAssociationMap, outputAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TransverseAssociationAlgorithm::LArTransverseCluster::LArTransverseCluster(Cluster *pSeedCluster, const ClusterVector &associatedClusters) :
    m_pSeedCluster(pSeedCluster),
    m_associatedClusters(associatedClusters),
    m_innerVertex(0.f, 0.f, 0.f),
    m_outerVertex(0.f, 0.f, 0.f),
    m_direction(0.f, 0.f, 0.f)
{
    float Swzz(0.f), Swxx(0.f), Swzx(0.f), Swz(0.f), Swx(0.f), Sw(0.f);
    float minX(std::numeric_limits<float>::max());
    float maxX(-std::numeric_limits<float>::max());

    ClusterList clusterList;
    clusterList.insert(pSeedCluster);
    clusterList.insert(associatedClusters.begin(), associatedClusters.end());

    for (ClusterList::const_iterator iterI = clusterList.begin(), iterEndI = clusterList.end(); iterI != iterEndI; ++iterI)
    {
        for (OrderedCaloHitList::const_iterator iterJ = (*iterI)->GetOrderedCaloHitList().begin(), iterEndJ = (*iterI)->GetOrderedCaloHitList().end(); iterJ != iterEndJ; ++iterJ)
        {
            for (CaloHitList::const_iterator iterK = iterJ->second->begin(), iterEndK = iterJ->second->end(); iterK != iterEndK; ++iterK)
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
        const float averageX(Swx / Sw);
        const float averageZ(Swz / Sw);

        if (Sw * Swxx - Swx * Swx > 0.f)
        {
            float m((Sw * Swzx - Swx * Swz) / (Sw * Swxx - Swx * Swx));
            float px(1.f / std::sqrt(1.f + m * m));
            float pz(m / std::sqrt(1.f + m * m));

            m_innerVertex.SetValues(minX, 0.f, averageZ + m * (minX - averageX));
            m_outerVertex.SetValues(maxX, 0.f, averageZ + m * (maxX - averageX));
            m_direction.SetValues(px, 0.f, pz);
        }
        else
        {
            m_innerVertex.SetValues(averageX, 0.f, averageZ);
            m_outerVertex.SetValues(averageX, 0.f, averageZ);
            m_direction.SetValues(1.f, 0.f, 0.f);
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
    m_firstLengthCut = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FirstLengthCut", m_firstLengthCut));

    m_secondLengthCut = 7.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SecondLengthCut", m_secondLengthCut));

    m_clusterWindow = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterWindow", m_clusterWindow));

    m_clusterAngle = 45.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "clusterAngle", m_clusterAngle));

    m_clusterCosAngle = std::cos(m_clusterAngle * M_PI / 180.f);
    m_clusterTanAngle = std::tan(m_clusterAngle * M_PI / 180.f);

    m_maxTransverseOverlap = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseOverlap", m_maxTransverseOverlap));

    m_maxProjectedOverlap = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxProjectedOverlap", m_maxProjectedOverlap));

    m_maxLongitudinalOverlap = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalOverlap", m_maxLongitudinalOverlap));

    m_transverseClusterMinCosTheta = 0.866f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TransverseClusterMinCosTheta", m_transverseClusterMinCosTheta));

    m_transverseClusterMinLength = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TransverseClusterMinLength", m_transverseClusterMinLength));

    m_transverseClusterMaxDisplacement = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TransverseClusterMaxDisplacement", m_transverseClusterMaxDisplacement));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
