/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the transverse association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

TransverseAssociationAlgorithm::TransverseAssociationAlgorithm() :
    m_firstLengthCut(1.5f),
    m_secondLengthCut(7.5f),
    m_clusterWindow(3.f),
    m_clusterAngle(45.f),
    m_clusterCosAngle(std::cos(m_clusterAngle * M_PI / 180.f)),
    m_clusterTanAngle(std::tan(m_clusterAngle * M_PI / 180.f)),
    m_maxTransverseOverlap(0.5f),
    m_maxProjectedOverlap(1.f),
    m_maxLongitudinalOverlap(1.5f),
    m_transverseClusterMinCosTheta(0.866f),
    m_transverseClusterMinLength(0.5f),
    m_transverseClusterMaxDisplacement(1.5f),
    m_searchRegionX(3.5f),
    m_searchRegionZ(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

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
        ClusterToClustersMap nearbyClusters;
        this->GetNearbyClusterMap(allClusters, nearbyClusters);

        if (nearbyClusters.empty())
            return;

        // Step 1: Sort the input clusters into sub-samples
        //         (a) shortClusters:  below first length cut
        //         (b) mediumClusters: between first and second length cuts (separated into transverse and longitudinal)
        //         (c) longClusters:   above second length cut
        ClusterVector shortClusters, transverseMediumClusters, longitudinalMediumClusters, longClusters;
        this->SortInputClusters(allClusters, shortClusters, transverseMediumClusters, longitudinalMediumClusters, longClusters);

        ClusterVector transverseClusters(shortClusters.begin(), shortClusters.end());
        transverseClusters.insert(transverseClusters.end(), transverseMediumClusters.begin(), transverseMediumClusters.end());

        ClusterVector establishedClusters(transverseMediumClusters.begin(), transverseMediumClusters.end());
        establishedClusters.insert(establishedClusters.end(), longitudinalMediumClusters.begin(), longitudinalMediumClusters.end());
        establishedClusters.insert(establishedClusters.end(), longClusters.begin(), longClusters.end());

        // Step 2: Form loose transverse associations between short clusters,
        //         without hopping over any established clusters
        ClusterAssociationMap firstAssociationMap;
        this->FillReducedAssociationMap(nearbyClusters, shortClusters, establishedClusters, firstAssociationMap);

        // Step 3: Form transverse cluster objects. Basically, try to assign a direction to each
        //         of the clusters in the 'transverseClusters' list. For the short clusters in
        //         this list, the direction is obtained from a straight line fit to its associated
        //         clusters as selected in the previous step.
        this->FillTransverseClusterList(nearbyClusters, transverseClusters, firstAssociationMap, transverseClusterList);

        // Step 4: Form loose transverse associations between transverse clusters
        //         (First, associate medium clusters, without hopping over long clusters
        //          Next, associate all transverse clusters, without hopping over any clusters)
        ClusterAssociationMap secondAssociationMap;
        this->FillReducedAssociationMap(nearbyClusters, transverseMediumClusters, longClusters, secondAssociationMap);
        this->FillReducedAssociationMap(nearbyClusters, transverseClusters, allClusters, secondAssociationMap);

        // Step 5: Form associations between transverse cluster objects
        //         (These transverse associations must already exist as loose associations
        //          between transverse clusters as identified in the previous step).
        ClusterAssociationMap transverseAssociationMap;
        this->FillTransverseAssociationMap(nearbyClusters, transverseClusterList, secondAssociationMap, transverseAssociationMap);

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

void TransverseAssociationAlgorithm::GetNearbyClusterMap(const ClusterVector &allClusters, ClusterToClustersMap &nearbyClusters) const
{
    HitToClusterMap hitToClusterMap;
    CaloHitList allCaloHits;

    for (const Cluster *const pCluster : allClusters)
    {
        CaloHitList daughterHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(daughterHits);
        allCaloHits.insert(allCaloHits.end(), daughterHits.begin(), daughterHits.end());

        for (const CaloHit *const pCaloHit : daughterHits)
            (void)hitToClusterMap.insert(HitToClusterMap::value_type(pCaloHit, pCluster));
    }

    HitKDTree2D kdTree;
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(allCaloHits, hitKDNode2DList));
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    for (const Cluster *const pCluster : allClusters)
    {
        CaloHitList daughterHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(daughterHits);

        for (const CaloHit *const pCaloHit : daughterHits)
        {
            KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegionX, m_searchRegionZ));

            HitKDNode2DList found;
            kdTree.search(searchRegionHits, found);

            for (const auto &hit : found)
                (void)nearbyClusters[pCluster].insert(hitToClusterMap.at(hit.data));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::SortInputClusters(const ClusterVector &inputVector, ClusterVector &shortVector,
    ClusterVector &transverseMediumVector, ClusterVector &longitudinalMediumVector, ClusterVector &longVector) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

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

    std::sort(shortVector.begin(), shortVector.end(), LArClusterHelper::SortByNHits);
    std::sort(transverseMediumVector.begin(), transverseMediumVector.end(), LArClusterHelper::SortByNHits);
    std::sort(longitudinalMediumVector.begin(), longitudinalMediumVector.end(), LArClusterHelper::SortByNHits);
    std::sort(longVector.begin(), longVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillReducedAssociationMap(const ClusterToClustersMap &nearbyClusters, const ClusterVector &firstVector,
    const ClusterVector &secondVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // To build a 'reduced' association map, form associations between clusters in the first cluster vector,
    // but prevent these associations from hopping over any clusters in the second cluster vector.
    // i.e. A->B from the first vector is forbidden if A->C->B exists with C from the second vector

    ClusterAssociationMap firstAssociationMap, firstAssociationMapSwapped;
    ClusterAssociationMap secondAssociationMap, secondAssociationMapSwapped;

    this->FillAssociationMap(nearbyClusters, firstVector, firstVector, firstAssociationMap, firstAssociationMapSwapped);
    this->FillAssociationMap(nearbyClusters, firstVector, secondVector, secondAssociationMap, secondAssociationMapSwapped);
    this->FillReducedAssociationMap(firstAssociationMap, secondAssociationMap, secondAssociationMapSwapped, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillAssociationMap(const ClusterToClustersMap &nearbyClusters, const ClusterVector &firstVector,
    const ClusterVector &secondVector, ClusterAssociationMap &firstAssociationMap, ClusterAssociationMap &secondAssociationMap) const
{
    for (ClusterVector::const_iterator iterI = firstVector.begin(), iterEndI = firstVector.end(); iterI != iterEndI; ++iterI)
    {
        const Cluster *const pClusterI = *iterI;

        for (ClusterVector::const_iterator iterJ = secondVector.begin(), iterEndJ = secondVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            const Cluster *const pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
                continue;

            if (this->IsAssociated(true, pClusterI, pClusterJ, nearbyClusters))
            {
                firstAssociationMap[pClusterI].m_forwardAssociations.insert(pClusterJ);
                secondAssociationMap[pClusterJ].m_backwardAssociations.insert(pClusterI);
            }

            if (this->IsAssociated(false, pClusterI, pClusterJ, nearbyClusters))
            {
                firstAssociationMap[pClusterI].m_backwardAssociations.insert(pClusterJ);
                secondAssociationMap[pClusterJ].m_forwardAssociations.insert(pClusterI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillTransverseClusterList(const ClusterToClustersMap &nearbyClusters, const ClusterVector &inputVector,
    const ClusterAssociationMap &inputAssociationMap, TransverseClusterList &transverseClusterList) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;
        const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster))};
        const float transverseClusterMinLengthAdjusted{ratio * m_transverseClusterMinLength};

        ClusterVector associatedClusters;

        this->GetAssociatedClusters(nearbyClusters, pCluster, inputAssociationMap, associatedClusters);

        if (this->GetTransverseSpan(pCluster, associatedClusters) < transverseClusterMinLengthAdjusted)
            continue;

        transverseClusterList.push_back(new LArTransverseCluster(pCluster, associatedClusters));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillTransverseAssociationMap(const ClusterToClustersMap &nearbyClusters,
    const TransverseClusterList &transverseClusterList, const ClusterAssociationMap &transverseAssociationMap,
    ClusterAssociationMap &clusterAssociationMap) const
{
    for (TransverseClusterList::const_iterator iter1 = transverseClusterList.begin(), iterEnd1 = transverseClusterList.end(); iter1 != iterEnd1; ++iter1)
    {
        LArTransverseCluster *const pInnerTransverseCluster = *iter1;
        const Cluster *const pInnerCluster(pInnerTransverseCluster->GetSeedCluster());

        ClusterAssociationMap::const_iterator iterInner = transverseAssociationMap.find(pInnerCluster);
        if (transverseAssociationMap.end() == iterInner)
            continue;

        for (TransverseClusterList::const_iterator iter2 = transverseClusterList.begin(), iterEnd2 = transverseClusterList.end();
             iter2 != iterEnd2; ++iter2)
        {
            LArTransverseCluster *const pOuterTransverseCluster = *iter2;
            const Cluster *const pOuterCluster(pOuterTransverseCluster->GetSeedCluster());

            ClusterAssociationMap::const_iterator iterOuter = transverseAssociationMap.find(pOuterCluster);
            if (transverseAssociationMap.end() == iterOuter)
                continue;

            if (pInnerCluster == pOuterCluster)
                continue;

            if (iterInner->second.m_forwardAssociations.count(pOuterCluster) == 0 || iterOuter->second.m_backwardAssociations.count(pInnerCluster) == 0)
                continue;

            if (!this->IsExtremalCluster(true, pInnerCluster, pOuterCluster) || !this->IsExtremalCluster(false, pOuterCluster, pInnerCluster))
                continue;

            if (!this->IsTransverseAssociated(pInnerTransverseCluster, pOuterTransverseCluster, nearbyClusters))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetAssociatedClusters(const ClusterToClustersMap &nearbyClusters, const Cluster *const pClusterI,
    const ClusterAssociationMap &associationMap, ClusterVector &associatedVector) const
{
    ClusterAssociationMap::const_iterator iterI = associationMap.find(pClusterI);
    if (associationMap.end() == iterI)
        return;

    for (ClusterSet::const_iterator iterJ = iterI->second.m_forwardAssociations.begin(), iterEndJ = iterI->second.m_forwardAssociations.end();
         iterJ != iterEndJ; ++iterJ)
    {
        const Cluster *const pClusterJ = *iterJ;

        if (this->IsTransverseAssociated(pClusterI, pClusterJ, nearbyClusters))
            associatedVector.push_back(pClusterJ);
    }

    for (ClusterSet::const_iterator iterJ = iterI->second.m_backwardAssociations.begin(), iterEndJ = iterI->second.m_backwardAssociations.end();
         iterJ != iterEndJ; ++iterJ)
    {
        const Cluster *const pClusterJ = *iterJ;

        if (this->IsTransverseAssociated(pClusterJ, pClusterI, nearbyClusters))
            associatedVector.push_back(pClusterJ);
    }

    std::sort(associatedVector.begin(), associatedVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsAssociated(const bool isForward, const Cluster *const pFirstCluster,
    const Cluster *const pSecondCluster, const ClusterToClustersMap &nearbyClusters) const
{
    if ((0 == nearbyClusters.at(pFirstCluster).count(pSecondCluster)) || (0 == nearbyClusters.at(pSecondCluster).count(pFirstCluster)))
        return false;

    CartesianVector firstInner(0.f, 0.f, 0.f), firstOuter(0.f, 0.f, 0.f);
    CartesianVector secondInner(0.f, 0.f, 0.f), secondOuter(0.f, 0.f, 0.f);
    this->GetExtremalCoordinatesX(pFirstCluster, firstInner, firstOuter);
    this->GetExtremalCoordinatesX(pSecondCluster, secondInner, secondOuter);

    const CartesianVector firstCoordinate(isForward ? firstOuter : firstInner);
    const CartesianVector secondCoordinate(isForward ? secondOuter : secondInner);

    if ((firstCoordinate.GetZ() > std::max(secondInner.GetZ(), secondOuter.GetZ()) + m_maxLongitudinalOverlap) ||
        (firstCoordinate.GetZ() < std::min(secondInner.GetZ(), secondOuter.GetZ()) - m_maxLongitudinalOverlap))
        return false;

    if ((isForward && secondCoordinate.GetX() < firstCoordinate.GetX()) || (!isForward && secondCoordinate.GetX() > firstCoordinate.GetX()))
        return false;

    const CartesianVector firstProjection(LArClusterHelper::GetClosestPosition(firstCoordinate, pSecondCluster));

    if ((isForward && firstProjection.GetX() < firstCoordinate.GetX() - m_maxTransverseOverlap) ||
        (!isForward && firstProjection.GetX() > firstCoordinate.GetX() + m_maxTransverseOverlap))
        return false;

    if ((isForward && firstProjection.GetX() > firstCoordinate.GetX() + m_clusterWindow) ||
        (!isForward && firstProjection.GetX() < firstCoordinate.GetX() - m_clusterWindow))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(
    const Cluster *const pInnerCluster, const Cluster *const pOuterCluster, const ClusterToClustersMap &nearbyClusters) const
{
    if ((0 == nearbyClusters.at(pInnerCluster).count(pOuterCluster)) || (0 == nearbyClusters.at(pOuterCluster).count(pInnerCluster)))
        return false;

    CartesianVector innerInner(0.f, 0.f, 0.f), innerOuter(0.f, 0.f, 0.f);
    CartesianVector outerInner(0.f, 0.f, 0.f), outerOuter(0.f, 0.f, 0.f);
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

    if (innerOuter.GetX() > innerProjection.GetX() + m_maxTransverseOverlap || outerInner.GetX() < outerProjection.GetX() - m_maxTransverseOverlap)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const LArTransverseCluster *const pInnerTransverseCluster,
    const LArTransverseCluster *const pOuterTransverseCluster, const ClusterToClustersMap &nearbyClusters) const
{
    if (pInnerTransverseCluster->GetDirection().GetDotProduct(pOuterTransverseCluster->GetDirection()) < m_transverseClusterMinCosTheta)
        return false;

    if (!this->IsTransverseAssociated(pInnerTransverseCluster, pOuterTransverseCluster->GetInnerVertex()))
        return false;

    if (!this->IsTransverseAssociated(pOuterTransverseCluster, pInnerTransverseCluster->GetOuterVertex()))
        return false;

    if (!this->IsTransverseAssociated(pInnerTransverseCluster->GetSeedCluster(), pOuterTransverseCluster->GetSeedCluster(), nearbyClusters))
        return false;

    if (this->IsOverlapping(pInnerTransverseCluster->GetSeedCluster(), pOuterTransverseCluster->GetSeedCluster()))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const LArTransverseCluster *const pTransverseCluster, const CartesianVector &testVertex) const
{
    const CartesianVector &innerVertex(pTransverseCluster->GetInnerVertex());
    const CartesianVector &outerVertex(pTransverseCluster->GetOuterVertex());
    const CartesianVector &direction(pTransverseCluster->GetDirection());

    const HitType view{LArClusterHelper::GetClusterHitType(pTransverseCluster->GetSeedCluster())};
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), view)};
    const float transverseMaxDisplacementAdjusted{ratio * m_transverseClusterMaxDisplacement};
    const float transverseMaxDisplacementSquaredAdjusted{transverseMaxDisplacementAdjusted * transverseMaxDisplacementAdjusted};

    if (direction.GetCrossProduct(testVertex - innerVertex).GetMagnitudeSquared() > transverseMaxDisplacementSquaredAdjusted)
        return false;

    if ((direction.GetDotProduct(testVertex - innerVertex) < -1.f * m_clusterWindow) ||
        (direction.GetDotProduct(testVertex - outerVertex) > +1.f * m_clusterWindow))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsOverlapping(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
{
    CartesianVector innerInner(0.f, 0.f, 0.f), innerOuter(0.f, 0.f, 0.f);
    CartesianVector outerInner(0.f, 0.f, 0.f), outerOuter(0.f, 0.f, 0.f);
    this->GetExtremalCoordinatesX(pInnerCluster, innerInner, innerOuter);
    this->GetExtremalCoordinatesX(pOuterCluster, outerInner, outerOuter);

    const CartesianVector innerProjection(LArClusterHelper::GetClosestPosition(outerInner, pInnerCluster));
    const CartesianVector outerProjection(LArClusterHelper::GetClosestPosition(innerOuter, pOuterCluster));

    const float innerOverlapSquared((innerProjection - innerOuter).GetMagnitudeSquared());
    const float outerOverlapSquared((outerProjection - outerInner).GetMagnitudeSquared());

    return (std::max(innerOverlapSquared, outerOverlapSquared) > m_maxProjectedOverlap * m_maxProjectedOverlap);
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
        const Cluster *const pAssociatedCluster = *iter;

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
        if (std::fabs(testMaxX - currentMaxX) > std::numeric_limits<float>::epsilon())
            return (testMaxX > currentMaxX);
    }
    else
    {
        if (std::fabs(testMinX - currentMaxX) > std::numeric_limits<float>::epsilon())
            return (testMinX < currentMinX);
    }

    return LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster);
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

void TransverseAssociationAlgorithm::GetExtremalCoordinatesX(
    const Cluster *const pCluster, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate) const
{
    CartesianVector firstCoordinate(0.f, 0.f, 0.f), secondCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pCluster, firstCoordinate, secondCoordinate);

    innerCoordinate = (firstCoordinate.GetX() < secondCoordinate.GetX() ? firstCoordinate : secondCoordinate);
    outerCoordinate = (firstCoordinate.GetX() > secondCoordinate.GetX() ? firstCoordinate : secondCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillReducedAssociationMap(
    const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const
{
    return this->FillReducedAssociationMap(inputAssociationMap, inputAssociationMap, inputAssociationMap, outputAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillReducedAssociationMap(const ClusterAssociationMap &firstAssociationMap,
    const ClusterAssociationMap &secondAssociationMap, const ClusterAssociationMap &secondAssociationMapSwapped,
    ClusterAssociationMap &clusterAssociationMap) const
{
    // Remove associations A->B from the first association map
    // if A->C exists in the second map and C->B exists in the reversed second map

    // Method can also be accessed through FillReducedAssociationMap(input,output) method,
    // which will remove association A->B from the input map if an association A->C and C->B
    // already exists in the map.

    ClusterVector sortedClusters;
    for (const auto &mapEntry : firstAssociationMap)
        sortedClusters.push_back(mapEntry.first);
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : sortedClusters)
    {
        const ClusterAssociation &firstAssociation(firstAssociationMap.at(pCluster));

        ClusterVector sortedOuterClusters(firstAssociation.m_forwardAssociations.begin(), firstAssociation.m_forwardAssociations.end());
        std::sort(sortedOuterClusters.begin(), sortedOuterClusters.end(), LArClusterHelper::SortByNHits);

        ClusterVector sortedInnerClusters(firstAssociation.m_backwardAssociations.begin(), firstAssociation.m_backwardAssociations.end());
        std::sort(sortedInnerClusters.begin(), sortedInnerClusters.end(), LArClusterHelper::SortByNHits);

        ClusterAssociationMap::const_iterator iterSecond = secondAssociationMap.find(pCluster);
        ClusterVector sortedMiddleClustersF, sortedMiddleClustersB;

        if (secondAssociationMap.end() != iterSecond)
        {
            sortedMiddleClustersF.insert(sortedMiddleClustersF.end(), iterSecond->second.m_forwardAssociations.begin(),
                iterSecond->second.m_forwardAssociations.end());
            sortedMiddleClustersB.insert(sortedMiddleClustersB.end(), iterSecond->second.m_backwardAssociations.begin(),
                iterSecond->second.m_backwardAssociations.end());
            std::sort(sortedMiddleClustersF.begin(), sortedMiddleClustersF.end(), LArClusterHelper::SortByNHits);
            std::sort(sortedMiddleClustersB.begin(), sortedMiddleClustersB.end(), LArClusterHelper::SortByNHits);
        }

        for (const Cluster *const pOuterCluster : sortedOuterClusters)
        {
            bool isNeighbouringCluster(true);

            for (const Cluster *const pMiddleCluster : sortedMiddleClustersF)
            {
                ClusterAssociationMap::const_iterator iterSecondCheck = secondAssociationMapSwapped.find(pMiddleCluster);
                if (secondAssociationMapSwapped.end() == iterSecondCheck)
                    continue;

                if (iterSecondCheck->second.m_forwardAssociations.count(pOuterCluster) > 0)
                {
                    isNeighbouringCluster = false;
                    break;
                }
            }

            if (isNeighbouringCluster)
                clusterAssociationMap[pCluster].m_forwardAssociations.insert(pOuterCluster);
        }

        for (const Cluster *const pInnerCluster : sortedInnerClusters)
        {
            bool isNeighbouringCluster(true);

            for (const Cluster *const pMiddleCluster : sortedMiddleClustersB)
            {
                ClusterAssociationMap::const_iterator iterSecondCheck = secondAssociationMapSwapped.find(pMiddleCluster);
                if (secondAssociationMapSwapped.end() == iterSecondCheck)
                    continue;

                if (iterSecondCheck->second.m_backwardAssociations.count(pInnerCluster) > 0)
                {
                    isNeighbouringCluster = false;
                    break;
                }
            }

            if (isNeighbouringCluster)
                clusterAssociationMap[pCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillSymmetricAssociationMap(
    const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const
{
    // Generate a symmetrised association map, so that both A--Fwd-->B and B--Bwd-->A both exist.
    // If A is associated to B through both a backward and forward association (very bad!),
    // try to rationalise this through majority voting, otherwise remove the association.

    ClusterVector sortedClusters;
    for (const auto &mapEntry : inputAssociationMap)
        sortedClusters.push_back(mapEntry.first);
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : sortedClusters)
    {
        const ClusterAssociation &inputAssociation(inputAssociationMap.at(pCluster));

        ClusterVector sortedForwardClusters(inputAssociation.m_forwardAssociations.begin(), inputAssociation.m_forwardAssociations.end());
        std::sort(sortedForwardClusters.begin(), sortedForwardClusters.end(), LArClusterHelper::SortByNHits);

        ClusterVector sortedBackwardClusters(inputAssociation.m_backwardAssociations.begin(), inputAssociation.m_backwardAssociations.end());
        std::sort(sortedBackwardClusters.begin(), sortedBackwardClusters.end(), LArClusterHelper::SortByNHits);

        // Symmetrise forward associations
        for (const Cluster *const pForwardCluster : sortedForwardClusters)
        {
            int nCounter(+1);

            if (inputAssociation.m_backwardAssociations.count(pForwardCluster))
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
                if (!(outputAssociationMap[pCluster].m_backwardAssociations.count(pForwardCluster) == 0 &&
                        outputAssociationMap[pForwardCluster].m_forwardAssociations.count(pCluster) == 0))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                outputAssociationMap[pCluster].m_forwardAssociations.insert(pForwardCluster);
                outputAssociationMap[pForwardCluster].m_backwardAssociations.insert(pCluster);
            }
        }

        // Symmetrise backward associations
        for (const Cluster *const pBackwardCluster : sortedBackwardClusters)
        {
            int nCounter(-1);

            if (inputAssociation.m_forwardAssociations.count(pBackwardCluster))
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
                if (!(outputAssociationMap[pCluster].m_forwardAssociations.count(pBackwardCluster) == 0 &&
                        outputAssociationMap[pBackwardCluster].m_backwardAssociations.count(pCluster) == 0))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                outputAssociationMap[pCluster].m_backwardAssociations.insert(pBackwardCluster);
                outputAssociationMap[pBackwardCluster].m_forwardAssociations.insert(pCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FinalizeClusterAssociationMap(
    const ClusterAssociationMap &inputAssociationMap, ClusterAssociationMap &outputAssociationMap) const
{
    ClusterAssociationMap intermediateAssociationMap;
    this->FillSymmetricAssociationMap(inputAssociationMap, intermediateAssociationMap);
    this->FillReducedAssociationMap(intermediateAssociationMap, outputAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TransverseAssociationAlgorithm::LArTransverseCluster::LArTransverseCluster(const Cluster *const pSeedCluster, const ClusterVector &associatedClusters) :
    m_pSeedCluster(pSeedCluster),
    m_associatedClusters(associatedClusters),
    m_innerVertex(0.f, 0.f, 0.f),
    m_outerVertex(0.f, 0.f, 0.f),
    m_direction(0.f, 0.f, 0.f)
{
    double Swxx(0.), Swzx(0.), Swz(0.), Swx(0.), Sw(0.);
    double minX(std::numeric_limits<double>::max());
    double maxX(-std::numeric_limits<double>::max());

    ClusterList clusterList(1, pSeedCluster);
    clusterList.insert(clusterList.end(), associatedClusters.begin(), associatedClusters.end());

    for (ClusterList::const_iterator iterI = clusterList.begin(), iterEndI = clusterList.end(); iterI != iterEndI; ++iterI)
    {
        for (OrderedCaloHitList::const_iterator iterJ = (*iterI)->GetOrderedCaloHitList().begin(),
                                                iterEndJ = (*iterI)->GetOrderedCaloHitList().end();
             iterJ != iterEndJ; ++iterJ)
        {
            for (CaloHitList::const_iterator iterK = iterJ->second->begin(), iterEndK = iterJ->second->end(); iterK != iterEndK; ++iterK)
            {
                const CaloHit *const pCaloHit = *iterK;

                if (pCaloHit->GetPositionVector().GetX() < minX)
                    minX = pCaloHit->GetPositionVector().GetX();

                if (pCaloHit->GetPositionVector().GetX() > maxX)
                    maxX = pCaloHit->GetPositionVector().GetX();

                Swxx += pCaloHit->GetPositionVector().GetX() * pCaloHit->GetPositionVector().GetX();
                Swzx += pCaloHit->GetPositionVector().GetZ() * pCaloHit->GetPositionVector().GetX();
                Swz += pCaloHit->GetPositionVector().GetZ();
                Swx += pCaloHit->GetPositionVector().GetX();
                Sw += 1.;
            }
        }
    }

    if (Sw > 0.f)
    {
        const double averageX(Swx / Sw);
        const double averageZ(Swz / Sw);

        if (Sw * Swxx - Swx * Swx > 0.)
        {
            double m((Sw * Swzx - Swx * Swz) / (Sw * Swxx - Swx * Swx));
            double px(1. / std::sqrt(1. + m * m));
            double pz(m / std::sqrt(1. + m * m));

            m_innerVertex.SetValues(static_cast<float>(minX), 0.f, static_cast<float>(averageZ + m * (minX - averageX)));
            m_outerVertex.SetValues(static_cast<float>(maxX), 0.f, static_cast<float>(averageZ + m * (maxX - averageX)));
            m_direction.SetValues(static_cast<float>(px), 0.f, static_cast<float>(pz));
        }
        else
        {
            m_innerVertex.SetValues(static_cast<float>(averageX), 0.f, static_cast<float>(averageZ));
            m_outerVertex.SetValues(static_cast<float>(averageX), 0.f, static_cast<float>(averageZ));
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FirstLengthCut", m_firstLengthCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SecondLengthCut", m_secondLengthCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterWindow", m_clusterWindow));

    const StatusCode angleStatusCode(XmlHelper::ReadValue(xmlHandle, "clusterAngle", m_clusterAngle));

    if (STATUS_CODE_SUCCESS == angleStatusCode)
    {
        m_clusterCosAngle = std::cos(m_clusterAngle * M_PI / 180.f);
        m_clusterTanAngle = std::tan(m_clusterAngle * M_PI / 180.f);
    }
    else if (STATUS_CODE_NOT_FOUND != angleStatusCode)
    {
        return angleStatusCode;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxTransverseOverlap", m_maxTransverseOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxProjectedOverlap", m_maxProjectedOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxLongitudinalOverlap", m_maxLongitudinalOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TransverseClusterMinCosTheta", m_transverseClusterMinCosTheta));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TransverseClusterMinLength", m_transverseClusterMinLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TransverseClusterMaxDisplacement", m_transverseClusterMaxDisplacement));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
