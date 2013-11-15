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
    clusterVector.clear();
    clusterVector.insert(clusterVector.begin(), pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &inputClusters, ClusterAssociationMap &clusterAssociationMap) const
{
    TransverseClusterList transverseClusterList;

    try
    {
        ClusterVector transverseClusters, longitudinalClusters;
        this->GetTransverseClusters(inputClusters, transverseClusters, longitudinalClusters);

        ClusterVector seedClusters, nonSeedClusters;
        this->GetSeedClusters(transverseClusters, longitudinalClusters, seedClusters, nonSeedClusters);

        this->FillTransverseClusterList(seedClusters, transverseClusters, transverseClusterList);

        LArClusterMergeMap forwardMergeMap, backwardMergeMap;
        this->FillClusterMergeMaps(transverseClusterList, forwardMergeMap, backwardMergeMap);
        this->FillClusterAssociationMap(forwardMergeMap, backwardMergeMap, clusterAssociationMap);
        this->PruneAssociations(clusterAssociationMap);
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

void TransverseAssociationAlgorithm::GetTransverseClusters(const ClusterVector &inputVector, ClusterVector &transverseVector, ClusterVector &longitudinalVector) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (pCluster->GetNCaloHits() <= m_clusterCaloHits)
        {
            transverseVector.push_back(pCluster);
        }
        else
        {
            static const CartesianVector longitudinalDirection(0.f, 0.f, 1.f);
            const ClusterHelper::ClusterFitResult &clusterFitResult(pCluster->GetFitToAllHitsResult());

            if (clusterFitResult.IsFitSuccessful())
            {
                const CartesianVector &fitDirection(clusterFitResult.GetDirection());

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

void TransverseAssociationAlgorithm::GetSeedClusters(const ClusterVector &transverseClusters, const ClusterVector &longitudinalClusters,
    ClusterVector &seedVector, ClusterVector &nonSeedVector) const
{
    for (ClusterVector::const_iterator iterT = transverseClusters.begin(), iterTEnd = transverseClusters.end(); iterT != iterTEnd; ++iterT)
    {
        Cluster *pClusterT = *iterT;
        bool isSeedCluster(true);

        if (pClusterT->GetNCaloHits() <= m_clusterCaloHits)
        {
            CartesianVector positionT((pClusterT->GetCentroid(pClusterT->GetInnerPseudoLayer()) + pClusterT->GetCentroid(pClusterT->GetOuterPseudoLayer())) * 0.5f);

            for (ClusterVector::const_iterator iterL = longitudinalClusters.begin(), iterLEnd = longitudinalClusters.end(); iterL != iterLEnd; ++iterL)
            {
                Cluster *pClusterL = *iterL;

                if (LArClusterHelper::GetClosestDistance(positionT, pClusterL) < m_clusterWindow)
                {
                    isSeedCluster = false;
                    break;
                }
            }
        }

        if (isSeedCluster)
        {
            seedVector.push_back(pClusterT);
        }
        else
        {
            nonSeedVector.push_back(pClusterT);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillTransverseClusterList(const ClusterVector &seedClusters, const ClusterVector &transverseClusters,
    TransverseClusterList &transverseClusterList) const
{
    for (ClusterVector::const_iterator iter = seedClusters.begin(), iterEnd = seedClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (pCluster->GetNCaloHits() > m_clusterCaloHits)
        {
            transverseClusterList.push_back(new LArTransverseCluster(pCluster, ClusterVector()));
        }
        else
        {
            ClusterVector associatedClusters;
            this->GetTransverseAssociatedClusters(pCluster, transverseClusters, associatedClusters);
            transverseClusterList.push_back(new LArTransverseCluster(pCluster, associatedClusters));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::GetTransverseAssociatedClusters(Cluster *const pCluster, const ClusterVector &inputClusterVector,
    ClusterVector &outputClusterVector) const
{
    for (ClusterVector::const_iterator iter = inputClusterVector.begin(), iterEnd = inputClusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pNearbyCluster = *iter;

        if (pCluster == pNearbyCluster)
            continue;

        if (pNearbyCluster->GetNCaloHits() > m_clusterCaloHits)
            continue;

        if (this->IsTransverseAssociated(pCluster, pNearbyCluster))
            outputClusterVector.push_back(pNearbyCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    const CartesianVector innerCentroid1(pCluster1->GetCentroid(pCluster1->GetInnerPseudoLayer()));
    const CartesianVector outerCentroid1(pCluster1->GetCentroid(pCluster1->GetOuterPseudoLayer()));
    const CartesianVector innerCentroid2(pCluster2->GetCentroid(pCluster2->GetInnerPseudoLayer()));
    const CartesianVector outerCentroid2(pCluster2->GetCentroid(pCluster2->GetOuterPseudoLayer()));

    const float averageX1(0.5f * (innerCentroid1.GetX() + outerCentroid1.GetX()));
    const float averageZ1(0.5f * (innerCentroid1.GetZ() + outerCentroid1.GetZ()));
    const float averageX2(0.5f * (innerCentroid2.GetX() + outerCentroid2.GetX()));
    const float averageZ2(0.5f * (innerCentroid2.GetZ() + outerCentroid2.GetZ()));

    if ((std::fabs(averageX2 - averageX1) < m_clusterWindow) && (std::fabs(averageZ2 - averageZ1) < m_clusterWindow) &&
        (std::fabs(averageZ2 - averageZ1) < std::fabs(averageX2 - averageX1) * std::fabs(m_clusterTanAngle)))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillClusterMergeMaps(const TransverseClusterList &transverseClusterList, LArClusterMergeMap &forwardMergeMap,
    LArClusterMergeMap &backwardMergeMap) const
{
    for (TransverseClusterList::const_iterator iter = transverseClusterList.begin(), iterEnd = transverseClusterList.end(); iter != iterEnd; ++iter)
    {
        LArTransverseCluster *pTransverseCluster = *iter;
        Cluster *pCentralCluster(pTransverseCluster->GetSeedCluster());
        const ClusterVector &associatedClusters(pTransverseCluster->GetAssociatedClusters());

        for (ClusterVector::const_iterator cIter = associatedClusters.begin(), cIterEnd = associatedClusters.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;

            if (pCentralCluster == pCluster)
                continue;

            if (this->IsExtremalCluster(true, pCentralCluster, pCluster) && !this->IsExtremalCluster(false, pCentralCluster, pCluster))
            {
                forwardMergeMap[pCentralCluster].insert(pCluster);
                backwardMergeMap[pCluster].insert(pCentralCluster);
            }

            if (this->IsExtremalCluster(false, pCentralCluster, pCluster) && !this->IsExtremalCluster(true, pCentralCluster, pCluster))
            {
                forwardMergeMap[pCluster].insert(pCentralCluster);
                backwardMergeMap[pCentralCluster].insert(pCluster);
            }
        }
    }

    for (TransverseClusterList::const_iterator iter1 = transverseClusterList.begin(), iterEnd1 = transverseClusterList.end(); iter1 != iterEnd1; ++iter1)
    {
        LArTransverseCluster *pInnerTransverseCluster = *iter1;
        Cluster *pInnerCluster(pInnerTransverseCluster->GetSeedCluster());

        for (TransverseClusterList::const_iterator iter2 = transverseClusterList.begin(), iterEnd2 = transverseClusterList.end(); iter2 != iterEnd2; ++iter2)
        {
            LArTransverseCluster *pOuterTransverseCluster = *iter2;
            Cluster *pOuterCluster(pOuterTransverseCluster->GetSeedCluster());

            if (pInnerCluster == pOuterCluster)
                continue;

            if (this->IsExtremalCluster(true, pInnerCluster, pOuterCluster) && this->IsExtremalCluster(false, pOuterCluster, pInnerCluster))
            {
                if (this->IsTransverseAssociated(pInnerTransverseCluster, pOuterTransverseCluster))
                {
                    forwardMergeMap[pInnerCluster].insert(pOuterCluster);
                    backwardMergeMap[pOuterCluster].insert(pInnerCluster);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
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
    minX = std::numeric_limits<float>::max();
    maxX = -std::numeric_limits<float>::max();
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const float caloHitX((*hitIter)->GetPositionVector().GetX());

            if (caloHitX < minX)
                minX = caloHitX;

            if (caloHitX > maxX)
                maxX = caloHitX;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const LArTransverseCluster *const pTransverseCluster1,
    const LArTransverseCluster *const pTransverseCluster2) const
{
    if (pTransverseCluster1->GetDirection().GetDotProduct(pTransverseCluster2->GetDirection()) < m_minCosRelativeAngle)
        return false;

    if (!this->IsTransverseAssociated(pTransverseCluster1, pTransverseCluster2->GetInnerVertex()))
        return false;

    if (!this->IsTransverseAssociated(pTransverseCluster2, pTransverseCluster1->GetOuterVertex()))
        return false;

    if ((pTransverseCluster1->GetDirection().GetDotProduct(pTransverseCluster2->GetInnerVertex() - pTransverseCluster1->GetInnerVertex()) > 0.f) &&
        (pTransverseCluster1->GetDirection().GetDotProduct(pTransverseCluster2->GetOuterVertex() - pTransverseCluster1->GetOuterVertex()) > 0.f) &&
        (pTransverseCluster2->GetDirection().GetDotProduct(pTransverseCluster1->GetInnerVertex() - pTransverseCluster2->GetInnerVertex()) < 0.f) &&
        (pTransverseCluster2->GetDirection().GetDotProduct(pTransverseCluster1->GetOuterVertex() - pTransverseCluster2->GetOuterVertex()) < 0.f))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseAssociationAlgorithm::IsTransverseAssociated(const LArTransverseCluster *const pTransverseCluster, const pandora::CartesianVector& testVertex) const
{
    const CartesianVector &innerVertex(pTransverseCluster->GetInnerVertex());
    const CartesianVector &outerVertex(pTransverseCluster->GetOuterVertex());
    const CartesianVector &direction(pTransverseCluster->GetDirection());

    if (std::fabs(direction.GetCrossProduct(testVertex - innerVertex).GetMagnitudeSquared()) > m_maxTransverseSeparation * m_maxTransverseSeparation)
        return false;

    if ((direction.GetDotProduct(testVertex - innerVertex) < -m_maxLongitudinalSeparation) || (direction.GetDotProduct(testVertex - outerVertex) > +m_maxLongitudinalSeparation))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::FillClusterAssociationMap(LArClusterMergeMap &forwardMergeMap, LArClusterMergeMap &backwardMergeMap,
    ClusterAssociationMap &clusterAssociationMap) const
{
    // Select neighbouring forward associations
    for (LArClusterMergeMap::const_iterator iter1 = forwardMergeMap.begin(), iterEnd1 = forwardMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pInnerCluster = iter1->first;
        const ClusterList &clusterMerges(iter1->second);

        for (ClusterList::iterator iter2 = clusterMerges.begin(), iterEnd2 = clusterMerges.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster *pOuterCluster = *iter2;

            if (pOuterCluster == pInnerCluster)
                throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            bool isNeighbouringCluster(true);

            for (ClusterList::iterator iter3 = clusterMerges.begin(), iterEnd3 = clusterMerges.end(); iter3 != iterEnd3; ++iter3)
            {
                Cluster *pMiddleCluster = *iter3;

                if (pMiddleCluster == pOuterCluster)
                    continue;

                if (forwardMergeMap[pMiddleCluster].count(pOuterCluster))
                {
                    isNeighbouringCluster = false;
                    break;
                }
            }

            if (isNeighbouringCluster)
            {
                clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
                clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
            }
        }
    }

    // Select neighbouring backward associations
    for (LArClusterMergeMap::const_iterator iter1 = backwardMergeMap.begin(), iterEnd1 = backwardMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pOuterCluster = iter1->first;
        const ClusterList &clusterMerges(iter1->second);

        for (ClusterList::iterator iter2 = clusterMerges.begin(), iterEnd2 = clusterMerges.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster *pInnerCluster = *iter2;

            if (pInnerCluster == pOuterCluster)
                throw pandora::StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            bool isNeighbouringCluster(true);

            for (ClusterList::iterator iter3 = clusterMerges.begin(), iterEnd3 = clusterMerges.end(); iter3 != iterEnd3; ++iter3)
            {
                Cluster *pMiddleCluster = *iter3;

                if (pMiddleCluster == pInnerCluster)
                    continue;

                if (backwardMergeMap[pMiddleCluster].count(pInnerCluster))
                {
                    isNeighbouringCluster = false;
                    break;
                }
            }

            if (isNeighbouringCluster)
            {
                clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
                clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseAssociationAlgorithm::PruneAssociations(ClusterAssociationMap &clusterAssociationMap) const
{
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
                        ClusterList::iterator eraseIter = iter2->second.m_backwardAssociations.find(pCluster);

                        if (iter2->second.m_backwardAssociations.end() != eraseIter)
                            iter2->second.m_backwardAssociations.erase(eraseIter);

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
                        ClusterList::iterator eraseIter = iter2->second.m_forwardAssociations.find(pCluster);

                        if (iter2->second.m_forwardAssociations.end() != eraseIter)
                            iter2->second.m_forwardAssociations.erase(eraseIter);

                        carryOn = true;
                    }
                }
            }
        }
    }
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
