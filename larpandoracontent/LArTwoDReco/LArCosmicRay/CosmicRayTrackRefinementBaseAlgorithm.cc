/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayTrackRefinementBaseAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray track refinement base class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayTrackRefinementBaseAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

CosmicRayTrackRefinementBaseAlgorithm::CosmicRayTrackRefinementBaseAlgorithm() :
    m_microSlidingFitWindow(10), //ISOBEL THIS WAS 20?
    m_macroSlidingFitWindow(1000),
    m_stableRegionClusterFraction(0.05),
    m_mergePointMinCosAngleDeviation(0.995),
    m_distanceFromLine(0.35f),
    m_minHitFractionForHitRemoval(0.05f),
    m_maxDistanceFromMainTrack(0.75f),
    m_maxHitDistanceFromCluster(4.f),
    m_maxHitSeparationForConnectedCluster(4.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool CosmicRayTrackRefinementBaseAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &clusterMicroFitResult, const TwoDSlidingFitResult &clusterMacroFitResult,
    const TwoDSlidingFitResult &associatedMacroFitResult, const bool isUpstream, CartesianVector &clusterMergePosition, CartesianVector &clusterMergeDirection) const
{
    CartesianVector clusterAverageDirection(0.f, 0.f, 0.f), associatedAverageDirection(0.f, 0.f, 0.f);
    clusterMacroFitResult.GetGlobalDirection(clusterMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), clusterAverageDirection);
    associatedMacroFitResult.GetGlobalDirection(associatedMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), associatedAverageDirection);

    const LayerFitResultMap &clusterMicroLayerFitResultMap(clusterMicroFitResult.GetLayerFitResultMap());
    const int startLayer(isUpstream ? clusterMicroFitResult.GetMaxLayer() : clusterMicroFitResult.GetMinLayer());
    const int endLayer(isUpstream ? clusterMicroFitResult.GetMinLayer() : clusterMicroFitResult.GetMaxLayer());
    const int loopTerminationLayer(endLayer + (isUpstream ? -1 : 1));
    const int step(isUpstream ? -1 : 1);

    // ATTN: Search for stable region for which the local layer gradient agrees well with associated cluster global gradient
    unsigned int gradientStabilityWindow(std::ceil(clusterMicroLayerFitResultMap.size() *  m_stableRegionClusterFraction));
    unsigned int goodLayerCount(0);

    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        clusterMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);

        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(associatedAverageDirection));
        if (cosDirectionOpeningAngle > m_mergePointMinCosAngleDeviation)
        {
            // ATTN: Set merge point and direction as that at the first layer in the stable region
            if (goodLayerCount == 0)
            {
                // ATTN: Cluster direction vectors must point to one another
                clusterMergeDirection = clusterAverageDirection * (isUpstream ? 1.f : -1.f);
                clusterMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), clusterMergePosition);
            }

            ++goodLayerCount;
        }
        else
        {
            goodLayerCount = 0;
        }

        if (goodLayerCount > gradientStabilityWindow)
            break;

        // ATTN: Abort merging process have not found a stable region
        if (i == endLayer)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRefinementBaseAlgorithm::GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const ClusterList *const pClusterList,
    ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());
    const float minX(std::min(upstreamPoint.GetX(), downstreamPoint.GetX())), maxX(std::max(upstreamPoint.GetX(), downstreamPoint.GetX()));

    for (const Cluster *const pCluster : *pClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < upstreamPoint.GetZ()) || (hitPosition.GetZ() > downstreamPoint.GetZ()))
                    continue;

                const float distanceFromLine(clusterAssociation.GetConnectingLineDirection().GetCrossProduct(hitPosition - upstreamPoint).GetMagnitude());
                if (distanceFromLine > m_distanceFromLine)
                    continue;

                clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
            }
        }
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *CosmicRayTrackRefinementBaseAlgorithm::RemoveOffAxisHitsFromTrack(const Cluster *const pCluster, const CartesianVector &splitPosition, const bool isUpstream,
    const ClusterToCaloHitListMap &clusterToCaloHitListMap, ClusterList &remnantClusterList, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    float rL(0.f), rT(0.f);
    const TwoDSlidingFitResult &microFitResult(microSlidingFitResultMap.at(pCluster));
    microFitResult.GetLocalPosition(splitPosition, rL, rT);

    const TwoDSlidingFitResult macroFitResult(macroSlidingFitResultMap.at(pCluster));
    CartesianVector averageDirection(0.f, 0.f, 0.f);
    macroFitResult.GetGlobalDirection(macroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), averageDirection);

    const bool isVertical(std::fabs(averageDirection.GetX()) < std::numeric_limits<float>::epsilon());
    const float clusterGradient(isVertical ? 0.f : averageDirection.GetZ() / averageDirection.GetX());
    const float clusterIntercept(isVertical ? splitPosition.GetX() : splitPosition.GetZ() - (clusterGradient * splitPosition.GetX()));

    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pCluster);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pMainTrackCluster(nullptr), *pAboveCluster(nullptr), *pBelowCluster(nullptr);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            float thisL(0.f), thisT(0.f);

            microFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            bool isAnExtrapolatedHit(false);
            const auto extrapolatedCaloHitIter(clusterToCaloHitListMap.find(pCluster));
            
            if (extrapolatedCaloHitIter != clusterToCaloHitListMap.end())
                isAnExtrapolatedHit = std::find(extrapolatedCaloHitIter->second.begin(), extrapolatedCaloHitIter->second.end(), pCaloHit) != extrapolatedCaloHitIter->second.end();
            
            const bool isAbove(((clusterGradient * hitPosition.GetX()) + clusterIntercept) < (isVertical ? hitPosition.GetX() : hitPosition.GetZ()));
            const bool isToRemove(!isAnExtrapolatedHit && (((thisL > rL) && isUpstream) || ((thisL < rL) && !isUpstream)));

            const Cluster *&pClusterToModify(isToRemove ? (isAbove ? pAboveCluster : pBelowCluster) : pMainTrackCluster);

            if (pClusterToModify)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
            }
            else
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));

                if (pClusterToModify != pMainTrackCluster)
                    remnantClusterList.push_back(pClusterToModify);
            }
        }
    }

    //this->UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    // End fragmentation
    if (pAboveCluster || pBelowCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
        return pMainTrackCluster;
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        return pCluster;
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRefinementBaseAlgorithm::AddHitsToMainTrack(const Cluster *const pShowerCluster, const Cluster *const pMainTrackCluster, const CaloHitList &caloHitsToMerge,
    const ClusterAssociation &clusterAssociation, ClusterList &remnantClusterList) const
{
    // To ignore crossing CR muon or test beam tracks
    if ((static_cast<float>(caloHitsToMerge.size()) / static_cast<float>(pShowerCluster->GetNCaloHits())) < m_minHitFractionForHitRemoval)
        return;

    if (pShowerCluster->GetNCaloHits() == caloHitsToMerge.size())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pMainTrackCluster, pShowerCluster));
        return;
    }

    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pShowerCluster);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pAboveCluster(nullptr), *pBelowCluster(nullptr);

    const bool isVertical(std::fabs(clusterAssociation.GetConnectingLineDirection().GetX()) < std::numeric_limits<float>::epsilon());
    const float connectingLineGradient(isVertical ? 0.f : clusterAssociation.GetConnectingLineDirection().GetZ() / clusterAssociation.GetConnectingLineDirection().GetX());
    const float connectingLineIntercept(isVertical ? clusterAssociation.GetUpstreamMergePoint().GetX() :
        clusterAssociation.GetUpstreamMergePoint().GetZ() - (connectingLineGradient * clusterAssociation.GetUpstreamMergePoint().GetX()));

    const OrderedCaloHitList orderedCaloHitList(pShowerCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            const bool isAnExtrapolatedHit(std::find(caloHitsToMerge.begin(), caloHitsToMerge.end(), pCaloHit) != caloHitsToMerge.end());
            if (isAnExtrapolatedHit)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pShowerCluster, pCaloHit));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pMainTrackCluster, pCaloHit));
            }
            else
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                const bool isAbove(((connectingLineGradient * hitPosition.GetX()) + connectingLineIntercept) < (isVertical ? hitPosition.GetX() : hitPosition.GetZ()));
                const Cluster *&pClusterToModify(isAbove ? pAboveCluster : pBelowCluster);

                if (pClusterToModify)
                {
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
                }
                else
                {
                    PandoraContentApi::Cluster::Parameters parameters;
                    parameters.m_caloHitList.push_back(pCaloHit);
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));

                    remnantClusterList.push_back(pClusterToModify);
                }
            }
        }
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
}

 

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRefinementBaseAlgorithm::ProcessRemnantClusters(const ClusterList &remnantClusterList, const Cluster *const pMainTrackCluster, const ClusterList *const pClusterList, ClusterList &createdClusters) const
{
    ClusterList fragmentedClusters;
    for (const Cluster *const pRemnantCluster : remnantClusterList)
    {
        if (this->IsClusterRemnantDisconnected(pRemnantCluster)) //ISOBEL: USE CLUSTER HELPER FUNCTION INSTEAD.
        {
            this->FragmentRemnantCluster(pRemnantCluster, fragmentedClusters);
        }
        else
        {
            fragmentedClusters.push_back(pRemnantCluster);
        }
    }

    for (const Cluster *const pFragmentedCluster : fragmentedClusters)
    {
        if ((pFragmentedCluster->GetNCaloHits() == 1) && (this->AddToNearestCluster(pFragmentedCluster, pMainTrackCluster, pClusterList)))
            continue;
            
        createdClusters.push_back(pFragmentedCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTrackRefinementBaseAlgorithm::AddToNearestCluster(const Cluster *const pClusterToMerge, const Cluster *const pMainTrackCluster, const ClusterList *const pClusterList) const
{
    const Cluster *pClosestCluster(nullptr);
    float closestDistance(std::numeric_limits<float>::max());

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster == pClusterToMerge)
            continue;

        const float separationDistance(LArClusterHelper::GetClosestDistance(pClusterToMerge, pCluster));

        if (separationDistance < closestDistance)
        {
            if ((pCluster == pMainTrackCluster) && (separationDistance > m_maxDistanceFromMainTrack))
                continue;

            pClosestCluster = pCluster;
            closestDistance = separationDistance;
        }
    }

    if (closestDistance < m_maxHitDistanceFromCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClosestCluster, pClusterToMerge));
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTrackRefinementBaseAlgorithm::IsClusterRemnantDisconnected(const Cluster *const pRemnantCluster) const
{
    if (pRemnantCluster->GetNCaloHits() == 1)
        return false;

    const OrderedCaloHitList &orderedCaloHitList(pRemnantCluster->GetOrderedCaloHitList());
    const CaloHit *pPreviousCaloHit(orderedCaloHitList.begin()->second->front());

    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            if (pCaloHit == pPreviousCaloHit)
                continue;

            const float separationDistanceSquared(pCaloHit->GetPositionVector().GetDistanceSquared(pPreviousCaloHit->GetPositionVector()));

            if (separationDistanceSquared > (m_maxHitSeparationForConnectedCluster * m_maxHitSeparationForConnectedCluster))
                return true;

            pPreviousCaloHit = pCaloHit;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRefinementBaseAlgorithm::FragmentRemnantCluster(const Cluster *const pRemnantCluster, ClusterList &fragmentedClusterList) const
{
    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pRemnantCluster);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    ClusterList createdClusters;

    const OrderedCaloHitList &orderedCaloHitList(pRemnantCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            const Cluster *pClosestCluster(nullptr);

            if (!createdClusters.empty())
            {
                float closestDistance(std::numeric_limits<float>::max());

                for (const Cluster *const pCreatedCluster : createdClusters)
                {
                    const float separationDistance(LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pCreatedCluster));
                    if ((separationDistance < closestDistance) && (separationDistance < m_maxHitDistanceFromCluster))
                    {
                        closestDistance = separationDistance;
                        pClosestCluster = pCreatedCluster;
                    }
                }
            }

            if (pClosestCluster)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClosestCluster, pCaloHit));
            }
            else
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClosestCluster));
                createdClusters.push_back(pClosestCluster);
            }
        }
    }

    // End fragmentation
    if (createdClusters.empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        fragmentedClusterList.push_back(pRemnantCluster);
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
        fragmentedClusterList.insert(fragmentedClusterList.begin(), createdClusters.begin(), createdClusters.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRefinementBaseAlgorithm::InitialiseContainers(const ClusterList *pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHits)
            continue;
    
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

         try
         {
             const TwoDSlidingFitResult microSlidingFitResult(pCluster, m_microSlidingFitWindow, slidingFitPitch);
             const TwoDSlidingFitResult macroSlidingFitResult(pCluster, m_macroSlidingFitWindow, slidingFitPitch);

             slidingFitResultMapPair.first->insert(TwoDSlidingFitResultMap::value_type(pCluster, microSlidingFitResult));
             slidingFitResultMapPair.second->insert(TwoDSlidingFitResultMap::value_type(pCluster, macroSlidingFitResult));
             clusterVector.push_back(pCluster);
         }
         catch (const StatusCodeException &) {}
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayEndpointCorrectionAlgorithm::UpdateContainers(const ClusterVector &clustersToDelete, const ClusterList &clustersToAdd, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    for (const Cluster *const pClusterToDelete : clustersToDelete)
    {
        const TwoDSlidingFitResultMap::const_iterator microFitToDelete(slidingFitResultMapPair.first->find(pClusterToDelete));
        if (microFitToDelete != slidingFitResultMapPair.first->end())
            slidingFitResultMapPair.first->erase(microFitToDelete);

        const TwoDSlidingFitResultMap::const_iterator macroFitToDelete(slidingFitResultMapPair.second->find(pClusterToDelete));
        if (macroFitToDelete != slidingFitResultMapPair.second->end())
            slidingFitResultMapPair.second->erase(macroFitToDelete);        

        const ClusterVector::const_iterator clusterToDelete(std::find(clusterVector.begin(), clusterVector.end(), pClusterToDelete));
        if (clusterToDelete != clusterVector.end())
            clusterVector.erase(clusterToDelete);
    }

    this->InitialiseContainers(clustersToAdd, clusterVector, slidingFitResultMapVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackRefinementBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MicroSlidingFitWindow", m_microSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MacroSlidingFitWindow", m_macroSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "StableRegionClusterFraction", m_stableRegionClusterFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergePointMinCosAngleDeviation", m_mergePointMinCosAngleDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFromLine", m_distanceFromLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitFractionForHitRemoval", m_minHitFractionForHitRemoval));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromMainTrack", m_maxDistanceFromMainTrack));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitDistanceFromCluster", m_maxHitDistanceFromCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitSeparationForConnectedCluster", m_maxHitSeparationForConnectedCluster));    
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayTrackRefinementBaseAlgorithm::ClusterAssociation::ClusterAssociation() :
    m_pUpstreamCluster(nullptr),
    m_pDownstreamCluster(nullptr),
    m_upstreamMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_upstreamMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_downstreamMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_downstreamMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_connectingLineDirection(CartesianVector(0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayTrackRefinementBaseAlgorithm::ClusterAssociation::ClusterAssociation(const Cluster *const pUpstreamCluster, const Cluster *const pDownstreamCluster,const CartesianVector &upstreamMergePoint,
        const CartesianVector &upstreamMergeDirection, const CartesianVector &downstreamMergePoint, const CartesianVector &downstreamMergeDirection) :
    m_pUpstreamCluster(pUpstreamCluster),
    m_pDownstreamCluster(pDownstreamCluster),
    m_upstreamMergePoint(upstreamMergePoint),
    m_upstreamMergeDirection(upstreamMergeDirection),
    m_downstreamMergePoint(downstreamMergePoint),
    m_downstreamMergeDirection(downstreamMergeDirection),
    m_connectingLineDirection(0.f, 0.f, 0.f)
{
    const CartesianVector connectingLineDirection(m_downstreamMergePoint.GetX() - m_upstreamMergePoint.GetX(), 0.f, m_downstreamMergePoint.GetZ() - m_upstreamMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

    

} // namespace lar_content
