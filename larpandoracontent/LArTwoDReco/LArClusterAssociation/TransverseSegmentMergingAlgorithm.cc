/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseSegmentMergingAlgorithm.cc
 *
 *  @brief  Implementation of the transverse segment merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/TransverseSegmentMergingAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

TransverseSegmentMergingAlgorithm::TransverseSegmentMergingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseSegmentMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    if (pClusterList->empty())
        return;
    const HitType view{LArClusterHelper::GetClusterHitType(pClusterList->front())};
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), view)};
    const float maxLength{ratio * 5.f};

    clusterVector.clear();
    for (const Cluster *pCluster : *pClusterList)
    {
        CartesianVector upstreamCoord(0, 0, 0), downstreamCoord(0, 0, 0);
        LArClusterHelper::GetExtremalCoordinates(pCluster, upstreamCoord, downstreamCoord);
        const float length{(downstreamCoord - upstreamCoord).GetMagnitude()};
        const float xSpan2{2 * std::abs(downstreamCoord.GetX() - upstreamCoord.GetX())};
        const float zSpan{std::abs(downstreamCoord.GetZ() - upstreamCoord.GetZ())};
        if (length < maxLength && (xSpan2 < zSpan))
            clusterVector.emplace_back(pCluster);
    }
    std::sort(clusterVector.begin(), clusterVector.end(), [](const Cluster *pCluster1, const Cluster *pCluster2) {
        const float x1{std::min(pCluster1->GetCentroid(pCluster1->GetInnerPseudoLayer()).GetX(),
            pCluster1->GetCentroid(pCluster1->GetOuterPseudoLayer()).GetX())};
        const float x2{std::min(pCluster2->GetCentroid(pCluster2->GetInnerPseudoLayer()).GetX(),
            pCluster2->GetCentroid(pCluster2->GetOuterPseudoLayer()).GetX())};

        return x1 <= x2; 
    });

    ClusterList clusterList;
    clusterList.insert(clusterList.end(), clusterVector.begin(), clusterVector.end());

    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "Clusters", RED));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseSegmentMergingAlgorithm::PopulateClusterAssociationMap(const ClusterVector &allClusters, ClusterAssociationMap &/*clusterAssociationMap*/) const
{
    ClusterToClustersMap nearbyClusters;
    this->GetNearbyClusterMap(allClusters, nearbyClusters);

    ClusterUseMap usedClusters;
    for (const Cluster *const pSeedCluster : allClusters)
    {
        if (usedClusters.find(pSeedCluster) != usedClusters.end())
            continue;
        ClusterSet associatedClusters;
        this->GetIndirectAssociations(pSeedCluster, nearbyClusters, usedClusters, associatedClusters);

        if (associatedClusters.empty())
            continue;
        ClusterList associatedList;
        associatedList.insert(associatedList.end(), associatedClusters.begin(), associatedClusters.end());

        // Currently assumnig all associations should be merged - this is not what we want ultimately
        // Then perform a weight PCA through each cluster's inner pseudo layers and outer pseudo layers
        // (i think this gives you the main direction and then two rough boundaries
        CaloHitList caloHits;
        for (const Cluster *const pCluster : associatedList)
        {
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
            const CaloHitList &isolatedHits(pCluster->GetIsolatedCaloHitList());
            caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
            usedClusters[pCluster] = true;
        }
        LArPcaHelper::WeightedPointVector weightedHits;
        for (const CaloHit *const pCaloHit : caloHits)
            weightedHits.emplace_back(std::make_pair(pCaloHit->GetPositionVector(), pCaloHit->GetInputEnergy()));

        // Spine
        CartesianVector centroid(0, 0, 0);
        LArPcaHelper::EigenValues pcaEigenValues(0, 0, 0);
        LArPcaHelper::EigenVectors pcaAxes;
        LArPcaHelper::RunPca(weightedHits, centroid, pcaEigenValues, pcaAxes);
        const Cluster *const pLeft{associatedList.front()};
        const Cluster *const pRight{associatedList.back()};
        const float length{(pLeft->GetCentroid(pLeft->GetInnerPseudoLayer()) - pRight->GetCentroid(pRight->GetInnerPseudoLayer())).GetMagnitude()};
        CartesianVector start(centroid - pcaAxes[0] * 0.5f * length);
        CartesianVector finish(centroid + pcaAxes[0] * 0.5f * length);
        CartesianVector startT(centroid - pcaAxes[1] * 0.5f * length);
        CartesianVector finishT(centroid + pcaAxes[1] * 0.5f * length);

        LArPcaHelper::WeightedPointVector weightedInner;
        LArPcaHelper::WeightedPointVector weightedOuter;
        CaloHitList innerHits, outerHits;
        for (const Cluster *const pCluster : associatedList)
        {
            CaloHitList localHits;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(localHits);
            const CaloHitList &isolatedHits(pCluster->GetIsolatedCaloHitList());
            localHits.insert(localHits.end(), isolatedHits.begin(), isolatedHits.end());
            const CaloHit *innerHit{nullptr}, *outerHit{nullptr};
            float innerProj{std::numeric_limits<float>::max()}, outerProj{-std::numeric_limits<float>::max()};
            for (const CaloHit *const pCaloHit : localHits)
            {
                CartesianVector hitVec{pCaloHit->GetPositionVector() - centroid};
                const float projection{pcaAxes[1].GetDotProduct(hitVec)};
                if (projection < innerProj)
                {
                    innerProj = projection;
                    innerHit = pCaloHit;
                }
                if (projection > outerProj)
                {
                    outerProj = projection;
                    outerHit = pCaloHit;
                }
            }
            // Ensure the hits are the right side of the spine axis
            if (innerHit && innerProj < 0.f)
            {
                weightedInner.emplace_back(std::make_pair(innerHit->GetPositionVector(), innerHit->GetInputEnergy()));
                innerHits.emplace_back(innerHit);
            }
            if (outerHit && outerProj > 0.f)
            {
                weightedOuter.emplace_back(std::make_pair(outerHit->GetPositionVector(), outerHit->GetInputEnergy()));
                outerHits.emplace_back(outerHit);
            }
        }

        // Inner
        CartesianVector centroidInner(0, 0, 0);
        LArPcaHelper::EigenValues pcaEigenValuesInner(0, 0, 0);
        LArPcaHelper::EigenVectors pcaAxesInner;
        LArPcaHelper::RunPca(weightedInner, centroidInner, pcaEigenValuesInner, pcaAxesInner);
        CartesianVector startInner(centroidInner - pcaAxesInner[0] * 0.5f * 30);
        CartesianVector finishInner(centroidInner + pcaAxesInner[0] * 0.5f * 30);

        // Outer
        CartesianVector centroidOuter(0, 0, 0);
        LArPcaHelper::EigenValues pcaEigenValuesOuter(0, 0, 0);
        LArPcaHelper::EigenVectors pcaAxesOuter;
        LArPcaHelper::RunPca(weightedOuter, centroidOuter, pcaEigenValuesOuter, pcaAxesOuter);
        CartesianVector startOuter(centroidOuter - pcaAxesOuter[0] * 0.5f * 30);
        CartesianVector finishOuter(centroidOuter + pcaAxesOuter[0] * 0.5f * 30);

        // Check which clusters have, say 2/3 of their hits inside the envelope defined by the inner and outer fits

        ClusterList mergedList;
        mergedList.insert(mergedList.end(), associatedList.begin(), associatedList.end());

        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &mergedList, "Merged", BLUE));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &finish, "Axis", RED, 5, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startT, &finishT, "T", BLACK, 5, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startInner, &finishInner, "Inner", GREEN, 5, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startOuter, &finishOuter, "Outer", GREEN, 5, 1));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &innerHits, "Inner Hits", BLACK));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &outerHits, "Outer Hits", ORANGE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    /*
            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
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
        this->FillTransverseSegmentMergingMap(nearbyClusters, transverseClusterList, secondAssociationMap, transverseAssociationMap);

        // Step 6: Finalise the forward/backward transverse associations by symmetrising the
        //         transverse association map and removing any double-counting
        this->FinalizeClusterAssociationMap(transverseAssociationMap, clusterAssociationMap);
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "TransverseSegmentMergingAlgorithm: exception " << statusCodeException.ToString() << std::endl;
    }

    for (TransverseClusterList::const_iterator iter = transverseClusterList.begin(), iterEnd = transverseClusterList.end(); iter != iterEnd; ++iter)
    {
        delete *iter;
    }*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseSegmentMergingAlgorithm::IsExtremalCluster(const bool /*isForward*/, const Cluster *const /*pCurrentCluster*/, const Cluster *const /*pTestCluster*/) const
{
/*    float currentMinX(0.f), currentMaxX(0.f);
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

    return LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster);*/
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseSegmentMergingAlgorithm::GetNearbyClusterMap(const ClusterVector &allClusters, ClusterToClustersMap &nearbyClusters) const
{
    const float drMax2{100.f}, dzMax{std::sqrt(33.f)};

    for (auto iter1 = allClusters.begin(); iter1 != allClusters.end(); ++iter1)
    {
        const Cluster *const pCluster1{*iter1};
        CartesianVector cluster1Pt1(0, 0, 0), cluster1Pt2(0, 0, 0);
        LArClusterHelper::GetExtremalCoordinates(pCluster1, cluster1Pt1, cluster1Pt2);
        for (auto iter2 = std::next(iter1); iter2 != allClusters.end(); ++iter2)
        {
            const Cluster *const pCluster2{*iter2};
            CartesianVector cluster2Pt1(0, 0, 0), cluster2Pt2(0, 0, 0);
            LArClusterHelper::GetExtremalCoordinates(pCluster2, cluster2Pt1, cluster2Pt2);

            bool near{false};
            const CartesianVector &diff1121{cluster1Pt1 - cluster2Pt1};
            const CartesianVector &diff1122{cluster1Pt1 - cluster2Pt2};
            const CartesianVector &diff1221{cluster1Pt2 - cluster2Pt1};
            const CartesianVector &diff1222{cluster1Pt2 - cluster2Pt2};
            if (diff1121.GetZ() < dzMax && diff1121.GetMagnitudeSquared() < drMax2)
                near = true;
            else if (diff1122.GetZ() < dzMax && diff1122.GetMagnitudeSquared() < drMax2)
                near = true;
            else if (diff1221.GetZ() < dzMax && diff1221.GetMagnitudeSquared() < drMax2)
                near = true;
            else if (diff1222.GetZ() < dzMax && diff1222.GetMagnitudeSquared() < drMax2)
                near = true;

            if (near)
            {
                nearbyClusters[pCluster1].insert(pCluster2);
                nearbyClusters[pCluster2].insert(pCluster1);
            }
        }
    }

/*    for (const Cluster *const pSeedCluster : allClusters)
    {
        if (nearbyClusters.find(pSeedCluster) == nearbyClusters.end())
            continue;
        const ClusterSet &assocClusters{nearbyClusters.at(pSeedCluster)};
        ClusterList seedList({pSeedCluster});
        ClusterList nearbyList;
        nearbyList.insert(nearbyList.end(), assocClusters.begin(), assocClusters.end());
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &seedList, "Seed", GREEN));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &nearbyList, "Near", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseSegmentMergingAlgorithm::GetIndirectAssociations(const Cluster *const pSeedCluster, const ClusterToClustersMap &nearbyClusters,
    ClusterUseMap &usedClusters, ClusterSet &associatedClusters) const
{
    /*
      Hoover up clusters based on chains of proximity. The try to get an overall direction fit through the collection of hits and also try
      to get "upper" and "lower" fits through the extremal points to create an appromxate envelope that can be used to decide what should
      be added. If we sort the clusters in x before we start, we know the seed cluster will be the leftmost point and so can be used as a kind
      of vertex (it might be an "endpoint" if the shower is actually propagating along negative x) which can be used to relate shower direction
      back to the neutrino vertex (we can similarly sort the clusters we've associated to find the other end)

      In this algorithm we're just getting the potential associations, the work of deciding if those associations are real will be done by
      another function.

      Need to ensure the original attempts to call this function are prevented if the seed is already used. Within this function the inner
      loop check prevents subsequent calls
     */

    if (nearbyClusters.find(pSeedCluster) == nearbyClusters.end())
        return;
    associatedClusters.insert(pSeedCluster);
    for (const Cluster *const pCluster : nearbyClusters.at(pSeedCluster))
    {
        if (usedClusters.find(pCluster) != usedClusters.end() || associatedClusters.find(pCluster) != associatedClusters.end())
            continue;
        associatedClusters.insert(pCluster);
        this->GetIndirectAssociations(pCluster, nearbyClusters, usedClusters, associatedClusters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseSegmentMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
