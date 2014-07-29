/**
 *  @file   LArContent/src/LArHelpers/LArClusterHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "Helpers/ClusterHelper.h"
#include "Helpers/XmlHelper.h"

#include "Pandora/PandoraSettings.h"

#include "LArCalculators/LArPseudoLayerCalculator.h"
#include "LArObjects/LArTwoDSlidingFitResult.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar
{

HitType LArClusterHelper::GetClusterHitType(const Cluster *const pCluster)
{
    if (0 == pCluster->GetNCaloHits())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (PandoraSettings::SingleHitTypeClusteringMode())
        return (*(pCluster->GetOrderedCaloHitList().begin()->second->begin()))->GetHitType();

    HitType hitType(CUSTOM);

    if (pCluster->ContainsHitType(TPC_VIEW_U))
    {
        if (CUSTOM == hitType) hitType = TPC_VIEW_U;
        else throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (pCluster->ContainsHitType(TPC_VIEW_V))
    {
        if (CUSTOM == hitType) hitType = TPC_VIEW_V;
        else throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (pCluster->ContainsHitType(TPC_VIEW_W))
    {
        if (CUSTOM == hitType) hitType = TPC_VIEW_W;
        else throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (CUSTOM == hitType)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return hitType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::LArTwoDSlidingFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, TwoDSlidingFitResult &twoDSlidingFitResult)
{
    // TODO: Look into more sophisticated treatment for sliding linear fit; lose layer structure and slide over pathlength-ordered hits
    CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster, innerCoordinate, outerCoordinate);

    LArClusterHelper::LArTwoDSlidingFit(pCluster, layerFitHalfWindow, innerCoordinate, (outerCoordinate - innerCoordinate).GetUnitVector(), twoDSlidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::LArTwoDSlidingXZFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, TwoDSlidingFitResult &twoDSlidingFitResult)
{
    const CartesianVector axisIntercept(0.f, 0.f, 0.f);
    const CartesianVector axisDirection(0.f, 0.f, 1.f);

    LArClusterHelper::LArTwoDSlidingFit(pCluster, layerFitHalfWindow, axisIntercept, axisDirection, twoDSlidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::LArTwoDSlidingFit(const Cluster *const pCluster, const unsigned int layerFitHalfWindow, const CartesianVector &axisIntercept,
    const CartesianVector &axisDirection, TwoDSlidingFitResult &twoDSlidingFitResult)
{
    if ((std::fabs(axisIntercept.GetY()) > std::numeric_limits<float>::epsilon()) || (std::fabs(axisDirection.GetY()) > std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    twoDSlidingFitResult.m_pCluster = pCluster;
    twoDSlidingFitResult.m_layerFitHalfWindow = layerFitHalfWindow;
    twoDSlidingFitResult.m_axisIntercept = axisIntercept;
    twoDSlidingFitResult.m_axisDirection = axisDirection;
    TwoDSlidingFitResult::LayerFitContributionMap &layerFitContributionMap(twoDSlidingFitResult.m_layerFitContributionMap);

    // Identify fit contributions
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            float rL(0.f), rT(0.f);
            twoDSlidingFitResult.GetLocalPosition((*hitIter)->GetPositionVector(), rL, rT);
            const int layer(twoDSlidingFitResult.GetLayer(rL));
            layerFitContributionMap[layer].AddPoint(rL, rT);
        }
    }

    LArClusterHelper::StoreSlidingFitResults(twoDSlidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::LArTwoDShowerEdgeFit(const TwoDSlidingFitResult &fullShowerFit, const ShowerEdge showerEdge, TwoDSlidingFitResult &twoDSlidingFitResult)
{
    twoDSlidingFitResult.m_pCluster = fullShowerFit.GetCluster();
    twoDSlidingFitResult.m_layerFitHalfWindow = fullShowerFit.GetLayerFitHalfWindow();
    twoDSlidingFitResult.m_axisIntercept = fullShowerFit.GetAxisIntercept();
    twoDSlidingFitResult.m_axisDirection = fullShowerFit.GetAxisDirection();
    TwoDSlidingFitResult::LayerFitContributionMap &layerFitContributionMap(twoDSlidingFitResult.m_layerFitContributionMap);

    // Examine all possible fit contributions
    typedef std::pair<float, float> FitCoordinate;
    typedef std::vector<FitCoordinate> FitCoordinateList;
    typedef std::map<int, FitCoordinateList> FitCoordinateMap;

    FitCoordinateMap fitCoordinateMap;
    const OrderedCaloHitList &orderedCaloHitList(fullShowerFit.GetCluster()->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            float rL(0.f), rT(0.f);
            twoDSlidingFitResult.GetLocalPosition((*hitIter)->GetPositionVector(), rL, rT);

            try
            {
                CartesianVector fullShowerFitPosition(0.f, 0.f, 0.f);
                fullShowerFit.GetGlobalFitPosition(rL, fullShowerFitPosition);

                float rLFit(0.f), rTFit(0.f);
                twoDSlidingFitResult.GetLocalPosition(fullShowerFitPosition, rLFit, rTFit);

                const float rTDiff(rT - rTFit);

                if (((POSITIVE_SHOWER_EDGE == showerEdge) && (rTDiff < 0.f)) || ((NEGATIVE_SHOWER_EDGE == showerEdge) && (rTDiff > 0.f)))
                    rT = rTFit;
            }
            catch (StatusCodeException &)
            {
                continue;
            }

            const int layer(twoDSlidingFitResult.GetLayer(rL));
            fitCoordinateMap[layer].push_back(FitCoordinate(rL, rT));
        }
    }

    // Select fit contributions representing relevant shower edge
    for (FitCoordinateMap::const_iterator iter = fitCoordinateMap.begin(), iterEnd = fitCoordinateMap.end(); iter != iterEnd; ++iter)
    {
        const int layer(iter->first);
        const FitCoordinateList &fitCoordinateList(iter->second);

        // TODO, improve this hit selection
        bool bestFitCoordinateFound(false);
        FitCoordinate bestFitCoordinate = (POSITIVE_SHOWER_EDGE == showerEdge) ?
            FitCoordinate(0.f, -std::numeric_limits<float>::max()) :
            FitCoordinate(0.f, +std::numeric_limits<float>::max());

        for (FitCoordinateList::const_iterator fIter = fitCoordinateList.begin(), fIterEnd = fitCoordinateList.end(); fIter != fIterEnd; ++fIter)
        {
            if (((POSITIVE_SHOWER_EDGE == showerEdge) && (fIter->second > bestFitCoordinate.second)) ||
                ((NEGATIVE_SHOWER_EDGE == showerEdge) && (fIter->second < bestFitCoordinate.second)))
            {
                bestFitCoordinate = *fIter;
                bestFitCoordinateFound = true;
            }
        }

        if (bestFitCoordinateFound)
            layerFitContributionMap[layer].AddPoint(bestFitCoordinate.first, bestFitCoordinate.second);
    }

    LArClusterHelper::StoreSlidingFitResults(twoDSlidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsMultivaluedInX(const TwoDSlidingFitResult &twoDSlidingFitResult)
{
    CartesianVector previousPosition(0.f, 0.f, 0.f);
    unsigned int nSteps(0), nPositiveSteps(0), nNegativeSteps(0), nUnchangedSteps(0);
    const TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap(twoDSlidingFitResult.GetLayerFitResultMap());

    for (TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        CartesianVector position(0.f, 0.f, 0.f);
        twoDSlidingFitResult.GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);

        const CartesianVector delta(position - previousPosition);
        previousPosition = position;
        ++nSteps;

        if (std::fabs(delta.GetX()) < std::fabs(delta.GetZ()) * m_multiValuedTanThetaCut)
        {
            ++nUnchangedSteps;
        }
        else if (delta.GetX() > 0.f)
        {
            ++nPositiveSteps;
        }
        else
        {
            ++nNegativeSteps;
        }
    }

    if (0 == nSteps)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    const float positiveStepFraction(static_cast<float>(nPositiveSteps) / static_cast<float>(nSteps));
    const float negativeStepFraction(static_cast<float>(nNegativeSteps) / static_cast<float>(nSteps));
    const float maxStepFraction(std::max(positiveStepFraction, negativeStepFraction));

    return (maxStepFraction < m_multiValuedStepFractionCut);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetSlidingFitWidth(const TwoDSlidingFitResult &twoDSlidingFitResult)
{
    FloatVector residuals;
    const OrderedCaloHitList &orderedCaloHitList(twoDSlidingFitResult.GetCluster()->GetOrderedCaloHitList());
    const TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap(twoDSlidingFitResult.GetLayerFitResultMap());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            float rL(0.f), rT(0.f);
            twoDSlidingFitResult.GetLocalPosition((*hitIter)->GetPositionVector(), rL, rT);
            const int layer(twoDSlidingFitResult.GetLayer(rL));

            TwoDSlidingFitResult::LayerFitResultMap::const_iterator fitResultIter = layerFitResultMap.find(layer);

            if (layerFitResultMap.end() == fitResultIter)
                continue;

            const double fitT(fitResultIter->second.GetFitT());
            const double gradient(fitResultIter->second.GetGradient());
            const double residualSquared((fitT - rT) * (fitT - rT) / (1. + gradient * gradient)); // angular correction (note: this is cheating!)
            residuals.push_back(residualSquared);
        }
    }

    if (residuals.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    std::sort(residuals.begin(), residuals.end());
    const float theQuantile(residuals[m_trackResidualQuantile * residuals.size()]);

    return std::sqrt(theQuantile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::LArTrackWidth(const Cluster *const pCluster)
{
    TwoDSlidingFitResult twoDSlidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pCluster, m_layerFitHalfWindow, twoDSlidingFitResult);
    return LArClusterHelper::GetSlidingFitWidth(twoDSlidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLengthSquared(const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // ATTN: in 2D case, we will actually calculate the quadrature sum of deltaX and deltaU/V/W
    float minX(std::numeric_limits<float>::max()), maxX(-std::numeric_limits<float>::max());
    float minY(std::numeric_limits<float>::max()), maxY(-std::numeric_limits<float>::max());
    float minZ(std::numeric_limits<float>::max()), maxZ(-std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CartesianVector &hitPosition((*hitIter)->GetPositionVector());
            minX = std::min(hitPosition.GetX(), minX);
            maxX = std::max(hitPosition.GetX(), maxX);
            minY = std::min(hitPosition.GetY(), minY);
            maxY = std::max(hitPosition.GetY(), maxY);
            minZ = std::min(hitPosition.GetZ(), minZ);
            maxZ = std::max(hitPosition.GetZ(), maxZ);
        }
    }

    const float deltaX(maxX - minX), deltaY(maxY - minY), deltaZ(maxZ - minZ);
    return (deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLength(const Cluster *const pCluster)
{
    return std::sqrt(LArClusterHelper::GetLengthSquared(pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetEnergyFromLength(const Cluster *const pCluster)
{
    static const float dEdX(0.002f); // approximately 2 MeV/cm

    return (dEdX * LArClusterHelper::GetLength(pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArClusterHelper::GetLayerSpan(const Cluster *const pCluster)
{
    return (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLayerOccupancy(const Cluster *const pCluster)
{
    const unsigned int nOccupiedLayers(pCluster->GetOrderedCaloHitList().size());
    const unsigned int nLayers(1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer());

    if (nLayers > 0)
        return (static_cast<float>(nOccupiedLayers) / static_cast<float>(nLayers));

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLayerOccupancy(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    const unsigned int nOccupiedLayers(pCluster1->GetOrderedCaloHitList().size() + pCluster2->GetOrderedCaloHitList().size());
    const unsigned int nLayers(1 + std::max(pCluster1->GetOuterPseudoLayer(), pCluster2->GetOuterPseudoLayer()) -
        std::min(pCluster1->GetInnerPseudoLayer(), pCluster2->GetInnerPseudoLayer()));

    if (nLayers > 0)
        return (static_cast<float>(nOccupiedLayers) / static_cast<float>(nLayers));

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const ClusterList &clusterList1, const ClusterList &clusterList2)
{
    if (clusterList1.empty() || clusterList2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float closestDistance(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iterEnd1 = clusterList1.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster *pCluster1 = *iter1;
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster1, clusterList2));
 
        if (thisDistance < closestDistance)
            closestDistance = thisDistance; 
    }

    return closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const Cluster *const pCluster, const ClusterList &clusterList)
{
    if (clusterList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float closestDistance(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *pTestCluster = *iter;
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster, pTestCluster));

        if (thisDistance < closestDistance)
            closestDistance = thisDistance; 
    }

    return closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    return ClusterHelper::GetDistanceToClosestHit(pCluster1, pCluster2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const CartesianVector &position, const Cluster *const pCluster)
{
    return (position - LArClusterHelper::GetClosestPosition(position, pCluster)).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::GetClosestPosition(const CartesianVector &position, const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    CaloHit *pClosestCaloHit(NULL);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit* pCaloHit = *hitIter;
            const float distanceSquared((pCaloHit->GetPositionVector() - position).GetMagnitudeSquared());

            if (distanceSquared < closestDistanceSquared)
            {
                closestDistanceSquared = distanceSquared;
                pClosestCaloHit = pCaloHit;
            }
        }
    }

    if (pClosestCaloHit)
        return pClosestCaloHit->GetPositionVector();

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClusterSpanXZ(const Cluster *const pCluster, CartesianVector &minimumCoordinate, CartesianVector &maximumCoordinate)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float xmin(std::numeric_limits<float>::max());
    float ymin(std::numeric_limits<float>::max());
    float zmin(std::numeric_limits<float>::max());
    float xmax(-std::numeric_limits<float>::max());
    float ymax(-std::numeric_limits<float>::max());
    float zmax(-std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *pCaloHit = *hIter;
            const CartesianVector &hit(pCaloHit->GetPositionVector());
            xmin = std::min(hit.GetX(), xmin);
            xmax = std::max(hit.GetX(), xmax);
            ymin = std::min(hit.GetY(), ymin);
            ymax = std::max(hit.GetY(), ymax);
            zmin = std::min(hit.GetZ(), zmin);
            zmax = std::max(hit.GetZ(), zmax);
        }
    }

    minimumCoordinate.SetValues(xmin, ymin, zmin);
    maximumCoordinate.SetValues(xmax, ymax, zmax);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClusterSpanX(const Cluster *const pCluster, float &xmin, float &xmax)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    xmin = std::numeric_limits<float>::max();
    xmax = -std::numeric_limits<float>::max();

    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *pCaloHit = *hIter;
            const CartesianVector &hit(pCaloHit->GetPositionVector());
            xmin = std::min(hit.GetX(), xmin);
            xmax = std::max(hit.GetX(), xmax);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClusterSpanZ(const Cluster *const pCluster, const float xmin, const float xmax, 
    float &zmin, float &zmax)
{
    if (xmin > xmax)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    zmin = std::numeric_limits<float>::max();
    zmax = -std::numeric_limits<float>::max();

    bool foundHits(false);
    
    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *pCaloHit = *hIter;
            const CartesianVector &hit(pCaloHit->GetPositionVector());

            if (hit.GetX() < xmin || hit.GetX() > xmax)
                continue;

            zmin = std::min(hit.GetZ(), zmin);
            zmax = std::max(hit.GetZ(), zmax);
            foundHits = true;
        }
    }

    if (!foundHits)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetAverageZ(const Cluster *const pCluster, const float xmin, const float xmax)
{
    if (xmin > xmax)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float zsum(0.f);
    int count(0);
    
    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *pCaloHit = *hIter;
            const CartesianVector &hit(pCaloHit->GetPositionVector());

            if (hit.GetX() < xmin || hit.GetX() > xmax)
                continue;

            zsum += hit.GetZ();
            ++count;
        }
    }

    if (count == 0)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return zsum / static_cast<float>(count);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinatesXZ(const Cluster *const pCluster, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate)
{
    CaloHitList candidateList;
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    // Transfer all inner layer hits
    OrderedCaloHitList::const_iterator iterInner = orderedCaloHitList.begin();
    for (CaloHitList::const_iterator hitIter = iterInner->second->begin(), hitIterEnd = iterInner->second->end(); hitIter != hitIterEnd; ++hitIter)
        candidateList.insert(*hitIter);

    // Transfer all outer layer hits
    OrderedCaloHitList::const_reverse_iterator iterOuter = orderedCaloHitList.rbegin();
    for (CaloHitList::const_iterator hitIter = iterOuter->second->begin(), hitIterEnd = iterOuter->second->end(); hitIter != hitIterEnd; ++hitIter)
        candidateList.insert(*hitIter);

    // Transfer the extremal hits in X (assume there are no ties)
    CaloHit *pFirstCaloHit(NULL);
    CaloHit *pSecondCaloHit(NULL);

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;

            if (NULL == pFirstCaloHit || pCaloHit->GetPositionVector().GetX() < pFirstCaloHit->GetPositionVector().GetX())
                pFirstCaloHit = pCaloHit;

            if (NULL == pSecondCaloHit || pCaloHit->GetPositionVector().GetX() > pSecondCaloHit->GetPositionVector().GetX())
                pSecondCaloHit = pCaloHit;
        }
    }

    if (NULL == pFirstCaloHit || NULL == pSecondCaloHit)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    candidateList.insert(pFirstCaloHit);
    candidateList.insert(pSecondCaloHit);

    // Find the two most separated hits
    float maxDistanceSquared(0.f);

    for (CaloHitList::const_iterator iterI = candidateList.begin(), iterEndI = candidateList.end(); iterI != iterEndI; ++iterI )
    {
        CaloHit *pCaloHitI = *iterI;

        for (CaloHitList::const_iterator iterJ = iterI, iterEndJ = candidateList.end(); iterJ != iterEndJ; ++iterJ )
        {
            CaloHit *pCaloHitJ = *iterJ;

            const float distanceSquared((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared());

            if (distanceSquared > maxDistanceSquared)
            {
                maxDistanceSquared = distanceSquared;
                pFirstCaloHit = pCaloHitI;
                pSecondCaloHit = pCaloHitJ;
            }
        }
    }

    // Set the inner and outer coordinates (Check Z first, then X in the event of a tie)
    const float deltaZ(pSecondCaloHit->GetPositionVector().GetZ() - pFirstCaloHit->GetPositionVector().GetZ());
    const float deltaX(pSecondCaloHit->GetPositionVector().GetX() - pFirstCaloHit->GetPositionVector().GetX());

    if ((deltaZ > 0.f) || ((std::fabs(deltaZ) < std::numeric_limits<float>::epsilon()) && (deltaX > 0.f)))
    {
        innerCoordinate = pFirstCaloHit->GetPositionVector();
        outerCoordinate = pSecondCaloHit->GetPositionVector();
    }
    else
    {
        innerCoordinate = pSecondCaloHit->GetPositionVector();
        outerCoordinate = pFirstCaloHit->GetPositionVector();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByInnerLayer(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int innerLayerLhs(pLhs->GetInnerPseudoLayer());
    const unsigned int innerLayerRhs(pRhs->GetInnerPseudoLayer());

    if( innerLayerLhs != innerLayerRhs )
      return (innerLayerLhs < innerLayerRhs);

    // Use SortByNOccupiedLayers method to resolve ties
    return SortByNOccupiedLayers(pLhs,pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNOccupiedLayers(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nOccupiedLayersLhs(pLhs->GetOrderedCaloHitList().size());
    const unsigned int nOccupiedLayersRhs(pRhs->GetOrderedCaloHitList().size());

    if (nOccupiedLayersLhs != nOccupiedLayersRhs)
    return (nOccupiedLayersLhs > nOccupiedLayersRhs);

    const unsigned int layerSpanLhs(pLhs->GetOuterPseudoLayer() - pLhs->GetInnerPseudoLayer());
    const unsigned int layerSpanRhs(pRhs->GetOuterPseudoLayer() - pRhs->GetInnerPseudoLayer());

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNHits(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nHitsLhs(pLhs->GetNCaloHits());
    const unsigned int nHitsRhs(pRhs->GetNCaloHits());

    if (nHitsLhs != nHitsRhs)
        return (nHitsLhs > nHitsRhs);

    const unsigned int layerSpanLhs(pLhs->GetOuterPseudoLayer() - pLhs->GetInnerPseudoLayer());
    const unsigned int layerSpanRhs(pRhs->GetOuterPseudoLayer() - pRhs->GetInnerPseudoLayer());

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::StoreSlidingFitResults(TwoDSlidingFitResult &twoDSlidingFitResult)
{
    const TwoDSlidingFitResult::LayerFitContributionMap &layerFitContributionMap(twoDSlidingFitResult.m_layerFitContributionMap);

    if (layerFitContributionMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    const int innerLayer(layerFitContributionMap.begin()->first);
    const int outerLayer(layerFitContributionMap.rbegin()->first);
    const int layerFitHalfWindow(twoDSlidingFitResult.m_layerFitHalfWindow);

    TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap(twoDSlidingFitResult.m_layerFitResultMap);

    if (!layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Use sliding fit to fill occupied layers
    unsigned int slidingNPoints(0);
    double slidingSumT(0.), slidingSumL(0.), slidingSumTT(0.), slidingSumLT(0.), slidingSumLL(0.);

    for (int iLayer = innerLayer; iLayer < static_cast<int>(innerLayer + layerFitHalfWindow); ++iLayer)
    {
        TwoDSlidingFitResult::LayerFitContributionMap::const_iterator lyrIter = layerFitContributionMap.find(iLayer);

        if (layerFitContributionMap.end() != lyrIter)
        {
            slidingSumT += lyrIter->second.GetSumT();
            slidingSumL += lyrIter->second.GetSumL();
            slidingSumTT += lyrIter->second.GetSumTT();
            slidingSumLT += lyrIter->second.GetSumLT();
            slidingSumLL += lyrIter->second.GetSumLL();
            slidingNPoints += lyrIter->second.GetNPoints();
        }
    }

    for (int iLayer = innerLayer; iLayer <= outerLayer; ++iLayer)
    {
        const int fwdLayer(iLayer + layerFitHalfWindow);
        TwoDSlidingFitResult::LayerFitContributionMap::const_iterator fwdIter = layerFitContributionMap.find(fwdLayer);

        if (layerFitContributionMap.end() != fwdIter)
        {
            slidingSumT += fwdIter->second.GetSumT();
            slidingSumL += fwdIter->second.GetSumL();
            slidingSumTT += fwdIter->second.GetSumTT();
            slidingSumLT += fwdIter->second.GetSumLT();
            slidingSumLL += fwdIter->second.GetSumLL();
            slidingNPoints += fwdIter->second.GetNPoints();
        }

        const int bwdLayer(iLayer - layerFitHalfWindow - 1);
        TwoDSlidingFitResult::LayerFitContributionMap::const_iterator bwdIter = layerFitContributionMap.find(bwdLayer);

        if (layerFitContributionMap.end() != bwdIter)
        {
            slidingSumT -= bwdIter->second.GetSumT();
            slidingSumL -= bwdIter->second.GetSumL();
            slidingSumTT -= bwdIter->second.GetSumTT();
            slidingSumLT -= bwdIter->second.GetSumLT();
            slidingSumLL -= bwdIter->second.GetSumLL();
            slidingNPoints -= bwdIter->second.GetNPoints();
        }

        // only fill the result map if there is an entry in the contribution map
        if (layerFitContributionMap.end() == layerFitContributionMap.find(iLayer))
            continue;

        // require three points for meaningful results
        if (slidingNPoints <=2)
            continue;

        const double denominator(slidingSumLL - slidingSumL * slidingSumL / static_cast<double>(slidingNPoints));

        if (std::fabs(denominator) < std::numeric_limits<double>::epsilon())
            continue;

        const double gradient((slidingSumLT - slidingSumL * slidingSumT / static_cast<double>(slidingNPoints)) / denominator);
        const double intercept((slidingSumLL * slidingSumT / static_cast<double>(slidingNPoints) - slidingSumL * slidingSumLT / static_cast<double>(slidingNPoints)) / denominator);

        const double l(twoDSlidingFitResult.GetL(iLayer));
        const double fitT(intercept + gradient * l);

        const double variance((slidingSumTT - 2. * intercept * slidingSumT - 2. * gradient * slidingSumLT + intercept * intercept * static_cast<double>(slidingNPoints) + 2. * gradient * intercept * slidingSumL + gradient * gradient * slidingSumLL) / (1. + gradient * gradient));

        if (variance < 0.)
            continue;

        const double rms(std::sqrt(variance / static_cast<double>(slidingNPoints)));

        const TwoDSlidingFitResult::TwoDSlidingFitResult::LayerFitResult layerFitResult(l, fitT, gradient, rms);
        (void) layerFitResultMap.insert(TwoDSlidingFitResult::LayerFitResultMap::value_type(iLayer, layerFitResult));
    }

    // TODO: Bail out if there are no results? This will break everything!
    // if (layerFitResultMap.empty())
    //     throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // Calculate sliding fit segments
    // TODO: Only calculate this information when required
    return LArClusterHelper::CalculateSlidingFitSegments(twoDSlidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::CalculateSlidingFitSegments(TwoDSlidingFitResult &twoDSlidingFitResult)
{
    unsigned int nSustainedSteps(0);
    float sustainedDirectionStartX(0.f), sustainedDirectionEndX(0.f);
    CartesianVector previousPosition(0.f, 0.f, 0.f);
    TransverseDirection previousDirection(UNKNOWN), sustainedDirection(UNKNOWN);
    TwoDSlidingFitResult::LayerFitResultMap::const_iterator sustainedDirectionStartIter, sustainedDirectionEndIter;

    const TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap(twoDSlidingFitResult.GetLayerFitResultMap());

    TwoDSlidingFitResult::FitSegmentList &fitSegmentList(twoDSlidingFitResult.m_fitSegmentList);

    if (!fitSegmentList.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        CartesianVector position(0.f, 0.f, 0.f);
        twoDSlidingFitResult.GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);

        const CartesianVector delta(position - previousPosition);
        const TransverseDirection currentDirection((delta.GetX() > 0.f) ? POSITIVE_IN_X : NEGATIVE_IN_X);
        // TODO: currentDirection could also be UNCHANGED_IN_X

        if (previousDirection == currentDirection)
        {
            ++nSustainedSteps;

            if (nSustainedSteps > 2)
            {
                sustainedDirection = currentDirection;
                sustainedDirectionEndIter = iter;
                sustainedDirectionEndX = position.GetX();
            }
        }
        else
        {
            if ((POSITIVE_IN_X == sustainedDirection) || (NEGATIVE_IN_X == sustainedDirection))
                fitSegmentList.push_back(TwoDSlidingFitResult::FitSegment(sustainedDirectionStartIter->first, sustainedDirectionEndIter->first,
                    sustainedDirectionStartX, sustainedDirectionEndX));

            nSustainedSteps = 0;
            sustainedDirection = UNKNOWN;
            sustainedDirectionStartIter = iter;
            sustainedDirectionStartX = position.GetX();
        }

        previousPosition = position;
        previousDirection = currentDirection;
    }

    if ((POSITIVE_IN_X == sustainedDirection) || (NEGATIVE_IN_X == sustainedDirection))
        fitSegmentList.push_back(TwoDSlidingFitResult::FitSegment(sustainedDirectionStartIter->first, sustainedDirectionEndIter->first,
            sustainedDirectionStartX, sustainedDirectionEndX));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArClusterHelper::m_layerFitHalfWindow = 20;
float LArClusterHelper::m_multiValuedTanThetaCut = 0.1f;
float LArClusterHelper::m_multiValuedStepFractionCut = 0.5f;
float LArClusterHelper::m_trackResidualQuantile = 0.8f;

StatusCode LArClusterHelper::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayerFitHalfWindow", m_layerFitHalfWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MultiValuedTanThetaCut", m_multiValuedTanThetaCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MultiValuedStepFractionCut", m_multiValuedStepFractionCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackResidualQuantile", m_trackResidualQuantile));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
