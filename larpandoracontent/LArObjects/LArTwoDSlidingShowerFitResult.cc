/**
 *  @file   larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.cc
 *
 *  @brief  Implementation of the lar two dimensional sliding shower fit result class.
 *
 *  $Log: $
 */

#include "Objects/Cluster.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

template <typename T>
TwoDSlidingShowerFitResult::TwoDSlidingShowerFitResult(const T *const pT, const unsigned int slidingFitWindow, const float slidingFitLayerPitch,
        const float showerEdgeMultiplier) :
    m_showerFitResult(TwoDSlidingFitResult(pT, slidingFitWindow, slidingFitLayerPitch)),
    m_negativeEdgeFitResult(TwoDSlidingShowerFitResult::LArTwoDShowerEdgeFit(pT, m_showerFitResult, NEGATIVE_SHOWER_EDGE, showerEdgeMultiplier)),
    m_positiveEdgeFitResult(TwoDSlidingShowerFitResult::LArTwoDShowerEdgeFit(pT, m_showerFitResult, POSITIVE_SHOWER_EDGE, showerEdgeMultiplier))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingShowerFitResult::GetShowerEdges(const float x, const bool widenIfAmbiguity, FloatVector &edgePositions) const
{
    edgePositions.clear();
    CartesianPointVector fitPositionVector;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->GetNegativeEdgeFitResult().GetGlobalFitPositionListAtX(x, fitPositionVector));
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->GetPositiveEdgeFitResult().GetGlobalFitPositionListAtX(x, fitPositionVector));

    if (fitPositionVector.size() < 2)
    {
        float minXn(0.f), maxXn(0.f), minXp(0.f), maxXp(0.f);
        this->GetNegativeEdgeFitResult().GetMinAndMaxX(minXn, maxXn);
        this->GetPositiveEdgeFitResult().GetMinAndMaxX(minXp, maxXp);
        const float minX(std::min(minXn, minXp)), maxX(std::max(maxXn, maxXp));

        if ((x < minX) || (x > maxX))
            return;

        float minZn(0.f), maxZn(0.f), minZp(0.f), maxZp(0.f);
        this->GetNegativeEdgeFitResult().GetMinAndMaxZ(minZn, maxZn);
        this->GetPositiveEdgeFitResult().GetMinAndMaxZ(minZp, maxZp);
        const float minZ(std::min(minZn, minZp)), maxZ(std::max(maxZn, maxZp));

        if (!widenIfAmbiguity)
        {
            return;
        }
        else if (fitPositionVector.empty())
        {
            fitPositionVector.push_back(CartesianVector(x, 0.f, minZ));
            fitPositionVector.push_back(CartesianVector(x, 0.f, maxZ));
        }
        else if (1 == fitPositionVector.size())
        {
            // ATTN Could improve sophistication of choice of second bounding edge
            const float existingEdge(fitPositionVector.front().GetZ());
            const float secondEdge((std::fabs(existingEdge - minZ) < std::fabs(existingEdge - maxZ)) ? minZ : maxZ);
            fitPositionVector.push_back(CartesianVector(x, 0.f, secondEdge));
        }
    }

    FloatVector localEdgePositions;
    for (const CartesianVector &fitPosition : fitPositionVector)
        localEdgePositions.push_back(fitPosition.GetZ());

    if (localEdgePositions.size() < 2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    std::sort(localEdgePositions.begin(), localEdgePositions.end());
    edgePositions.push_back(localEdgePositions.front());
    edgePositions.push_back(localEdgePositions.back());
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingFitResult TwoDSlidingShowerFitResult::LArTwoDShowerEdgeFit(const Cluster *const pCluster, const TwoDSlidingFitResult &fullShowerFit,
    const ShowerEdge showerEdge, const float showerEdgeMultiplier)
{
    CartesianPointVector pointVector;
    LArClusterHelper::GetCoordinateVector(pCluster, pointVector);
    return TwoDSlidingShowerFitResult::LArTwoDShowerEdgeFit(&pointVector, fullShowerFit, showerEdge, showerEdgeMultiplier);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingFitResult TwoDSlidingShowerFitResult::LArTwoDShowerEdgeFit(const CartesianPointVector *const pPointVector,
    const TwoDSlidingFitResult &fullShowerFit, const ShowerEdge showerEdge, const float showerEdgeMultiplier)
{
    // Examine all possible fit contributions
    FitCoordinateMap fitCoordinateMap;

    for (const CartesianVector &hitPosition : *pPointVector)
    {
        float rL(0.f), rT(0.f);
        fullShowerFit.GetLocalPosition(hitPosition, rL, rT);
        rT *= showerEdgeMultiplier;

        CartesianVector fullShowerFitPosition(0.f, 0.f, 0.f);
        if (STATUS_CODE_SUCCESS != fullShowerFit.GetGlobalFitPosition(rL, fullShowerFitPosition))
            continue;

        float rLFit(0.f), rTFit(0.f);
        fullShowerFit.GetLocalPosition(fullShowerFitPosition, rLFit, rTFit);

        const float rTDiff(rT - rTFit);
        if (((POSITIVE_SHOWER_EDGE == showerEdge) && (rTDiff < 0.f)) || ((NEGATIVE_SHOWER_EDGE == showerEdge) && (rTDiff > 0.f)))
            rT = rTFit;

        const int layer(fullShowerFit.GetLayer(rL));
        fitCoordinateMap[layer].push_back(FitCoordinate(rL, rT));
    }

    // Select fit contributions representing relevant shower edge
    LayerFitContributionMap layerFitContributionMap;

    for (const FitCoordinateMap::value_type &mapEntry : fitCoordinateMap)
    {
        // ATTN Could modify this hit selection, e.g. add inertia to edge positions
        bool bestFitCoordinateFound(false);
        FitCoordinate bestFitCoordinate = (POSITIVE_SHOWER_EDGE == showerEdge) ?
            FitCoordinate(0.f, -std::numeric_limits<float>::max()) :
            FitCoordinate(0.f, +std::numeric_limits<float>::max());

        for (const FitCoordinate &fitCoordinate : mapEntry.second)
        {
            if (((POSITIVE_SHOWER_EDGE == showerEdge) && (fitCoordinate.second > bestFitCoordinate.second)) ||
                ((NEGATIVE_SHOWER_EDGE == showerEdge) && (fitCoordinate.second < bestFitCoordinate.second)))
            {
                bestFitCoordinate = fitCoordinate;
                bestFitCoordinateFound = true;
            }
        }

        if (bestFitCoordinateFound)
            layerFitContributionMap[mapEntry.first].AddPoint(bestFitCoordinate.first, bestFitCoordinate.second);
    }

    return TwoDSlidingFitResult(fullShowerFit.GetLayerFitHalfWindow(), fullShowerFit.GetLayerPitch(), fullShowerFit.GetAxisIntercept(),
        fullShowerFit.GetAxisDirection(), fullShowerFit.GetOrthoDirection(), layerFitContributionMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template TwoDSlidingShowerFitResult::TwoDSlidingShowerFitResult(const pandora::Cluster *const, const unsigned int, const float, const float);
template TwoDSlidingShowerFitResult::TwoDSlidingShowerFitResult(const pandora::CartesianPointVector *const, const unsigned int, const float, const float);

} // namespace lar_content
