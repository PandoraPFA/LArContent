/**
 *  @file   larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.cc
 *
 *  @brief  Implementation of the lar two dimensional sliding shower fit result class.
 *
 *  $Log: $
 */

#include "Objects/Cluster.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

TwoDSlidingShowerFitResult::TwoDSlidingShowerFitResult(const Cluster *const pCluster, const unsigned int slidingFitWindow,
        const float slidingFitLayerPitch, const float showerEdgeMultiplier) :
    m_showerFitResult(TwoDSlidingFitResult(pCluster, slidingFitWindow, slidingFitLayerPitch)),
    m_negativeEdgeFitResult(TwoDSlidingShowerFitResult::LArTwoDShowerEdgeFit(m_showerFitResult, NEGATIVE_SHOWER_EDGE, showerEdgeMultiplier)),
    m_positiveEdgeFitResult(TwoDSlidingShowerFitResult::LArTwoDShowerEdgeFit(m_showerFitResult, POSITIVE_SHOWER_EDGE, showerEdgeMultiplier))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingShowerFitResult::GetShowerEdges(const float x, const bool widenIfAmbiguity, FloatVector &edgePositions) const
{
    edgePositions.clear();
    CartesianPointList fitPositionList;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->GetNegativeEdgeFitResult().GetGlobalFitPositionListAtX(x, fitPositionList));
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->GetPositiveEdgeFitResult().GetGlobalFitPositionListAtX(x, fitPositionList));

    if (fitPositionList.size() < 2)
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
        else if (fitPositionList.empty())
        {
            fitPositionList.push_back(CartesianVector(x, 0.f, minZ));
            fitPositionList.push_back(CartesianVector(x, 0.f, maxZ));
        }
        else if (1 == fitPositionList.size())
        {
            // ATTN Could improve sophistication of choice of second bounding edge
            const float existingEdge(fitPositionList.front().GetZ());
            const float secondEdge((std::fabs(existingEdge - minZ) < std::fabs(existingEdge - maxZ)) ? minZ : maxZ);
            fitPositionList.push_back(CartesianVector(x, 0.f, secondEdge));
        }
    }

    FloatVector localEdgePositions;
    for (CartesianPointList::const_iterator iter = fitPositionList.begin(), iterEnd = fitPositionList.end(); iter != iterEnd; ++iter)
        localEdgePositions.push_back(iter->GetZ());

    if (localEdgePositions.size() < 2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    std::sort(localEdgePositions.begin(), localEdgePositions.end());
    edgePositions.push_back(localEdgePositions.front());
    edgePositions.push_back(localEdgePositions.back());
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingFitResult TwoDSlidingShowerFitResult::LArTwoDShowerEdgeFit(const TwoDSlidingFitResult &fullShowerFit, const ShowerEdge showerEdge,
    const float showerEdgeMultiplier)
{
    LayerFitContributionMap layerFitContributionMap;

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
            fullShowerFit.GetLocalPosition((*hitIter)->GetPositionVector(), rL, rT);
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
    }

    // Select fit contributions representing relevant shower edge
    for (FitCoordinateMap::const_iterator iter = fitCoordinateMap.begin(), iterEnd = fitCoordinateMap.end(); iter != iterEnd; ++iter)
    {
        const int layer(iter->first);
        const FitCoordinateList &fitCoordinateList(iter->second);

        // ATTN Could modify this hit selection, e.g. add inertia to edge positions
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

    return TwoDSlidingFitResult(fullShowerFit.GetCluster(), fullShowerFit.GetLayerFitHalfWindow(), fullShowerFit.GetLayerPitch(),
        fullShowerFit.GetAxisIntercept(), fullShowerFit.GetAxisDirection(), fullShowerFit.GetOrthoDirection(), layerFitContributionMap);
}

} // namespace lar_content
