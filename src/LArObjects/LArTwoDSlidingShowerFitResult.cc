/**
 *  @file   LArContent/src/LArObjects/LArTwoDSlidingShowerFitResult.cc
 *
 *  @brief  Implementation of the lar two dimensional sliding shower fit result class.
 *
 *  $Log: $
 */

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArTwoDSlidingShowerFitResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar
{

TwoDSlidingShowerFitResult::TwoDSlidingShowerFitResult(const Cluster *const pCluster, const unsigned int slidingFitWindow)
{
    LArClusterHelper::LArTwoDSlidingFit(pCluster, slidingFitWindow, m_showerFitResult);
    LArClusterHelper::LArTwoDShowerEdgeFit(m_showerFitResult, NEGATIVE_SHOWER_EDGE, m_negativeEdgeFitResult);
    LArClusterHelper::LArTwoDShowerEdgeFit(m_showerFitResult, POSITIVE_SHOWER_EDGE, m_positiveEdgeFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingShowerFitResult::GetShowerEdges(const float x, FloatVector &edgePositions) const
{
    edgePositions.clear();
    CartesianPointList fitPositionList;
    try {this->GetNegativeEdgeFitResult().GetGlobalFitPositionListAtX(x, fitPositionList);} catch (StatusCodeException &) {}
    try {this->GetPositiveEdgeFitResult().GetGlobalFitPositionListAtX(x, fitPositionList);} catch (StatusCodeException &) {}

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

        if (fitPositionList.empty())
        {
            fitPositionList.push_back(CartesianVector(x, 0.f, minZ));
            fitPositionList.push_back(CartesianVector(x, 0.f, maxZ));
        }
        else if (1 == fitPositionList.size())
        {
            // TODO More sophisticated choice of second bounding edge
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

} // namespace lar
