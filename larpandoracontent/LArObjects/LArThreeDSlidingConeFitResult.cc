/**
 *  @file   larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.cc
 *
 *  @brief  Implementation of the lar three dimensional sliding cone fit result class.
 *
 *  $Log: $
 */

#include "Objects/Cluster.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"

#include <iterator>

using namespace pandora;

namespace lar_content
{

float SimpleCone::GetBoundedHitFraction(const Cluster *const pCluster) const
{
    unsigned int nMatchedHits(0);
    const unsigned int nClusterHits(pCluster->GetNCaloHits());

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (const auto &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *pCaloHit : *(mapEntry.second))
        {
            const CartesianVector displacement(pCaloHit->GetPositionVector() - this->GetConeApex());
            const float rL(displacement.GetDotProduct(this->GetConeDirection()));

            if ((rL < 0.f) || (rL > this->GetConeLength()))
                continue;

            const float rT(displacement.GetCrossProduct(this->GetConeDirection()).GetMagnitude());
            
            if (rL * this->GetConeTanHalfAngle() > rT)
                ++nMatchedHits;
        }
    }

    const float boundedFraction((nClusterHits > 0) ? static_cast<float>(nMatchedHits) / static_cast<float>(nClusterHits) : 0.f);
    return boundedFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDSlidingConeFitResult::ThreeDSlidingConeFitResult(const Cluster *const pCluster, const unsigned int slidingFitWindow,
        const float slidingFitLayerPitch) :
    m_slidingFitResult(ThreeDSlidingFitResult(pCluster, slidingFitWindow, slidingFitLayerPitch))
{
    const CartesianVector &minLayerPosition3D(m_slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector &maxLayerPosition3D(m_slidingFitResult.GetGlobalMaxLayerPosition());

    const unsigned int nSteps(static_cast<unsigned int>((maxLayerPosition3D - minLayerPosition3D).GetMagnitude() / slidingFitLayerPitch));

    for (unsigned int iStep = 0; iStep < nSteps; ++iStep)
    {
        try
        {
            const float rL((static_cast<float>(iStep) + 0.5f) * slidingFitLayerPitch);
            CartesianVector fitPosition3D(0.f, 0.f, 0.f), fitDirection3D(0.f, 0.f, 0.f);

            if ((STATUS_CODE_SUCCESS != m_slidingFitResult.GetGlobalFitPosition(rL, fitPosition3D)) ||
                (STATUS_CODE_SUCCESS != m_slidingFitResult.GetGlobalFitDirection(rL, fitDirection3D)))
            {
                continue;
            }

            (void) m_trackStateMap.insert(TrackStateMap::value_type(iStep, TrackState(fitPosition3D, fitDirection3D)));
        }
        catch (const StatusCodeException &)
        {
            /* Deliberately empty */
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDSlidingConeFitResult::GetSimpleConeList(const unsigned int nLayersForConeFit, const unsigned int nCones, SimpleConeList &simpleConeList) const
{
    const TrackStateMap &trackStateMap(this->GetTrackStateMap());
    const unsigned int nLayers(trackStateMap.size());

    if (nLayers < nLayersForConeFit)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    //const unsigned int segmentLayers();
    const unsigned int coneInterval((nCones <= 2) ? 1 : nLayers / (nCones - 2)); // TODO
    unsigned int nConeSamplingSteps(0);

    TrackStateLinkedList trackStateList;
    CartesianVector directionSum(0.f, 0.f, 0.f);
    const float clusterLength((trackStateMap.begin()->second.GetPosition() - trackStateMap.rbegin()->second.GetPosition()).GetMagnitude());

    for (TrackStateMap::const_iterator iter = trackStateMap.begin(), iterEnd = trackStateMap.end(); iter != iterEnd; ++iter)
    {
        if (nConeSamplingSteps >= nCones)
            return;

        trackStateList.push_back(iter->second);
        directionSum += iter->second.GetMomentum();

        const unsigned int beginDistance(static_cast<unsigned int>(std::distance(trackStateMap.begin(), iter)));
        const unsigned int endDistance(static_cast<unsigned int>(std::distance(iter, trackStateMap.end())));

        if ((beginDistance < nLayersForConeFit) || (endDistance < nLayersForConeFit))
            continue;

        const TrackState &maxLayerTrackState(trackStateList.back());
        const TrackState &minLayerTrackState(trackStateList.front());

        if ((beginDistance == nLayersForConeFit) || (endDistance == nLayersForConeFit) || (beginDistance % coneInterval == 0))
        {
            const CartesianVector &minLayerApex(minLayerTrackState.GetPosition());
            const CartesianVector &maxLayerApex(maxLayerTrackState.GetPosition());

            const CartesianVector minLayerDirection(directionSum.GetUnitVector());
            const CartesianVector maxLayerDirection(directionSum.GetUnitVector() * -1.f);

            // TODO Estimate cone length and angle here too
            ++nConeSamplingSteps;
            simpleConeList.push_back(SimpleCone(minLayerApex, minLayerDirection, clusterLength * 3.f, 0.5f)); // TODO
            simpleConeList.push_back(SimpleCone(maxLayerApex, maxLayerDirection, clusterLength * 3.f, 0.5f));
         }

        directionSum -= minLayerTrackState.GetMomentum();
        trackStateList.pop_front();
    }
}

} // namespace lar_content
