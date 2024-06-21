/**
 *  @file   larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.cc
 *
 *  @brief  Implementation of the lar three dimensional sliding cone fit result class.
 *
 *  $Log: $
 */

#include "Objects/Cluster.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"

#include <iterator>

using namespace pandora;

namespace lar_content
{

float SimpleCone::GetMeanRT(const Cluster *const pCluster) const
{
    CartesianPointVector hitPositionVector;
    LArClusterHelper::GetCoordinateVector(pCluster, hitPositionVector);

    float rTSum(0.f);
    const unsigned int nClusterHits(pCluster->GetNCaloHits());

    for (const CartesianVector &hitPosition : hitPositionVector)
    {
        const CartesianVector displacement(hitPosition - this->GetConeApex());
        const float rT(displacement.GetCrossProduct(this->GetConeDirection()).GetMagnitude());
        rTSum += rT;
    }

    return ((nClusterHits > 0) ? rTSum / static_cast<float>(nClusterHits) : 0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float SimpleCone::GetBoundedHitFraction(const Cluster *const pCluster, const float coneLength, const float coneTanHalfAngle) const
{
    CartesianPointVector hitPositionVector;
    LArClusterHelper::GetCoordinateVector(pCluster, hitPositionVector);

    unsigned int nMatchedHits(0);
    const unsigned int nClusterHits(pCluster->GetNCaloHits());

    for (const CartesianVector &hitPosition : hitPositionVector)
    {
        const CartesianVector displacement(hitPosition - this->GetConeApex());
        const float rL(displacement.GetDotProduct(this->GetConeDirection()));

        if ((rL < 0.f) || (rL > coneLength))
            continue;

        const float rT(displacement.GetCrossProduct(this->GetConeDirection()).GetMagnitude());

        if (rL * coneTanHalfAngle > rT)
            ++nMatchedHits;
    }

    return ((nClusterHits > 0) ? static_cast<float>(nMatchedHits) / static_cast<float>(nClusterHits) : 0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
ThreeDSlidingConeFitResult::ThreeDSlidingConeFitResult(const T *const pT, const unsigned int slidingFitWindow, const float slidingFitLayerPitch) :
    m_slidingFitResult(ThreeDSlidingFitResult(pT, slidingFitWindow, slidingFitLayerPitch))
{
    const CartesianVector &minLayerPosition3D(m_slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector &maxLayerPosition3D(m_slidingFitResult.GetGlobalMaxLayerPosition());

    const TwoDSlidingFitResult &fitResult1(m_slidingFitResult.GetFirstFitResult());
    const TwoDSlidingFitResult &fitResult2(m_slidingFitResult.GetSecondFitResult());

    const LayerFitContributionMap &contributionMap1(fitResult1.GetLayerFitContributionMap());
    const LayerFitContributionMap &contributionMap2(fitResult2.GetLayerFitContributionMap());

    const int nSteps(static_cast<int>((maxLayerPosition3D - minLayerPosition3D).GetMagnitude() / slidingFitLayerPitch));

    for (int iStep = 0; iStep <= nSteps; ++iStep)
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

            // ATTN Do not include layers without contributions in track state map (contain only interpolations)
            if (!contributionMap1.count(fitResult1.GetLayer(rL)) && !contributionMap2.count(fitResult2.GetLayer(rL)))
                continue;

            (void)m_trackStateMap.insert(TrackStateMap::value_type(iStep, TrackState(fitPosition3D, fitDirection3D)));
        }
        catch (const StatusCodeException &)
        {
            /* Deliberately empty */
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDSlidingConeFitResult::GetSimpleConeList(
    const unsigned int nLayersForConeFit, const unsigned int nCones, const ConeSelection coneSelection, SimpleConeList &simpleConeList,
    const float tanHalfAngle, const bool legacyMode) const
{
    const TrackStateMap &trackStateMap(this->GetTrackStateMap());
    const unsigned int nLayers(trackStateMap.size());

    if (nLayers + 1 < nLayersForConeFit + nCones)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Calculate intervals and offsets such that cones are evenly spaced along sliding fit and are equidistant from the extremal layers
    const unsigned int coneInterval((nCones > 1) ? (nLayers - nLayersForConeFit) / (nCones - 1) : 1);
    const ThreeDSlidingFitResult &slidingFitResult(this->GetSlidingFitResult());
    const CartesianVector coneDisplacement(slidingFitResult.GetGlobalMaxLayerPosition() - slidingFitResult.GetGlobalMinLayerPosition());
    const bool isForward(coneDisplacement.GetZ() > std::numeric_limits<float>::epsilon());
    const unsigned int coneOffset1((nLayers - nLayersForConeFit - (nCones - 1) * coneInterval) / 2);
    const unsigned int coneOffset2((1 == nLayers % 2) && isForward ? 1 : 0);
    const unsigned int coneOffset(coneOffset1 + coneOffset2);

    TrackStateLinkedList trackStateList;
    CartesianVector directionSum(0.f, 0.f, 0.f);
    const float clusterLength((trackStateMap.begin()->second.GetPosition() - trackStateMap.rbegin()->second.GetPosition()).GetMagnitude());
    const float lengthStep{(nCones > 1) ? 0.5f * clusterLength / (nCones - 1) : 0.5f * clusterLength};
    const float angleStep{(nCones > 1) ? 0.5f * tanHalfAngle / (nCones - 1) : 0.5f * tanHalfAngle};

    unsigned int nConeSamplingSteps(0);

    for (TrackStateMap::const_iterator iter = trackStateMap.begin(), iterEnd = trackStateMap.end(); iter != iterEnd; ++iter)
    {
        if (nConeSamplingSteps >= nCones)
            return;

        trackStateList.push_back(iter->second);
        directionSum += iter->second.GetMomentum();

        const unsigned int beginDistance(static_cast<unsigned int>(std::distance(trackStateMap.begin(), iter)));

        if (beginDistance + 1 < nLayersForConeFit)
            continue;

        const TrackState &maxLayerTrackState(trackStateList.back());
        const TrackState &minLayerTrackState(trackStateList.front());

        if ((beginDistance + 1 >= nLayersForConeFit + coneOffset) && (beginDistance + 1 - nLayersForConeFit - coneOffset) % coneInterval == 0)
        {
            const CartesianVector &minLayerApex(minLayerTrackState.GetPosition());
            const CartesianVector &maxLayerApex(maxLayerTrackState.GetPosition());

            const CartesianVector minLayerDirection(directionSum.GetUnitVector());
            const CartesianVector maxLayerDirection(directionSum.GetUnitVector() * -1.f);

            // TODO Estimate cone length and angle here too, maybe by projecting positions onto direction and looking at rT distribution?
            ++nConeSamplingSteps;

            if ((CONE_FORWARD_ONLY == coneSelection) || (CONE_BOTH_DIRECTIONS == coneSelection))
            {
                if (legacyMode)
                {
                    simpleConeList.push_back(SimpleCone(minLayerApex, minLayerDirection, clusterLength, tanHalfAngle));
                }
                else
                {
                    const float stepConeLength{clusterLength - (nConeSamplingSteps - 1) * lengthStep};
                    const float stepTanHalfAngle{tanHalfAngle - (nCones - nConeSamplingSteps) * angleStep};
                    simpleConeList.push_back(SimpleCone(minLayerApex, minLayerDirection, stepConeLength, stepTanHalfAngle));
                }
            }

            if ((CONE_BACKWARD_ONLY == coneSelection) || (CONE_BOTH_DIRECTIONS == coneSelection))
            {
                if (legacyMode)
                {
                    simpleConeList.push_back(SimpleCone(maxLayerApex, maxLayerDirection, clusterLength, tanHalfAngle));
                }
                else
                {
                    const float stepConeLength{clusterLength - (nCones - nConeSamplingSteps) * lengthStep};
                    const float stepTanHalfAngle{tanHalfAngle - (nConeSamplingSteps - 1) * angleStep};
                    simpleConeList.push_back(SimpleCone(maxLayerApex, maxLayerDirection, stepConeLength, stepTanHalfAngle));
                }
            }
        }

        directionSum -= minLayerTrackState.GetMomentum();
        trackStateList.pop_front();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template ThreeDSlidingConeFitResult::ThreeDSlidingConeFitResult(const pandora::Cluster *const, const unsigned int, const float);
template ThreeDSlidingConeFitResult::ThreeDSlidingConeFitResult(const pandora::CartesianPointVector *const, const unsigned int, const float);

} // namespace lar_content
