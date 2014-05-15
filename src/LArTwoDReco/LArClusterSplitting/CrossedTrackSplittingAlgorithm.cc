/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CrossedTrackSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar
{

StatusCode CrossedTrackSplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2,
    CartesianVector &splitPosition, CartesianVector &firstDirection, CartesianVector &secondDirection) const
{
    //
    // TODO: SPEED UP!
    //

    const float m_maxSeparationSquared = 2.5f * 2.5f;

    bool foundSplit(false);
    float closestSeparationSquared(m_maxSeparationSquared);

    const float halfWindowLength1(slidingFitResult1.GetLayerFitHalfWindowLength());
    const float halfWindowLength2(slidingFitResult2.GetLayerFitHalfWindowLength());

    const TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap1(slidingFitResult1.GetLayerFitResultMap());
    const TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap2(slidingFitResult2.GetLayerFitResultMap());

    for (TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter1 = layerFitResultMap1.begin(), iterEnd1 = layerFitResultMap1.end();
        iter1 != iterEnd1; ++iter1)
    {
        const float rL1(iter1->second.GetL());
        const float rT1(iter1->second.GetFitT());

        CartesianVector R1(0.f, 0.f, 0.f);
        CartesianVector F1(0.f, 0.f, 0.f);
        CartesianVector B1(0.f, 0.f, 0.f);

        try
        {
            slidingFitResult1.GetGlobalPosition(rL1, rT1, R1);
            slidingFitResult1.GetGlobalFitPosition(rL1 + halfWindowLength1, F1);
            slidingFitResult1.GetGlobalFitPosition(rL1 - halfWindowLength1, B1);
        }
        catch (StatusCodeException &)
        {
            continue;
        }

        for (TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter2 = layerFitResultMap2.begin(), iterEnd2 = layerFitResultMap2.end();
            iter2 != iterEnd2; ++iter2)
        {
            const float rL2(iter2->second.GetL());
            const float rT2(iter2->second.GetFitT());

            CartesianVector R2(0.f, 0.f, 0.f);
            CartesianVector F2(0.f, 0.f, 0.f);
            CartesianVector B2(0.f, 0.f, 0.f);

            try
            {
                slidingFitResult2.GetGlobalPosition(rL2, rT2, R2);
                slidingFitResult2.GetGlobalFitPosition(rL2 + halfWindowLength2, F2);
                slidingFitResult2.GetGlobalFitPosition(rL2 - halfWindowLength2, B2);
            }
            catch (StatusCodeException &)
            {
                continue;
            }

            if ((R1 - R2).GetMagnitudeSquared() > m_maxSeparationSquared)
                continue;

            const CartesianVector C0((R1 + R2) * 0.5);

            const CartesianVector a1(B1);
            const CartesianVector a2(F1);

            for (unsigned int iForward = 0; iForward<2; ++iForward)
            {
                CartesianVector b1((0 == iForward) ? F2 : B2);
                CartesianVector b2((0 == iForward) ? B2 : F2);

                if ((b1 - C0).GetDotProduct(a1 - C0) > 0 || (b2 - C0).GetDotProduct(a2 - C0) > 0)
                    continue;

                //
                // ADD MORE SELECTION HERE
                //

                try{
                    float mu1(0.f), mu2(0.f);
                    CartesianVector C1(0.f,0.f,0.f);

                    const CartesianVector p1((b1 - a1).GetUnitVector());
                    const CartesianVector p2((b2 - a2).GetUnitVector());
                    LArPointingClusterHelper::GetIntersection(a1, p1, a2, p2, C1, mu1, mu2);

                    const float thisSeparationSquared((C0 - C1).GetMagnitudeSquared());

                    if (thisSeparationSquared < closestSeparationSquared)
                    {
                        closestSeparationSquared = thisSeparationSquared;
                        splitPosition = (C0 + C1) * 0.5;
                        firstDirection = (b1 - C0).GetUnitVector();
                        secondDirection = (C0 - a1).GetUnitVector();
                        foundSplit = true;
                    }
                }
                catch (StatusCodeException &)
                {

                }

            }
        }
    }

    if (foundSplit)
        return STATUS_CODE_SUCCESS;

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CrossedTrackSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{


    return TwoDSlidingFitSplittingAndSwitchingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
