/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CrossedTrackSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the crossed track splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar_content
{

StatusCode CrossedTrackSplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2,
    CartesianVector &splitPosition, CartesianVector &firstDirection, CartesianVector &secondDirection) const
{
    // Identify crossed-track topology and find candidate intersection positions
    const CartesianVector& minPosition1(slidingFitResult1.GetGlobalMinLayerPosition());
    const CartesianVector& maxPosition1(slidingFitResult1.GetGlobalMaxLayerPosition());

    const CartesianVector& minPosition2(slidingFitResult2.GetGlobalMinLayerPosition());
    const CartesianVector& maxPosition2(slidingFitResult2.GetGlobalMaxLayerPosition());

    if (LArClusterHelper::GetClosestDistance(minPosition1, slidingFitResult2.GetCluster()) < 2.f * m_maxClusterSeparation ||
        LArClusterHelper::GetClosestDistance(maxPosition1, slidingFitResult2.GetCluster()) < 2.f * m_maxClusterSeparation ||
        LArClusterHelper::GetClosestDistance(minPosition2, slidingFitResult1.GetCluster()) < 2.f * m_maxClusterSeparation ||
        LArClusterHelper::GetClosestDistance(maxPosition2, slidingFitResult1.GetCluster()) < 2.f * m_maxClusterSeparation)
        return STATUS_CODE_NOT_FOUND;

    CartesianPointList candidateList;
    this->FindCandidateSplitPositions(slidingFitResult1.GetCluster(), slidingFitResult2.GetCluster(), candidateList);

    if (candidateList.empty())
        return STATUS_CODE_NOT_FOUND;


    // Loop over candidate positions and find best split position
    bool foundSplit(false);
    float closestSeparationSquared(std::numeric_limits<float>::max());

    const float halfWindowLength1(slidingFitResult1.GetLayerFitHalfWindowLength());
    const float halfWindowLength2(slidingFitResult2.GetLayerFitHalfWindowLength());

    for (CartesianPointList::const_iterator iter = candidateList.begin(), iterEnd = candidateList.end(); iter != iterEnd; ++iter)
    {
        const CartesianVector &candidatePosition(*iter);

        try
        {
            // Projections onto first cluster
            float rL1(0.f), rT1(0.f);
            CartesianVector R1(0.f, 0.f, 0.f);
            CartesianVector F1(0.f, 0.f, 0.f);
            CartesianVector B1(0.f, 0.f, 0.f);

            slidingFitResult1.GetGlobalFitProjection(candidatePosition, R1);
            slidingFitResult1.GetLocalPosition(R1, rL1, rT1);
            slidingFitResult1.GetGlobalFitPosition(rL1 + halfWindowLength1, F1);
            slidingFitResult1.GetGlobalFitPosition(rL1 - halfWindowLength1, B1);

            // Projections onto second cluster
            float rL2(0.f), rT2(0.f);
            CartesianVector R2(0.f, 0.f, 0.f);
            CartesianVector F2(0.f, 0.f, 0.f);
            CartesianVector B2(0.f, 0.f, 0.f);

            slidingFitResult2.GetGlobalFitProjection(candidatePosition, R2);
            slidingFitResult2.GetLocalPosition(R2, rL2, rT2);
            slidingFitResult2.GetGlobalFitPosition(rL2 + halfWindowLength2, F2);
            slidingFitResult2.GetGlobalFitPosition(rL2 - halfWindowLength2, B2);

            // Calculate average position
            const CartesianVector C0((R1 + R2) * 0.5);

            // Calculate intersected position:
            // ==============================
            // First cluster gives set of points: B1->R1->F1
            // Second cluster gives set of points: B2->R2->F2
            //
            // Try swapping B1 with B2 to see if this gives intersecting straight lines:
            //
            //   F1   F2     a2   b1
            //    |   |       |   |
            //     |  |        |  |
            //     R1 R2       R1 R2
            //      | |         |  |
            //      |  |        |   |
            //     B1   B2     a1    b2

            // First straight line is a1->R1->b1
            // Second straight line is a2->R2->b2

            const CartesianVector a1(B1);
            const CartesianVector a2(F1);

            for (unsigned int iForward = 0; iForward<2; ++iForward)
            {
                const CartesianVector b1((0 == iForward) ? F2 : B2);
                const CartesianVector b2((0 == iForward) ? B2 : F2);

                const CartesianVector s1((b1 - R2).GetUnitVector());
                const CartesianVector t1((R1 - a1).GetUnitVector());
                const CartesianVector s2((b2 - R2).GetUnitVector());
                const CartesianVector t2((R1 - a2).GetUnitVector());

                if (s1.GetDotProduct(t1) < std::max(m_minCosRelativeAngle,-s1.GetDotProduct(s2)) ||
                    s2.GetDotProduct(t2) < std::max(m_minCosRelativeAngle,-t1.GetDotProduct(t2)))
                    continue;

                const CartesianVector p1((b1 - a1).GetUnitVector());
                const CartesianVector p2((b2 - a2).GetUnitVector());

                float mu1(0.f), mu2(0.f);
                CartesianVector C1(0.f,0.f,0.f);

                LArPointingClusterHelper::GetIntersection(a1, p1, a2, p2, C1, mu1, mu2);

                if (mu1 < 0.f || mu2 < 0.f || mu1 > (b1 - a1).GetMagnitude() || mu2 > (b2 - a2).GetMagnitude())
                    continue;

                const float thisSeparationSquared((C0 - C1).GetMagnitudeSquared());

                if (thisSeparationSquared < closestSeparationSquared)
                {
                    closestSeparationSquared = thisSeparationSquared;
                    splitPosition = (C0 + C1) * 0.5;
                    firstDirection = t2 * -1.f;
                    secondDirection = t1;
                    foundSplit = true;
                }
            }
        }
        catch (StatusCodeException &)
        {

        }
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossedTrackSplittingAlgorithm::FindCandidateSplitPositions(const Cluster *const pCluster1, const Cluster *const pCluster2,
    CartesianPointList &candidateList) const
{
    // ATTN The following is double-double counting
    CaloHitList caloHitList1, caloHitList2;
    pCluster1->GetOrderedCaloHitList().GetCaloHitList(caloHitList1);
    pCluster2->GetOrderedCaloHitList().GetCaloHitList(caloHitList2);

    for (CaloHitList::const_iterator iter1 = caloHitList1.begin(), iterEnd1 = caloHitList1.end(); iter1 != iterEnd1; ++iter1)
    {
        CaloHit *pCaloHit = *iter1;

        const CartesianVector position1(pCaloHit->GetPositionVector());
        const CartesianVector position2(LArClusterHelper::GetClosestPosition(position1, pCluster2));

        if ((position1 - position2).GetMagnitudeSquared() < m_maxClusterSeparationSquared)
          candidateList.push_back((position1 + position2) * 0.5);
    }

    for (CaloHitList::const_iterator iter2 = caloHitList2.begin(), iterEnd2 = caloHitList2.end(); iter2 != iterEnd2; ++iter2)
    {
        CaloHit *pCaloHit = *iter2;

        const CartesianVector position2(pCaloHit->GetPositionVector());
        const CartesianVector position1(LArClusterHelper::GetClosestPosition(position2, pCluster1));

        if ((position2 - position1).GetMagnitudeSquared() < m_maxClusterSeparationSquared)
          candidateList.push_back((position2 + position1) * 0.5);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CrossedTrackSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_maxClusterSeparation = 2.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterSeparation", m_maxClusterSeparation));
    m_maxClusterSeparationSquared = m_maxClusterSeparation * m_maxClusterSeparation;

    m_minCosRelativeAngle = 0.966;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    return TwoDSlidingFitSplittingAndSwitchingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
