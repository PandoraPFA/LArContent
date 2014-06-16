/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the track consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"

using namespace pandora;

namespace lar
{

void TrackConsolidationAlgorithm::GetReclusteredHits(const TwoDSlidingFitResultList &slidingFitResultListI,
    const ClusterVector &showerClustersJ, ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
    for (TwoDSlidingFitResultList::const_iterator iterI = slidingFitResultListI.begin(), iterEndI = slidingFitResultListI.end(); iterI != iterEndI; ++iterI)
    {
        const TwoDSlidingFitResult &slidingFitResultI = *iterI;

        const Cluster* pClusterI = slidingFitResultI.GetCluster();
        const float thisLengthSquaredI(LArClusterHelper::GetLengthSquared(pClusterI));

        for (ClusterVector::const_iterator iterJ = showerClustersJ.begin(), iterEndJ = showerClustersJ.end(); iterJ != iterEndJ; ++iterJ)
        {
            const Cluster* pClusterJ = *iterJ;
            const float thisLengthSquaredJ(LArClusterHelper::GetLengthSquared(pClusterJ));

            if (pClusterI == pClusterJ)
                continue;

            if (2.0 * thisLengthSquaredJ > thisLengthSquaredI)
                continue;

            this->GetReclusteredHits(slidingFitResultI, pClusterJ, caloHitsToAddI, caloHitsToRemoveJ);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackConsolidationAlgorithm::GetReclusteredHits(const TwoDSlidingFitResult& slidingFitResultI, const Cluster* pClusterJ,
    ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
    const Cluster* pClusterI(slidingFitResultI.GetCluster());

    CaloHitList associatedHits, caloHitListJ;
    pClusterJ->GetOrderedCaloHitList().GetCaloHitList(caloHitListJ);

    float minL(std::numeric_limits<float>::max());
    float maxL(std::numeric_limits<float>::max());

    // Loop over hits from shower clusters, and make associations with track clusters
    // (Determine if hits from shower clusters can be used to fill gaps in track cluster)
    //
    // Apply the following selection:
    //  rJ = candidate hit from shower cluster
    //  rI = nearest hit on track cluster
    //  rK = projection of shower hit onto track cluster
    //
    //                  o rJ
    //  o o o o o o - - x - - - - o o o o o o o
    //           rI    rK
    //
    //  Require: rJK < std::min(rCut, rIJ, rKI)

    for (CaloHitList::const_iterator iterJ = caloHitListJ.begin(), iterEndJ = caloHitListJ.end(); iterJ != iterEndJ; ++iterJ)
    {
        CaloHit* pCaloHitJ = *iterJ;

        try
        {
            const CartesianVector positionJ(pCaloHitJ->GetPositionVector());
            const CartesianVector positionI(LArClusterHelper::GetClosestPosition(positionJ, pClusterI));

            float rL(0.f), rT(0.f);
            CartesianVector positionK(0.f, 0.f, 0.f);
            slidingFitResultI.GetGlobalFitProjection(positionJ, positionK);
            slidingFitResultI.GetLocalPosition(positionK, rL, rT);

            const float rsqIJ((positionI - positionJ).GetMagnitudeSquared());
            const float rsqJK((positionJ - positionK).GetMagnitudeSquared());
            const float rsqKI((positionK - positionI).GetMagnitudeSquared());

            if (rsqJK < std::min(m_maxTransverseDisplacement * m_maxTransverseDisplacement, std::min(rsqIJ, rsqKI)))
            {
                if (associatedHits.empty())
                {
                    minL = rL;
                    maxL = rL;
                }
                else
                {
                    minL = std::min(minL, rL);
                    maxL = std::max(maxL, rL);
                }

                associatedHits.insert(pCaloHitJ);
            }
        }
        catch (StatusCodeException &)
        {
        }
    }

    const float associatedSpan(maxL - minL);
    const float associatedFraction(associatedHits.empty() ? 0.f : static_cast<float>(associatedHits.size()) / static_cast<float>(pClusterJ->GetNCaloHits()));

    if (associatedSpan > m_minAssociatedSpan || associatedFraction > m_minAssociatedFraction)
    {
        for (CaloHitList::const_iterator iterK = associatedHits.begin(), iterEndK = associatedHits.end(); iterK != iterEndK; ++iterK)
        {
            CaloHit* pCaloHit = *iterK;

            if (caloHitsToRemoveJ[pClusterJ].count(pCaloHit))
                continue;

            caloHitsToAddI[pClusterI].insert(pCaloHit);
            caloHitsToRemoveJ[pClusterJ].insert(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_maxTransverseDisplacement = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_minAssociatedSpan = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedSpan", m_minAssociatedSpan));

    m_minAssociatedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedFraction", m_minAssociatedFraction));

    return TwoDSlidingFitConsolidationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
