/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the track consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackConsolidationAlgorithm::TrackConsolidationAlgorithm() :
    m_maxTransverseDisplacement(1.f), m_minAssociatedSpan(1.f), m_minAssociatedFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackConsolidationAlgorithm::GetReclusteredHits(const TwoDSlidingFitResultList &slidingFitResultListI,
    const ClusterVector &showerClustersJ, ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
    for (const TwoDSlidingFitResult &slidingFitResultI : slidingFitResultListI)
    {
        const Cluster *const pClusterI = slidingFitResultI.GetCluster();
        const float thisLengthSquaredI(LArClusterHelper::GetLengthSquared(pClusterI));

        for (const Cluster *const pClusterJ : showerClustersJ)
        {
            const float thisLengthSquaredJ(LArClusterHelper::GetLengthSquared(pClusterJ));

            if (pClusterI == pClusterJ)
                continue;

            if (2.f * thisLengthSquaredJ > thisLengthSquaredI)
                continue;

            this->GetReclusteredHits(slidingFitResultI, pClusterJ, caloHitsToAddI, caloHitsToRemoveJ);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackConsolidationAlgorithm::GetReclusteredHits(const TwoDSlidingFitResult &slidingFitResultI, const Cluster *const pClusterJ,
    ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
    const Cluster *const pClusterI(slidingFitResultI.GetCluster());
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), LArClusterHelper::GetClusterHitType(pClusterI))};
    const float maxTransverseDisplacementAdjusted{ratio * m_maxTransverseDisplacement};
    const float minAssociatedSpanAdjusted{ratio * m_minAssociatedSpan};

    CaloHitList associatedHits, caloHitListJ;
    pClusterJ->GetOrderedCaloHitList().FillCaloHitList(caloHitListJ);

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
        const CaloHit *const pCaloHitJ = *iterJ;

        const CartesianVector positionJ(pCaloHitJ->GetPositionVector());
        const CartesianVector positionI(LArClusterHelper::GetClosestPosition(positionJ, pClusterI));

        float rL(0.f), rT(0.f);
        CartesianVector positionK(0.f, 0.f, 0.f);

        if (STATUS_CODE_SUCCESS != slidingFitResultI.GetGlobalFitProjection(positionJ, positionK))
            continue;

        slidingFitResultI.GetLocalPosition(positionK, rL, rT);

        const float rsqIJ((positionI - positionJ).GetMagnitudeSquared());
        const float rsqJK((positionJ - positionK).GetMagnitudeSquared());
        const float rsqKI((positionK - positionI).GetMagnitudeSquared());

        if (rsqJK < std::min(maxTransverseDisplacementAdjusted * maxTransverseDisplacementAdjusted, std::min(rsqIJ, rsqKI)))
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

            associatedHits.push_back(pCaloHitJ);
        }
    }

    const float associatedSpan(maxL - minL);
    const float associatedFraction(
        associatedHits.empty() ? 0.f : static_cast<float>(associatedHits.size()) / static_cast<float>(pClusterJ->GetNCaloHits()));

    if (associatedSpan > minAssociatedSpanAdjusted || associatedFraction > m_minAssociatedFraction)
    {
        for (CaloHitList::const_iterator iterK = associatedHits.begin(), iterEndK = associatedHits.end(); iterK != iterEndK; ++iterK)
        {
            const CaloHit *const pCaloHit = *iterK;
            const CaloHitList &caloHitList(caloHitsToRemoveJ[pClusterJ]);

            if (caloHitList.end() != std::find(caloHitList.begin(), caloHitList.end(), pCaloHit))
                continue;

            caloHitsToAddI[pClusterI].push_back(pCaloHit);
            caloHitsToRemoveJ[pClusterJ].push_back(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinAssociatedSpan", m_minAssociatedSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinAssociatedFraction", m_minAssociatedFraction));

    return TwoDSlidingFitConsolidationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
