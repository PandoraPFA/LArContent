/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray track consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.h"



using namespace pandora;

namespace lar
{

StatusCode CosmicRayTrackConsolidationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));



    ClusterVector trackClusters, shortClusters;

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        const float thisLengthSquared(LArClusterHelper::GetLengthSquared(pCluster));

        if (thisLengthSquared < m_maxClusterLength * m_maxClusterLength)
            shortClusters.push_back(pCluster);

        if (thisLengthSquared > m_minTrackLength * m_minTrackLength)
            trackClusters.push_back(pCluster);
    }

    std::sort(trackClusters.begin(), trackClusters.end(), LArClusterHelper::SortByNHits);






    TwoDSlidingFitResultList slidingFitResultList;

    for (ClusterVector::const_iterator iter = trackClusters.begin(), iterEnd = trackClusters.end(); iter != iterEnd; ++iter)
    {
         TwoDSlidingFitResult slidingFitResult;
         LArClusterHelper::LArTwoDSlidingFit(*iter, m_halfWindowLayers, slidingFitResult);

         slidingFitResultList.push_back(slidingFitResult);
    }







    for (TwoDSlidingFitResultList::const_iterator iterI = slidingFitResultList.begin(), iterEndI = slidingFitResultList.end(); iterI != iterEndI; ++iterI)
    {
        const TwoDSlidingFitResult slidingFitResultI = *iterI;

        const Cluster* pClusterI = slidingFitResultI.GetCluster();
        const float thisLengthSquaredI(LArClusterHelper::GetLengthSquared(pClusterI));

        HitAssociationMap hitAssociationMap;

        for (ClusterVector::const_iterator iterJ = shortClusters.begin(), iterEndJ = shortClusters.end(); iterJ != iterEndJ; ++iterJ)
        {
            const Cluster* pClusterJ = *iterJ;
            const float thisLengthSquaredJ(LArClusterHelper::GetLengthSquared(pClusterJ));

            if (pClusterI == pClusterJ)
                continue;

            if (2.0 * thisLengthSquaredJ > thisLengthSquaredI)
                continue;

            this->GetAssociatedHits(slidingFitResultI, pClusterJ, hitAssociationMap);
        }

        ClusterList tempListI;  tempListI.insert((Cluster*)slidingFitResultI.GetCluster());
        ClusterList tempListJ;
        CaloHitList tempListK;

        for(HitAssociationMap::const_iterator iterJ = hitAssociationMap.begin(), iterEndJ = hitAssociationMap.end(); iterJ != iterEndJ; ++iterJ)
        {
            const Cluster* pClusterJ = iterJ->first;
            const CaloHitList &caloHitList = iterJ->second;

            tempListJ.insert((Cluster*)pClusterJ);
            tempListK.insert(caloHitList.begin(), caloHitList.end());
        }

        if (!tempListJ.empty())
        {
        PandoraMonitoringApi::VisualizeClusters(&tempListI, "Cluster", BLUE);
        PandoraMonitoringApi::VisualizeCaloHits(&tempListK, "AssociatedHits", RED);
        PandoraMonitoringApi::ViewEvent();
        }

    }


    // RE-CLUSTERING GOES HERE


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackConsolidationAlgorithm::GetAssociatedHits(const TwoDSlidingFitResult& slidingFitResultI, const Cluster* pClusterJ,
    HitAssociationMap &hitAssociationMap) const
{

    const Cluster* pClusterI(slidingFitResultI.GetCluster());
    const OrderedCaloHitList &orderedCaloHitList(pClusterJ->GetOrderedCaloHitList());

    CaloHitList associatedHits;

    float minL(std::numeric_limits<float>::max());
    float maxL(std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter1 = orderedCaloHitList.begin(), iterEnd1 = orderedCaloHitList.end(); iter1 != iterEnd1; ++iter1)
    {
        for (CaloHitList::const_iterator iter2 = iter1->second->begin(), iterEnd2 = iter1->second->end(); iter2 != iterEnd2; ++iter2)
        {
            CaloHit* pCaloHitJ = *iter2;

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
    }

    const float associatedSpan(maxL - minL);
    const float associatedFraction(associatedHits.empty() ? 0.f : static_cast<float>(associatedHits.size()) / static_cast<float>(pClusterJ->GetNCaloHits()));

    if (associatedSpan > m_minAssociatedSpan || associatedFraction > m_minAssociatedFraction)
        hitAssociationMap[pClusterJ].insert(associatedHits.begin(), associatedHits.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minTrackLength = 7.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLength", m_minTrackLength));

    m_maxClusterLength = 15.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterLength", m_maxClusterLength));

    m_maxTransverseDisplacement = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_minAssociatedSpan = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedSpan", m_minAssociatedSpan));

    m_minAssociatedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedFraction", m_minAssociatedFraction));

    m_halfWindowLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
