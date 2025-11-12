/**
 *  @file   larpandoradlcontent/LArTwoDReco/DlTrackShowerStreamSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArTwoDReco/DlTrackShowerStreamSelectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <numeric>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

StatusCode DlTrackShowerStreamSelectionAlgorithm::AllocateToStreams(const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList{pCluster->GetOrderedCaloHitList()};
    CaloHitList caloHits;
    orderedCaloHitList.FillCaloHitList(caloHits);
    const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
    caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
    FloatVector trackLikelihoods;
    try
    {
        for (const CaloHit *pCaloHit : caloHits)
        {
            const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
            if (pLArCaloHit)
            {
                const float pTrack{pLArCaloHit->GetTrackProbability()};
                const float pShower{pLArCaloHit->GetShowerProbability()};
                if ((pTrack + pShower) > std::numeric_limits<float>::epsilon())
                    trackLikelihoods.emplace_back(pTrack / (pTrack + pShower));
            }
        }

        const unsigned long N{trackLikelihoods.size()};
        if (N > 0)
        {
            float mean{std::accumulate(std::begin(trackLikelihoods), std::end(trackLikelihoods), 0.f) / N};
            if (mean >= 0.5f)
                m_clusterListMap.at(m_trackListName).emplace_back(pCluster);
            else
                m_clusterListMap.at(m_showerListName).emplace_back(pCluster);
        }
    }
    catch (const StatusCodeException &)
    {
    }

    return STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlTrackShowerStreamSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, StreamSelectionAlgorithm::ReadSettings(xmlHandle));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackListName", m_trackListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerListName", m_showerListName));

    m_listNames.emplace_back(m_trackListName);
    m_listNames.emplace_back(m_showerListName);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
