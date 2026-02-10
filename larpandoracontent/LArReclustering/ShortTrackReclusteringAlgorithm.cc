/**
 *  @file   larpandoracontent/LArReclustering/ShortTrackReclusteringAlgorithm.cc
 *
 *  @brief  Tries to identify short tracks that are either lost, or under-clustered and recover missing hits.
 *          This algorithm looks at all three views of a PFO to try to identify if one or more views exhibit step changes in ADC values that
 *          might be indicative of a decay or intelastic interaction. If such a change occurs it looks for consistency across views and also
 *          examines nearby unclustered hits, or hits in nearby clusters to see if a more cohenrent clustering can be identified.
 *          If these various conditions are met, then reclustering is performed.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArReclustering/ShortTrackReclusteringAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

ShortTrackReclusteringAlgorithm::ShortTrackReclusteringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShortTrackReclusteringAlgorithm::Run()
{
    // Note: This may become a tool within the reclustering paradigm, but I'll need to chat with Alex about that

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_IF(STATUS_CODE_SUCCESS, !this->GetList(m_caloHitListName, pCaloHitList));

    ViewToHitsMap viewToUnclusteredHitsMap;
    this->CollectUnclusteredHits(*pCaloHitList, viewToUnclusteredHitsMap);

    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_IF(STATUS_CODE_SUCCESS, !this->GetList(m_pfoListName, pPfoList));

    ViewToClustersMap viewToClustersMap;
    ClusterToPfoMap clusterToPfoMap;
    this->CollectClusters(*pPfoList, viewToClustersMap, clusterToPfoMap);

    // Identify the end points of clusters in PFOs
    for (const auto &[pCluster, pPfo] : clusterToPfoMap)
    {
        CaloHitList clusterHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHitList);
        if (clusterHitList.empty())
            continue;

        CaloHitVector clusterHits(clusterHitList.begin(), clusterHitList.end());
        HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
        const float nHalfWindow{2};
        const float pitch{LArGeometryHelper::GetWirePitch(this->GetPandora(), view)};
        LArPointingCluster pointingCluster(pCluster, nHalfWindow, pitch);

        CaloHitVector forwardHits, backwardHits;
        this->OrderHitsRelativeToVertex(clusterHits, pointingCluster.GetInnerVertex(), forwardHits);
        this->OrderHitsRelativeToVertex(clusterHits, pointingCluster.GetOuterVertex(), backwardHits);

        IntVector discontinuities;
        this->GetStableAdcDiscontinuities(forwardHits, discontinuities);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool ShortTrackReclusteringAlgorithm::GetList(const std::string &listName, const T *&pList) const
{
    PandoraContentApi::GetList(*this, listName, pList);
    return pList && !pList->empty();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::CollectUnclusteredHits(const CaloHitList &caloHitList, ViewToHitsMap &viewToUnclusteredHitsMap) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            viewToUnclusteredHitsMap[pCaloHit->GetHitType()].insert(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::CollectClusters(const PfoList &pfoList, ViewToClustersMap &viewToClustersMap, ClusterToPfoMap &clusterToPfoMap) const
{
    for (const Pfo *const pPfo : pfoList)
    {
        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            ClusterList pfoClusterList;
            LArPfoHelper::GetClusters(pPfo, view, pfoClusterList);
            for (const Cluster *const pCluster : pfoClusterList)
            {
                viewToClustersMap[view].emplace_back(pCluster);
                clusterToPfoMap[pCluster] = pPfo;
            }
        }
    }
}

void ShortTrackReclusteringAlgorithm::OrderHitsRelativeToVertex(const CaloHitVector &clusterHits, const LArPointingCluster::Vertex &vertex,
    CaloHitVector &orderedHits) const
{
    const CartesianVector &position(vertex.GetPosition());
    const CartesianVector &direction(vertex.GetDirection());

    // Get the hit-vertex vector and project it onto the pointing cluster vector
    FloatVector projections;
    for (const CaloHit *const pCaloHit : clusterHits)
    {
        const CartesianVector hitDirection(pCaloHit->GetPositionVector() - position);
        projections.emplace_back(direction.GetDotProduct(hitDirection));
    }

    // Sort the hits according to their projections
    std::vector<size_t> hitIndices(clusterHits.size());
    std::iota(hitIndices.begin(), hitIndices.end(), 0);
    std::sort(hitIndices.begin(), hitIndices.end(), [&projections](size_t i, size_t j) { return projections[i] < projections[j]; });

    for (size_t i = 0; i < hitIndices.size(); ++i)
        orderedHits.emplace_back(clusterHits[hitIndices[i]]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::GetAdcMovingAverage(const FloatVector &adcs, FloatVector &movingAverage, FloatVector &movingVariance,
    const size_t window) const
{
    // Compute backward looking moving average and variance. Peforms rolling update for efficiency
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, window == 0);
    const size_t nHits{adcs.size()};
    movingAverage.resize(nHits);
    movingVariance.resize(nHits);

    float S{0.f}, V{0.f};
    size_t count{0};
    // initial window
    for (size_t i = 0; i < std::min(nHits, window); ++i)
    {
        S += adcs[i];
        V += adcs[i] * adcs[i];
        ++count;

        movingAverage[i] = S / count;
        movingVariance[i] = count > 1 ? (V - S*S / count) / (count - 1) : 0.0; // sample variance
    }
    // rolling update
    for (size_t i = window; i < nHits; ++i)
    {
        const float o{adcs[i - window]};
        const float n{adcs[i]};

        S += n - o;
        V += n*n - o*o;

        movingAverage[i] = S / count;
        movingVariance[i] = (V - S*S / count) / (count - 1); // sample variance
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::GetStableAdcDiscontinuities(const pandora::CaloHitVector &hits, pandora::IntVector &discontinuities,
    const size_t window) const
{
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, window == 0);
    FloatVector normalizedAdc, movingAverage, movingVariance;
    this->NormalizeAdc(hits, normalizedAdc);
    this->GetAdcMovingAverage(normalizedAdc, movingAverage, movingVariance, window);
    (void)discontinuities;

    /*
    // Look for step changes in ADC
    std::cout << "Cluster with " << clusterHits.size() << std::endl;
    int i{0};
    for (const CaloHit *const pCaloHit : hits)
    {
        const CartesianVector position(pCaloHit->GetPositionVector());
        const float adc{pCaloHit->GetInputEnergy()};
        std::cout << "Hit (" << position.GetX() << ", " << position.GetZ() << ") = " << adc << "(" << forwardAdcs[i] << ")" << std::endl;
        if (i >= 3 && (forwardAdcs[i] > 2 * forwardAdcs[i - 3] || forwardAdcs[i] < 0.5f * forwardAdcs[i - 3]))
            std::cout << "  Potential step change in ADC values!" << std::endl;
        ++i;
    }
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::NormalizeAdc(const pandora::CaloHitVector &hits, pandora::FloatVector &normalizedAdc) const
{
    normalizedAdc.clear();
    const size_t mid{hits.size() / 2};

    FloatVector adcs;
    for (const CaloHit *const pCaloHit : hits)
        adcs.emplace_back(pCaloHit->GetInputEnergy());
    // Find median via nth_element - generally faster than sorting (O(n) vs O(n log n))
    float median{0.f};
    if (mid % 2 == 0)
    {
        std::nth_element(adcs.begin(), adcs.begin() + mid, adcs.end());
        const float upper{adcs[mid]};
        std::nth_element(adcs.begin(), adcs.begin() + mid - 1, adcs.end());
        const float lower{adcs[mid - 1]};
        median = 0.5f * (lower + upper);
    }
    else
    {
        std::nth_element(adcs.begin(), adcs.begin() + mid, adcs.end());
        median = adcs[mid];
    }

    // This shouldn't be an issue, but just in case
    if (median < 1.f)
        median = 1.f;
    for (const CaloHit *const pCaloHit : hits)
        normalizedAdc.emplace_back(pCaloHit->GetInputEnergy() / median);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShortTrackReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
