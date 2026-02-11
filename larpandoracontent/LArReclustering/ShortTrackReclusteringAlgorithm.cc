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

    FloatVector rollingWindow;
    float S{0.f}, V{0.f};
    size_t count{0};
    // initial window
    for (size_t i = 0; i < std::min(nHits, window); ++i)
    {
        rollingWindow.emplace_back(adcs[i]);
        S += adcs[i];
        V += adcs[i] * adcs[i];
        ++count;

        movingAverage[i] = static_cast<float>(this->GetMedian(rollingWindow));
        movingVariance[i] = count > 1 ? (V - S*S / count) / (count - 1) : 0.0; // sample variance
    }
    // rolling update
    for (size_t i = window; i < nHits; ++i)
    {
        rollingWindow.erase(rollingWindow.begin());
        rollingWindow.emplace_back(adcs[i]);
        const float o{adcs[i - window]};
        const float n{adcs[i]};

        S += n - o;
        V += n*n - o*o;

        movingAverage[i] = static_cast<float>(this->GetMedian(rollingWindow));
        movingVariance[i] = (V - S*S / count) / (count - 1); // sample variance
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
double ShortTrackReclusteringAlgorithm::GetMedian(const std::vector<T> &values) const
{
    std::vector<T> copy(values);
    const size_t mid{copy.size() / 2};

    if (mid % 2 == 0)
    {
        std::nth_element(copy.begin(), copy.begin() + mid, copy.end());
        const double upper{copy[mid]};
        std::nth_element(copy.begin(), copy.begin() + mid - 1, copy.end());
        const double lower{copy[mid - 1]};
        return 0.5 * (lower + upper);
    }
    else
    {
        std::nth_element(copy.begin(), copy.begin() + mid, copy.end());
        return copy[mid];
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::GetStableAdcDiscontinuities(const pandora::CaloHitVector &hits, pandora::IntVector &discontinuities,
    const size_t window) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, -1, 1));
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, window == 0);
    FloatVector normalizedAdc, movingAverage, movingVariance;
    this->NormalizeAdc(hits, normalizedAdc);
    this->GetAdcMovingAverage(normalizedAdc, movingAverage, movingVariance, window);

    // Look for step changes in ADC
    for (size_t i = window; i < hits.size(); ++i)
    {
        const float current{movingAverage[i]}, previous{movingAverage[i - 1] > 0 ? movingAverage[i - 1] : 1.f},
            previous2{((i >= 2) && (movingAverage[i - 2] > 0)) ? movingAverage[i - 2] : 1.f};
        const float ratio{std::max(current / previous, current / previous2)};

        if (ratio > 2.f)
        {
            // Possible discontinuity to more highly ionising particle, check local variance before step
            int nLow{0}, nHigh{0};
            for (size_t j = i - window; j < i; ++j)
            {
                if (movingVariance[j] < 0.3f) // Arbitrary threshold for now, but probably the right scale
                    ++nLow;
                else
                    ++nHigh;
            }
            if (nHigh >= nLow)
                continue;
            // Check for possible Bragg peak if the step is to more highly ionising particle
            const size_t start{i}, end{std::min(i + 2 * window - 1, hits.size() - 1)};
            if (this->IsBraggPeak(hits, start, end))
                continue;
            discontinuities.emplace_back(i);
        }
        else
        {
            continue;
        }
    }
    std::cout << "Found " << discontinuities.size() << " discontinuities" << std::endl;
    for (const int index : discontinuities)
    {
        const CartesianVector &position(hits[index]->GetPositionVector());
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "d", RED, 2));
    }
    if (!discontinuities.empty())
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortTrackReclusteringAlgorithm::NormalizeAdc(const pandora::CaloHitVector &hits, pandora::FloatVector &normalizedAdc) const
{
    normalizedAdc.clear();
    FloatVector adcs;
    for (const CaloHit *const pCaloHit : hits)
        adcs.emplace_back(pCaloHit->GetInputEnergy());
    // Find median via nth_element - generally faster than sorting (O(n) vs O(n log n))
    float median{static_cast<float>(this->GetMedian(adcs))};

    // This shouldn't be an issue, but just in case
    if (median < 1.)
        median = 1.;
    for (const CaloHit *const pCaloHit : hits)
        normalizedAdc.emplace_back(pCaloHit->GetInputEnergy() / median);
}

bool ShortTrackReclusteringAlgorithm::IsBraggPeak(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    const float linearSlopeScore{this->GetLinearSlopeScore(hits, start, end)};
    // Curvature needs at least 5 hits for calculation
    const float curvatureScore{(end - start) >= 4 ? this->GetQuadraticCurvatureScore(hits, start, end) : 0.f};
    const float contrastScore{this->GetContrastScore(hits, start, end)};
    const float monotonicityScore{this->GetMonotonicityScore(hits, start, end)};

    int consensus{0};
    consensus += linearSlopeScore > 0.16f;
    consensus += curvatureScore > 0.f;
    consensus += contrastScore > 1.4f;
    consensus += monotonicityScore > 0.5f;

    return consensus >= 3;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetLinearSlopeScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    FloatVector adcs;
    for (size_t i = start; i <= end; ++i)
        adcs.emplace_back(hits[i]->GetInputEnergy());
    const size_t mid{adcs.size() / 2};
    std::nth_element(adcs.begin(), adcs.begin() + mid, adcs.end());
    const float median{adcs[mid]};

    const float dx{hits[end]->GetPositionVector().GetX() - hits[start]->GetPositionVector().GetX()};
    const float dz{hits[end]->GetPositionVector().GetZ() - hits[start]->GetPositionVector().GetZ()};
    const float dr{std::sqrt(dx*dx + dz*dz)};

    return (dr > 0) ? (hits[end]->GetInputEnergy() - hits[start]->GetInputEnergy()) / (median * dr) : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetQuadraticCurvatureScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    // Use standard three-point non-uniform finite difference for second derivative
    float curvature{0.f};
    for (size_t i = start + 1; i < end; ++i)
    {
        const float a_m{hits[i - 1]->GetInputEnergy()};
        const float a_0{hits[i]->GetInputEnergy()};
        const float a_p{hits[i + 1]->GetInputEnergy()};
        const float dx_m{hits[i]->GetPositionVector().GetX() - hits[i - 1]->GetPositionVector().GetX()};
        const float dz_m{hits[i]->GetPositionVector().GetZ() - hits[i - 1]->GetPositionVector().GetZ()};
        const float dr_m{std::sqrt(dx_m*dx_m + dz_m*dz_m)};
        const float dx_p{hits[i + 1]->GetPositionVector().GetX() - hits[i]->GetPositionVector().GetX()};
        const float dz_p{hits[i + 1]->GetPositionVector().GetZ() - hits[i]->GetPositionVector().GetZ()};
        const float dr_p{std::sqrt(dx_p*dx_p + dz_p*dz_p)};

        curvature += (2.f / (dr_m + dr_p)) * (((a_p - a_0) / dr_p) - ((a_0 - a_m) / dr_m));
    }

    return curvature / (end - start - 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetContrastScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    const size_t mid{(start + end) / 2};
    float muHead{0}, muTail{0};
    for (size_t i = start; i <= mid; ++i)
        muHead += hits[i]->GetInputEnergy();
    muHead /= (mid - start + 1);

    for (size_t i = mid + 1; i <= end; ++i)
        muTail += hits[i]->GetInputEnergy();
    muTail /= (end - mid);
    
    return (muHead > 0.f) ? muTail / muHead : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShortTrackReclusteringAlgorithm::GetMonotonicityScore(const pandora::CaloHitVector &hits, const size_t start, const size_t end) const
{
    if (start == end)
        return 0.f;
    float monotonicity{0.f};
    for (size_t i = start + 1; i <= end; ++i)
    {
        const float delta{hits[i]->GetInputEnergy() - hits[i - 1]->GetInputEnergy()};
        // We only care about monotonically increasing cases here
        monotonicity += delta > 0.f ? 1.f : 0.f;
    }
    
    return monotonicity / (end - start);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShortTrackReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
