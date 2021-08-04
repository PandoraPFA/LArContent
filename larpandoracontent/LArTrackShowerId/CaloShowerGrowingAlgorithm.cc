/**
 *  @file   larpandoracontent/LArTrackShowerId/CaloShowerGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the shower growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArTrackShowerId/CaloShowerGrowingAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

CaloShowerGrowingAlgorithm::CaloShowerGrowingAlgorithm() :
    m_minCaloHitsForSeed(5),
    m_radiationLength(14.f),
    m_moliereRadius(9.043f),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    std::copy(pClusterList->begin(), pClusterList->end(), std::back_inserter(clusterVector));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterList seedClusterList;
    this->GetSeedClusters(clusterVector, seedClusterList);
    std::map<const Cluster *, float> seedToChi2Map;
    std::map<const Cluster *, ClusterList> seedToAssociationMap;
    std::map<const Cluster *, ClusterList> seedIntersectionMap;
    for (const Cluster *pSeed : seedClusterList)
    {
        ClusterList associatedClusterList;
        const Bounds &bounds{this->GetSeedBounds(pSeed)};
        if (m_visualize)
        {
            CaloHitList caloHits;
            LArClusterHelper::GetAllHits(pSeed, caloHits);
            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "seed", GRAY));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &bounds.m_a, &bounds.m_b, "ab", RED, 1, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &bounds.m_b, &bounds.m_c, "bc", GREEN, 1, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &bounds.m_c, &bounds.m_d, "cd", BLUE, 1, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &bounds.m_d, &bounds.m_a, "da", BLACK, 1, 1));
        }
        for (const Cluster *pCluster : clusterVector)
        {
            if (pSeed == pCluster)
                continue;
            CaloHitList testCaloHits;
            LArClusterHelper::GetAllHits(pCluster, testCaloHits);
            bool isContained{false};
            for (const CaloHit *pCaloHit : testCaloHits)
            {
                isContained = bounds.Contains(pCaloHit->GetPositionVector());
                if (isContained)
                    break;
            }
            if (isContained)
                associatedClusterList.emplace_back(pCluster);
        }
        if (m_visualize)
        {
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &associatedClusterList, "associated", BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
        seedToAssociationMap[pSeed] = ClusterList();
        seedToChi2Map[pSeed] = this->AssessAssociation(pSeed, associatedClusterList, seedToAssociationMap[pSeed]);
        // Need to think about whether or not an additional pass to potentially merge seeds is wanted
    }

    for (auto iter1 = seedClusterList.begin(); iter1 != seedClusterList.end(); ++iter1)
    {
        const Cluster *pSeed1{*iter1};
        // Ensure all seeds have an entry in the intersection map, even if they have exclusive associations
        if (seedIntersectionMap.find(pSeed1) == seedIntersectionMap.end())
            seedIntersectionMap[pSeed1] = ClusterList();
        const ClusterList &list1{seedToAssociationMap[pSeed1]};
        for (auto iter2 = std::next(iter1); iter2 != seedClusterList.end(); ++iter2)
        {
            const Cluster *pSeed2{*iter2};
            const ClusterList &list2{seedToAssociationMap[pSeed2]};
            for (const Cluster *pCluster : list1)
            {
                if (std::find(list2.begin(), list2.end(), pCluster) != list2.end())
                {   // Shared cluster association between seeds, can't merge both
                    seedIntersectionMap[pSeed1].emplace_back(pSeed2);
                    seedIntersectionMap[pSeed2].emplace_back(pSeed1);
                    break;
                }
            }
        }
    }

    ClusterList consideredSeeds;
    for (const auto &[ key, intersections ] : seedIntersectionMap)
    {
        const Cluster *pBestSeed{nullptr};
        float bestChi2{std::numeric_limits<float>::max()};
        if (std::find(consideredSeeds.begin(), consideredSeeds.end(), key) == consideredSeeds.end())
        {
            pBestSeed = key;
            bestChi2 = seedToChi2Map[key];
            consideredSeeds.emplace_back(key);
        }

        for (const Cluster *pSeed : intersections)
        {   // Check if another seed with shared associations would be better
            if (std::find(consideredSeeds.begin(), consideredSeeds.end(), pSeed) == consideredSeeds.end())
            {
                if (seedToChi2Map[pSeed] < bestChi2)
                {
                    pBestSeed = pSeed;
                    bestChi2 = seedToChi2Map[pSeed];
                }
                consideredSeeds.emplace_back(pSeed);
            }
        }

        if (pBestSeed)
        {   // Populate the merge map with the best seed and its associations
            for (const Cluster *pCluster : seedToAssociationMap[pBestSeed])
            {
                clusterMergeMap[pBestSeed].emplace_back(pCluster);
                clusterMergeMap[pCluster].emplace_back(pBestSeed);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GetSeedClusters(const ClusterVector &clusterVector, ClusterList &seedClusterList) const
{
    for (const Cluster *const pCluster : clusterVector)
    {
        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsForSeed)
            continue;

        seedClusterList.emplace_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloShowerGrowingAlgorithm::Bounds CaloShowerGrowingAlgorithm::GetSeedBounds(const Cluster *pSeed) const
{
    CaloHitList caloHits;
    pSeed->GetOrderedCaloHitList().FillCaloHitList(caloHits);

    // Get the eigen vectors for this cluster
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(caloHits, centroid, eigenValues, eigenVectors);

    // Project the extremal hits of the cluster onto the principal axis
    CartesianVector axis0(0.f, 0.f, 0.f), axis1(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pSeed, axis0, axis1);
    axis0 -= centroid;
    axis1 -= centroid;
    // ATTN: PCA axis is a unit vector, so no need to correct for length in dot products
    const CartesianVector lAxis(eigenVectors[0]);
    const CartesianVector tAxis(-lAxis.GetZ() * m_moliereRadius, lAxis.GetY() * m_moliereRadius, lAxis.GetX() * m_moliereRadius);
    const float proj0{axis0.GetDotProduct(lAxis)};
    const float proj1{axis1.GetDotProduct(lAxis)};
    const float ext0{proj0 >= proj1 ? proj0 + m_radiationLength : proj0 - m_radiationLength};
    const float ext1{proj0 >= proj1 ? proj1 - m_radiationLength : proj1 + m_radiationLength};
    const CartesianVector a{lAxis * ext0 + centroid - tAxis};
    const CartesianVector b{lAxis * ext1 + centroid - tAxis};
    const CartesianVector c{lAxis * ext1 + centroid + tAxis};
    const CartesianVector d{lAxis * ext0 + centroid + tAxis};

    return Bounds(a, b, c, d);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CaloShowerGrowingAlgorithm::AssessAssociation(const Cluster *pSeed, const ClusterList &associatedClusterList, ClusterList &showerClusterList) const
{
    CaloHitList caloHits;
    LArClusterHelper::GetAllHits(pSeed, caloHits);
    for (const Cluster *pCluster : associatedClusterList)
        LArClusterHelper::GetAllHits(pCluster, caloHits);

    float chi2Baseline{this->GetShowerProfileChi2(caloHits)};

    std::map<const Cluster*, bool> dropoutClusterMap;
    for (const Cluster *pDropout : associatedClusterList)
    {
        caloHits.clear();
        LArClusterHelper::GetAllHits(pSeed, caloHits);
        for (const Cluster *pCluster : associatedClusterList)
        {
            if (pCluster == pDropout || dropoutClusterMap.find(pCluster) != dropoutClusterMap.end())
                continue;
            LArClusterHelper::GetAllHits(pCluster, caloHits);
        }
        const float chi2{this->GetShowerProfileChi2(caloHits)};
        if (chi2 < chi2Baseline * 0.8f)
        {
            dropoutClusterMap[pDropout] = true;
            chi2Baseline = chi2;
        }
    }

    for (const Cluster *pCluster : associatedClusterList)
    {
        if (dropoutClusterMap.find(pCluster) != dropoutClusterMap.end())
            continue;
        showerClusterList.emplace_back(pCluster);
    }

    return chi2Baseline;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CaloShowerGrowingAlgorithm::GetShowerProfileChi2(const pandora::CaloHitList &caloHitList) const
{
    // Get the eigen vectors for this collection of hits
    CartesianVector origin(0.f, 0.f, 0.f);
    CartesianVector dir(0.f, 0.f, 0.f);
    // ATTN: PCA axis is a unit vector, so no need to correct for length in dot products
    this->GetProjectionAxis(caloHitList, origin, dir);

    if (m_visualize)
    {
        CartesianVector source{origin - dir * 20.f}, sink{origin + dir * 20.f};
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "associated", BLUE));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &source, &sink, "axis", RED, 1, 1));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    // Project all hits onto the direction axis
    FloatVector projVector, energies;
    float min{std::numeric_limits<float>::max()}, max{-std::numeric_limits<float>::max()};
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float proj{dir.GetDotProduct(pCaloHit->GetPositionVector() - origin)};
        if (proj < min)
            min = proj;
        if (proj > max)
            max = proj;
        projVector.emplace_back(proj);
        energies.emplace_back(pCaloHit->GetElectromagneticEnergy());
    }
    // Adjust min and max to avoid underflow/overflow issues
    min -= std::numeric_limits<float>::epsilon();
    max += std::numeric_limits<float>::epsilon();
    float binSize{0.5f * m_radiationLength};
    if (binSize < std::numeric_limits<float>::epsilon())
    {
        std::cout << "Error: radiation length implausibly small" << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
    const int N{static_cast<int>(std::ceil((max - min) / binSize))};
    if (N <= 1)
        return std::numeric_limits<float>::max();
    FloatVector binCentres(N);
    for (int i = 0; i < N; ++i)
        binCentres[i] = 0.5f * (2 * i + 1) * binSize;

    // Reconstructed energy profile
    FloatVector binEnergiesForward(N), binEnergiesBackward(N);
    for (unsigned int i = 0; i < projVector.size(); ++i)
    {
        const int bin{static_cast<int>((projVector[i] - min) / binSize)};
        binEnergiesForward[bin] += energies[i];
        binEnergiesBackward[N - 1 - bin] += energies[i];
    }

    // Theoretical energy profiles
    const float e0{1000.f * std::accumulate(binEnergiesForward.begin(), binEnergiesForward.end(), 0.f)};
    FloatVector photonEnergies(N), electronEnergies(N);
    this->GetPhotonLongitudinalEnergyProfile(e0, binCentres, binSize, photonEnergies);
    this->GetElectronLongitudinalEnergyProfile(e0, binCentres, binSize, electronEnergies);

    // Assess profiles
    float photonChi2Forward{0.f}, photonChi2Backward{0.f}, electronChi2Forward{0.f}, electronChi2Backward{0.f};
    for (int i = 0; i < N; ++i)
    {
        float dE{binEnergiesForward[i] - photonEnergies[i]};
        photonChi2Forward += (dE * dE) / photonEnergies[i];
        dE = binEnergiesBackward[i] - photonEnergies[i];
        photonChi2Backward += dE * dE / photonEnergies[i];
        dE = binEnergiesForward[i] - electronEnergies[i];
        electronChi2Forward += (dE * dE) / electronEnergies[i];
        dE = binEnergiesBackward[i] - electronEnergies[i];
        electronChi2Backward += dE * dE / electronEnergies[i];
    }

    photonChi2Forward /= N - 1;
    photonChi2Backward /= N - 1;
    electronChi2Forward /= N - 1;
    electronChi2Backward /= N - 1;

    return std::min({photonChi2Forward, photonChi2Backward, electronChi2Forward, electronChi2Backward});
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GetProjectionAxis(const CaloHitList &caloHitList, CartesianVector &origin, CartesianVector &dir) const
{
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(caloHitList, origin, eigenValues, eigenVectors);

    dir = eigenVectors[0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GetPhotonLongitudinalEnergyProfile(const float e0, const FloatVector &positions, const float binSize,
    FloatVector &fractionalEnergies) const
{
    this->GetLongitudinalEnergyProfile(e0, positions, +0.5f, binSize, fractionalEnergies);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GetElectronLongitudinalEnergyProfile(const float e0, const FloatVector &positions, const float binSize,
    FloatVector &fractionalEnergies) const
{
    this->GetLongitudinalEnergyProfile(e0, positions, -0.5f, binSize, fractionalEnergies);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GetLongitudinalEnergyProfile(const float e0, const FloatVector &positions, const float cj, const float binSize,
    FloatVector &fractionalEnergies) const
{
    const unsigned int N{static_cast<unsigned int>(positions.size())};
    const float eCrit{32.84f};
    const float b{0.5f};
    float a{0.f};
    if (e0 > eCrit)
        a = (1 + b * cj) + b * std::log(e0 / eCrit);
    else
        a = 1 + b * cj;
    float gammaA{std::tgamma(a)};
    for (unsigned int i = 0; i < N; ++i)
    {
        const float t{positions[i]};
        fractionalEnergies[i] = b * std::pow(b * t * binSize, a - 1) * std::exp(-b * t) / gammaA;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsForSeed", m_minCaloHitsForSeed));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RadiationLength", m_radiationLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MoliereRadius", m_moliereRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CaloShowerGrowingAlgorithm::Bounds::Bounds(const CartesianVector &tl, const CartesianVector &bl, const CartesianVector &br, const CartesianVector &tr) :
    m_a{tl},
    m_b{bl},
    m_c{br},
    m_d{tr},
    m_ab{m_b - m_a},
    m_ac{m_c - m_a},
    m_ad{m_d - m_a}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CaloShowerGrowingAlgorithm::Bounds::Contains(const CartesianVector &point) const
{
    const CartesianVector &ap(point - m_a);

    return this->Contains(m_ab, m_ac, ap) || this->Contains(m_ac, m_ad, ap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CaloShowerGrowingAlgorithm::Bounds::Contains(const CartesianVector &ab, const CartesianVector &ac, const CartesianVector &ap) const
{
    const float dotABAB{ab.GetDotProduct(ab)};
    const float dotABAP{ab.GetDotProduct(ap)};
    const float dotACAC{ac.GetDotProduct(ac)};
    const float dotACAB{ac.GetDotProduct(ab)};
    const float dotACAP{ac.GetDotProduct(ap)};

    const float denom{dotACAC * dotABAB - dotACAB * dotACAB};
    if (std::fabs(denom) < std::numeric_limits<float>::epsilon())
        return false;
    const float invDenom{1.f / denom};
    const float u{(dotABAB * dotACAP - dotACAB * dotABAP) * invDenom};
    const float v{(dotACAC * dotABAP - dotACAB * dotACAP) * invDenom};

    return (u >= 0) && (v >= 0) && ((u + v) <= 1);
}

} // namespace lar_content

