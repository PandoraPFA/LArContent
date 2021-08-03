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

StatusCode CaloShowerGrowingAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        try
        {
            const ClusterList *pClusterList{nullptr};
            PANDORA_RETURN_RESULT_IF_AND_IF(
                STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

            if (!pClusterList || pClusterList->empty())
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "CaloShowerGrowingAlgorithm: unable to find cluster list " << clusterListName << std::endl;

                continue;
            }

            this->GrowShowers(*pClusterList);
        }
        catch (StatusCodeException &statusCodeException)
        {
            throw statusCodeException;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GrowShowers(const pandora::ClusterList &clusterList) const
{
    ClusterList seedClusterList;
    this->GetSeedClusters(clusterList, seedClusterList);
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
        for (const Cluster *pCluster : clusterList)
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
        ClusterList showerClusterList;
        this->AssessAssociation(pSeed, associatedClusterList, showerClusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GetSeedClusters(const ClusterList &clusterList, ClusterList &seedClusterList) const
{
    for (const Cluster *const pCluster : clusterList)
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

void CaloShowerGrowingAlgorithm::AssessAssociation(const Cluster *pSeed, const ClusterList &associatedClusterlist, ClusterList &showerClusterList) const
{
    CaloHitList caloHits;
    LArClusterHelper::GetAllHits(pSeed, caloHits);
    for (const Cluster *pCluster : associatedClusterlist)
        LArClusterHelper::GetAllHits(pCluster, caloHits);

    // Get the eigen vectors for this cluster
    CartesianVector origin(0.f, 0.f, 0.f);
    CartesianVector dir(0.f, 0.f, 0.f);
    // ATTN: PCA axis is a unit vector, so no need to correct for length in dot products
    this->GetProjectionAxis(caloHits, origin, dir);

    if (m_visualize)
    {
        CartesianVector source{origin - dir * 20.f}, sink{origin + dir * 20.f};
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "associated", BLUE));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &source, &sink, "axis", RED, 1, 1));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    // Project all hits onto the direction axis
    FloatVector projVector, energies;
    float min{std::numeric_limits<float>::max()}, max{-std::numeric_limits<float>::max()};
    for (const CaloHit *pCaloHit : caloHits)
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
    if (N < 1)
        return;
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
    std::cout << "Photon Forward: " << photonChi2Forward << " Backward: " << photonChi2Backward << std::endl;
    std::cout << "Photon R Forward: " << (photonChi2Forward / N) << " R Backward: " << (photonChi2Backward / N) << std::endl;
    std::cout << "Electron Forward: " << electronChi2Forward << " Backward: " << electronChi2Backward << std::endl;
    std::cout << "Electron R Forward: " << (electronChi2Forward / N) << " R Backward: " << (electronChi2Backward / N) << std::endl;

    (void)showerClusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CaloShowerGrowingAlgorithm::GetFigureOfMerit(const ClusterAssociationMap &clusterAssociationMap) const
{
    float figureOfMerit{0.f};
    (void)clusterAssociationMap;
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloShowerGrowingAlgorithm::GetProjectionAxis(const CaloHitList &caloHitList, CartesianVector &origin, CartesianVector &dir) const
{
    // Get the eigen vectors for this cluster
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(caloHitList, origin, eigenValues, eigenVectors);

    // Project the extremal hits of the cluster onto the principal axis
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsForSeed", m_minCaloHitsForSeed));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RadiationLength", m_radiationLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MoliereRadius", m_moliereRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
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

