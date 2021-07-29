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

using namespace pandora;

namespace lar_content
{

CaloShowerGrowingAlgorithm::CaloShowerGrowingAlgorithm() :
    m_minCaloHitsForSeed(5),
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
        // End temp stuff
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
            if (m_visualize)
            {
                if (isContained)
                {
                    ClusterList containedClusterList{pCluster};
                    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &containedClusterList, "contained", BLUE));
                }
                else
                {
                    ClusterList uncontainedClusterList{pCluster};
                    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &uncontainedClusterList, "uncontained", RED));
                }
            }
        }
        if (m_visualize)
        {
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
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
    const CartesianVector tAxis(-lAxis.GetZ() * 9.f, lAxis.GetY() * 9.f, lAxis.GetX() * 9.f);
    const float proj0{axis0.GetDotProduct(lAxis)};
    const float proj1{axis1.GetDotProduct(lAxis)};
    const float ext0{proj0 >= proj1 ? proj0 + 14.f : proj0 - 14.f};
    const float ext1{proj0 >= proj1 ? proj1 - 14.f : proj1 + 14.f};
    const CartesianVector a{lAxis * ext0 + centroid - tAxis};
    const CartesianVector b{lAxis * ext1 + centroid - tAxis};
    const CartesianVector c{lAxis * ext1 + centroid + tAxis};
    const CartesianVector d{lAxis * ext0 + centroid + tAxis};

    return Bounds(a, b, c, d);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CaloShowerGrowingAlgorithm::AreClustersAssociated(const Cluster *const pClusterSeed, const Cluster *const pCluster) const
{
    // Look at primary and transverse axes of the seed and determine if the cluster falls within a suitable bounding region
    (void)pClusterSeed; (void)pCluster;


    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CaloShowerGrowingAlgorithm::GetFigureOfMerit(const ClusterAssociationMap &clusterAssociationMap) const
{
    float figureOfMerit{0.f};
    (void)clusterAssociationMap;
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsForSeed", m_minCaloHitsForSeed));

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

