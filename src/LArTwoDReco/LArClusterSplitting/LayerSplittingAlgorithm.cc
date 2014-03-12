/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/LayerSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the layer splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/LayerSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode LayerSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    unsigned int splitLayer(0);

    if (STATUS_CODE_SUCCESS == this->FindBestSplitLayer(pCluster, splitLayer))
        return this->SplitCluster(pCluster, splitLayer, firstHitList, secondHitList);

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LayerSplittingAlgorithm::FindBestSplitLayer(const Cluster *const pCluster, unsigned int &splitLayer) const
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.size() < m_minClusterLayers)
        return STATUS_CODE_NOT_FOUND;

    bool foundSplit(false);

    float bestCosTheta(1.f);
    CartesianVector bestPosition(0.f,0.f,0.f);

    for (unsigned int iLayer = pCluster->GetInnerPseudoLayer() + 4; iLayer + 4 <= pCluster->GetOuterPseudoLayer(); ++iLayer)
    {
        if (orderedCaloHitList.find(iLayer) == orderedCaloHitList.end())
            continue;

        unsigned int innerLayer((pCluster->GetInnerPseudoLayer() + m_layerWindow > iLayer) ? pCluster->GetInnerPseudoLayer() : iLayer - m_layerWindow);
        unsigned int outerLayer((iLayer + m_layerWindow > pCluster->GetOuterPseudoLayer()) ? pCluster->GetOuterPseudoLayer() : iLayer + m_layerWindow);

        for ( ; innerLayer >= pCluster->GetInnerPseudoLayer(); --innerLayer)
        {
            if (orderedCaloHitList.find(innerLayer) != orderedCaloHitList.end())
                break;
        }

        for ( ; outerLayer <= pCluster->GetOuterPseudoLayer(); ++outerLayer)
        {
            if (orderedCaloHitList.find(outerLayer) != orderedCaloHitList.end())
                break;
        }

        const CartesianVector splitPosition(pCluster->GetCentroid(iLayer));
        const CartesianVector innerPosition(pCluster->GetCentroid(innerLayer));
        const CartesianVector outerPosition(pCluster->GetCentroid(outerLayer));

        const CartesianVector r1(innerPosition - splitPosition);
        const CartesianVector r2(outerPosition - splitPosition);
        const CartesianVector p1(r1.GetUnitVector());
        const CartesianVector p2(r2.GetUnitVector());

        const float cosTheta(-p1.GetDotProduct(p2));
        const float rms1(this->CalculateRms(pCluster, innerLayer, iLayer));
        const float rms2(this->CalculateRms(pCluster, outerLayer, iLayer));
        const float rms(std::max(rms1, rms2));

        float rmsCut(std::numeric_limits<float>::max());

        if (cosTheta > 0.f)
        {
            rmsCut = m_maxScatterRms;

            if (cosTheta > m_maxScatterCosTheta)
            {
            rmsCut *= ((m_maxSlidingCosTheta > cosTheta) ? (m_maxSlidingCosTheta - cosTheta) /
                    (m_maxSlidingCosTheta - m_maxScatterCosTheta) : 0.f);
            }
        }

        if (rms < rmsCut && cosTheta < bestCosTheta)
        {
            bestCosTheta = cosTheta;
            bestPosition = splitPosition;

            splitLayer = iLayer;
            foundSplit = true;
        }
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

// --- BEGIN DISPLAY ---
// ClusterList tempList;
// tempList.insert((Cluster*)pCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList, "Cluster", GREEN));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&bestPosition, "SplitPosition", RED, 2.75));
// PANDORA_MONITORING_API(ViewEvent());
// --- END DISPLAY ---
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LayerSplittingAlgorithm::CalculateRms(const Cluster *const pCluster, const unsigned int &firstLayer, const unsigned int& secondLayer) const
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    const unsigned int innerLayer(std::min(firstLayer, secondLayer));
    const unsigned int outerLayer(std::max(firstLayer, secondLayer));

    const CartesianVector innerPosition(pCluster->GetCentroid(innerLayer));
    const CartesianVector outerPosition(pCluster->GetCentroid(outerLayer));
    const CartesianVector predictedDirection((outerPosition - innerPosition).GetUnitVector());

    float totalChi2(0.f);
    float totalLayers(0.f);

    for (unsigned int iLayer = innerLayer + 1; iLayer + 1 < outerLayer; ++iLayer)
    {
        if (orderedCaloHitList.find(iLayer) == orderedCaloHitList.end())
            continue;

        const CartesianVector hitPosition(pCluster->GetCentroid(iLayer));
        const CartesianVector predictedPosition(innerPosition + predictedDirection * predictedDirection.GetDotProduct(hitPosition - innerPosition));

        totalChi2 += (predictedPosition - hitPosition).GetMagnitudeSquared();
        totalLayers += 1.f;
    }

    if (totalLayers > 0.f)
        return std::sqrt(totalChi2/totalLayers);

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LayerSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, const unsigned int &splitLayer, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        const unsigned int thisLayer(iter->first);

        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;

            if (thisLayer < splitLayer)
            {
                firstHitList.insert(pCaloHit);
            }
            else
            {
                secondHitList.insert(pCaloHit);
            }
        }
    }

    if (firstHitList.empty() || secondHitList.empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LayerSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    m_layerWindow = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayerWindow", m_layerWindow));

    m_maxScatterRms = 0.35f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxScatterRms", m_maxScatterRms));

    m_maxScatterCosTheta = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxScatterCosTheta", m_maxScatterCosTheta));

    m_maxSlidingCosTheta = 0.866f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSlidingCosTheta", m_maxSlidingCosTheta));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
