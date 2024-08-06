/**
 *  @file   larpandoracontent/LArVertex/AsymmetryFeatureBaseTool.cc
 *
 *  @brief  Implementation of the  asymmetry feature tool class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArVertex/AsymmetryFeatureBaseTool.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

AsymmetryFeatureBaseTool::AsymmetryFeatureBaseTool() :
    m_maxAsymmetryDistance(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AsymmetryFeatureBaseTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm,
    const Vertex *const pVertex, const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &slidingFitDataListMap,
    const VertexSelectionBaseAlgorithm::ClusterListMap &, const VertexSelectionBaseAlgorithm::KDTreeMap &,
    const VertexSelectionBaseAlgorithm::ShowerClusterListMap &showerClusterListMap, const float, float &)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float asymmetry(0.f);

    asymmetry += this->GetAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U),
        slidingFitDataListMap.at(TPC_VIEW_U),
        showerClusterListMap.empty() ? VertexSelectionBaseAlgorithm::ShowerClusterList() : showerClusterListMap.at(TPC_VIEW_U));

    asymmetry += this->GetAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V),
        slidingFitDataListMap.at(TPC_VIEW_V),
        showerClusterListMap.empty() ? VertexSelectionBaseAlgorithm::ShowerClusterList() : showerClusterListMap.at(TPC_VIEW_V));

    asymmetry += this->GetAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W),
        slidingFitDataListMap.at(TPC_VIEW_W),
        showerClusterListMap.empty() ? VertexSelectionBaseAlgorithm::ShowerClusterList() : showerClusterListMap.at(TPC_VIEW_W));

    featureVector.push_back(asymmetry);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AsymmetryFeatureBaseTool::IncrementAsymmetryParameters(
    const float weight, const CartesianVector &clusterDirection, CartesianVector &localWeightedDirectionSum) const
{
    // If the new axis direction is at an angle of greater than 90 deg to the current axis direction, flip it 180 degs.
    CartesianVector newDirection(clusterDirection);

    if (localWeightedDirectionSum.GetMagnitudeSquared() > std::numeric_limits<float>::epsilon())
    {
        if (localWeightedDirectionSum.GetCosOpeningAngle(clusterDirection) < 0.f)
            newDirection *= -1.f;
    }

    localWeightedDirectionSum += newDirection * weight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float AsymmetryFeatureBaseTool::CalculateAsymmetry(const bool useEnergyMetrics, const CartesianVector &vertexPosition2D,
    const ClusterVector &asymmetryClusters, const CartesianVector &localWeightedDirectionSum) const
{
    // Project every hit onto local event axis direction and record side of the projected vtx position on which it falls
    float beforeVtxHitEnergy(0.f), afterVtxHitEnergy(0.f);
    unsigned int beforeVtxHitCount(0), afterVtxHitCount(0);

    const CartesianVector localWeightedDirection(localWeightedDirectionSum.GetUnitVector());
    const float evtProjectedVtxPos(vertexPosition2D.GetDotProduct(localWeightedDirection));

    for (const Cluster *const pCluster : asymmetryClusters)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHitVector)
        {
            if (pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection) < evtProjectedVtxPos)
            {
                beforeVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
                ++beforeVtxHitCount;
            }
            else
            {
                afterVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
                ++afterVtxHitCount;
            }
        }
    }

    // Use energy metrics if possible, otherwise fall back on hit counting.
    const float totHitEnergy(beforeVtxHitEnergy + afterVtxHitEnergy);
    const unsigned int totHitCount(beforeVtxHitCount + afterVtxHitCount);

    if (useEnergyMetrics && (totHitEnergy > std::numeric_limits<float>::epsilon()))
        return std::fabs((afterVtxHitEnergy - beforeVtxHitEnergy)) / totHitEnergy;

    if (0 == totHitCount)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return std::fabs((static_cast<float>(afterVtxHitCount) - static_cast<float>(beforeVtxHitCount))) / static_cast<float>(totHitCount);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AsymmetryFeatureBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAsymmetryDistance", m_maxAsymmetryDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
