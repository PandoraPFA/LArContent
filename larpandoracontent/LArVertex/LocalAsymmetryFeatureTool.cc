/**
 *  @file   larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.cc
 *
 *  @brief  Implementation of the local asymmetry feature tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h"

using namespace pandora;

namespace lar_content
{

LocalAsymmetryFeatureTool::LocalAsymmetryFeatureTool() :
    AsymmetryFeatureBaseTool(), m_minAsymmetryCosAngle(0.9962), m_maxAsymmetryNClusters(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LocalAsymmetryFeatureTool::GetAsymmetryForView(const CartesianVector &vertexPosition2D,
    const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList, const VertexSelectionBaseAlgorithm::ShowerClusterList &) const
{
    bool useEnergy(true), useAsymmetry(true);
    CartesianVector energyWeightedDirectionSum(0.f, 0.f, 0.f), hitWeightedDirectionSum(0.f, 0.f, 0.f);
    ClusterVector asymmetryClusters;

    for (const VertexSelectionBaseAlgorithm::SlidingFitData &slidingFitData : slidingFitDataList)
    {
        const Cluster *const pCluster(slidingFitData.GetCluster());

        if (pCluster->GetElectromagneticEnergy() < std::numeric_limits<float>::epsilon())
            useEnergy = false;

        const CartesianVector vertexToMinLayer(slidingFitData.GetMinLayerPosition() - vertexPosition2D);
        const CartesianVector vertexToMaxLayer(slidingFitData.GetMaxLayerPosition() - vertexPosition2D);

        const bool minLayerClosest(vertexToMinLayer.GetMagnitudeSquared() < vertexToMaxLayer.GetMagnitudeSquared());
        const CartesianVector &clusterDirection((minLayerClosest) ? slidingFitData.GetMinLayerDirection() : slidingFitData.GetMaxLayerDirection());

        if (useAsymmetry && (LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster) < m_maxAsymmetryDistance))
        {
            useAsymmetry &= this->CheckAngle(energyWeightedDirectionSum, clusterDirection);
            this->IncrementAsymmetryParameters(pCluster->GetElectromagneticEnergy(), clusterDirection, energyWeightedDirectionSum);

            useAsymmetry &= this->CheckAngle(hitWeightedDirectionSum, clusterDirection);
            this->IncrementAsymmetryParameters(static_cast<float>(pCluster->GetNCaloHits()), clusterDirection, hitWeightedDirectionSum);

            asymmetryClusters.push_back(pCluster);
        }

        if (!useAsymmetry)
            return 1.f;
    }

    // Default: maximum asymmetry (i.e. not suppressed), zero for energy kick (i.e. not suppressed)
    if ((useEnergy && energyWeightedDirectionSum == CartesianVector(0.f, 0.f, 0.f)) ||
        (!useEnergy && hitWeightedDirectionSum == CartesianVector(0.f, 0.f, 0.f)))
        return 1.f;

    if (asymmetryClusters.empty() || (asymmetryClusters.size() > m_maxAsymmetryNClusters))
        return 1.f;

    const CartesianVector &localWeightedDirectionSum(useEnergy ? energyWeightedDirectionSum : hitWeightedDirectionSum);
    return this->CalculateAsymmetry(useEnergy, vertexPosition2D, asymmetryClusters, localWeightedDirectionSum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LocalAsymmetryFeatureTool::CheckAngle(const CartesianVector &weightedDirectionSum, const CartesianVector &clusterDirection) const
{
    if (!(weightedDirectionSum.GetMagnitudeSquared() > std::numeric_limits<float>::epsilon()))
        return true;

    const float cosOpeningAngle(weightedDirectionSum.GetCosOpeningAngle(clusterDirection));
    return std::fabs(cosOpeningAngle) > m_minAsymmetryCosAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LocalAsymmetryFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinAsymmetryCosAngle", m_minAsymmetryCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAsymmetryNClusters", m_maxAsymmetryNClusters));

    return AsymmetryFeatureBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
