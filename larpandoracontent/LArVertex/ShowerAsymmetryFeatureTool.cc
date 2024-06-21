/**
 *  @file   larpandoracontent/LArVertex/ShowerAsymmetryFeatureTool.cc
 *
 *  @brief  Implementation of the shower asymmetry feature tool class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArVertex/ShowerAsymmetryFeatureTool.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

ShowerAsymmetryFeatureTool::ShowerAsymmetryFeatureTool() :
    AsymmetryFeatureBaseTool(),
    m_vertexClusterDistance(4.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShowerAsymmetryFeatureTool::GetAsymmetryForView(const CartesianVector &vertexPosition2D,
    const VertexSelectionBaseAlgorithm::SlidingFitDataList &, const VertexSelectionBaseAlgorithm::ShowerClusterList &showerClusterList) const
{
    float showerAsymmetry(1.f);

    for (const VertexSelectionBaseAlgorithm::ShowerCluster &showerCluster : showerClusterList)
    {
        if (this->ShouldUseShowerCluster(vertexPosition2D, showerCluster))
        {
            const TwoDSlidingFitResult &showerFit = showerCluster.GetFit();

            float rL(0.f), rT(0.f);
            showerFit.GetLocalPosition(vertexPosition2D, rL, rT);

            CartesianVector showerDirection(0.f, 0.f, 0.f);
            if (STATUS_CODE_SUCCESS != showerFit.GetGlobalFitDirection(rL, showerDirection))
                continue;

            ClusterVector asymmetryClusters;
            std::copy(showerCluster.GetClusters().begin(), showerCluster.GetClusters().end(), std::back_inserter(asymmetryClusters));

            showerAsymmetry = this->CalculateAsymmetry(true, vertexPosition2D, asymmetryClusters, showerDirection);

            break;
        }
    }

    return showerAsymmetry;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerAsymmetryFeatureTool::ShouldUseShowerCluster(
    const CartesianVector &vertexPosition, const VertexSelectionBaseAlgorithm::ShowerCluster &showerCluster) const
{
    for (const Cluster *const pCluster : showerCluster.GetClusters())
    {
        if (LArClusterHelper::GetClosestDistance(vertexPosition, pCluster) < m_vertexClusterDistance)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerAsymmetryFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexClusterDistance", m_vertexClusterDistance));

    return AsymmetryFeatureBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
