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

ShowerAsymmetryFeatureTool::ShowerAsymmetryFeatureTool() : m_vertexClusterDistance(4.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerAsymmetryFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm,
    const Vertex *const pVertex, const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &,
    const VertexSelectionBaseAlgorithm::ClusterListMap &, const VertexSelectionBaseAlgorithm::KDTreeMap &,
    const VertexSelectionBaseAlgorithm::ShowerClusterListMap &showerClusterListMap, const float, float &)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float showerAsymmetry(0.f);

    showerAsymmetry += this->GetShowerAsymmetryForView(
        LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U), showerClusterListMap.at(TPC_VIEW_U));

    showerAsymmetry += this->GetShowerAsymmetryForView(
        LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V), showerClusterListMap.at(TPC_VIEW_V));

    showerAsymmetry += this->GetShowerAsymmetryForView(
        LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W), showerClusterListMap.at(TPC_VIEW_W));

    featureVector.push_back(showerAsymmetry);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShowerAsymmetryFeatureTool::GetShowerAsymmetryForView(
    const CartesianVector &vertexPosition2D, const VertexSelectionBaseAlgorithm::ShowerClusterList &showerClusterList) const
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

            const float projectedVtxPosition = vertexPosition2D.GetDotProduct(showerDirection);

            float beforeVtxEnergy(0.f), afterVtxEnergy(0.f);
            this->CalculateAsymmetryParameters(showerCluster, projectedVtxPosition, showerDirection, beforeVtxEnergy, afterVtxEnergy);

            if (beforeVtxEnergy + afterVtxEnergy > 0.f)
                showerAsymmetry = std::fabs(afterVtxEnergy - beforeVtxEnergy) / (afterVtxEnergy + beforeVtxEnergy);

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

void ShowerAsymmetryFeatureTool::CalculateAsymmetryParameters(const VertexSelectionBaseAlgorithm::ShowerCluster &showerCluster,
    const float projectedVtxPosition, const CartesianVector &showerDirection, float &beforeVtxEnergy, float &afterVtxEnergy) const
{
    beforeVtxEnergy = 0.f;
    afterVtxEnergy = 0.f;

    for (const Cluster *const pCluster : showerCluster.GetClusters())
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHitVector)
        {
            if (pCaloHit->GetPositionVector().GetDotProduct(showerDirection) < projectedVtxPosition)
                beforeVtxEnergy += pCaloHit->GetElectromagneticEnergy();

            else if (pCaloHit->GetPositionVector().GetDotProduct(showerDirection) > projectedVtxPosition)
                afterVtxEnergy += pCaloHit->GetElectromagneticEnergy();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerAsymmetryFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexClusterDistance", m_vertexClusterDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
