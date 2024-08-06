/**
 *  @file   larpandoracontent/LArVertex/EnergyKickFeatureTool.cc
 *
 *  @brief  Implementation of the energy kick feature tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArVertex/EnergyKickFeatureTool.h"

using namespace pandora;

namespace lar_content
{

EnergyKickFeatureTool::EnergyKickFeatureTool() :
    m_rOffset(10.f),
    m_xOffset(0.06f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm,
    const Vertex *const pVertex, const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &slidingFitDataListMap,
    const VertexSelectionBaseAlgorithm::ClusterListMap &, const VertexSelectionBaseAlgorithm::KDTreeMap &,
    const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float energyKick(0.f);

    energyKick += this->GetEnergyKickForView(
        LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U), slidingFitDataListMap.at(TPC_VIEW_U));

    energyKick += this->GetEnergyKickForView(
        LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V), slidingFitDataListMap.at(TPC_VIEW_V));

    energyKick += this->GetEnergyKickForView(
        LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W), slidingFitDataListMap.at(TPC_VIEW_W));

    featureVector.push_back(energyKick);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickFeatureTool::GetEnergyKickForView(
    const CartesianVector &vertexPosition2D, const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList) const
{
    unsigned int totHits(0);
    bool useEnergy(true);
    float totEnergy(0.f), totEnergyKick(0.f), totHitKick(0.f);

    for (const VertexSelectionBaseAlgorithm::SlidingFitData &slidingFitData : slidingFitDataList)
    {
        const Cluster *const pCluster(slidingFitData.GetCluster());

        if (pCluster->GetElectromagneticEnergy() < std::numeric_limits<float>::epsilon())
            useEnergy = false;

        const CartesianVector vertexToMinLayer(slidingFitData.GetMinLayerPosition() - vertexPosition2D);
        const CartesianVector vertexToMaxLayer(slidingFitData.GetMaxLayerPosition() - vertexPosition2D);

        const bool minLayerClosest(vertexToMinLayer.GetMagnitudeSquared() < vertexToMaxLayer.GetMagnitudeSquared());
        const CartesianVector &clusterDisplacement((minLayerClosest) ? vertexToMinLayer : vertexToMaxLayer);
        const CartesianVector &clusterDirection((minLayerClosest) ? slidingFitData.GetMinLayerDirection() : slidingFitData.GetMaxLayerDirection());

        this->IncrementEnergyKickParameters(pCluster, clusterDisplacement, clusterDirection, totEnergyKick, totEnergy, totHitKick, totHits);
    }

    float energyKick(0.f);
    if (useEnergy && totEnergy > std::numeric_limits<float>::epsilon())
        energyKick = totEnergyKick / totEnergy;

    else if (!useEnergy && totHits > 0)
        energyKick = totHitKick / static_cast<float>(totHits);

    return energyKick;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickFeatureTool::IncrementEnergyKickParameters(const Cluster *const pCluster, const CartesianVector &clusterDisplacement,
    const CartesianVector &clusterDirection, float &totEnergyKick, float &totEnergy, float &totHitKick, unsigned int &totHits) const
{
    const float impactParameter(clusterDisplacement.GetCrossProduct(clusterDirection).GetMagnitude());
    const float displacement(clusterDisplacement.GetMagnitude());

    totEnergyKick += pCluster->GetElectromagneticEnergy() * (impactParameter + m_xOffset) / (displacement + m_rOffset);
    totEnergy += pCluster->GetElectromagneticEnergy();

    totHitKick += static_cast<float>(pCluster->GetNCaloHits()) * (impactParameter + m_xOffset) / (displacement + m_rOffset);
    totHits += pCluster->GetNCaloHits();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyKickFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ROffset", m_rOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "XOffset", m_xOffset));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
