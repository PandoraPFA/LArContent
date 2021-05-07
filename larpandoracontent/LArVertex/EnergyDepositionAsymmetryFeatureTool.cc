/**
 *  @file   larpandoracontent/LArVertex/EnergyDepositionAsymmetryFeatureTool.cc
 *
 *  @brief  Implementation of the energy deposition asymmetry feature tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArVertex/EnergyDepositionAsymmetryFeatureTool.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

EnergyDepositionAsymmetryFeatureTool::EnergyDepositionAsymmetryFeatureTool() :
    m_maxAsymmetryDistance(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyDepositionAsymmetryFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm, const Vertex * const pVertex,
    const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &slidingFitDataListMap, const VertexSelectionBaseAlgorithm::ClusterListMap &,
    const VertexSelectionBaseAlgorithm::KDTreeMap &, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float dEdxAsymmetry(0.f);

    dEdxAsymmetry += this->GetEnergyDepositionAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U),
        slidingFitDataListMap.at(TPC_VIEW_U));

    dEdxAsymmetry += this->GetEnergyDepositionAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V),
        slidingFitDataListMap.at(TPC_VIEW_V));

    dEdxAsymmetry += this->GetEnergyDepositionAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W),
        slidingFitDataListMap.at(TPC_VIEW_W));

    featureVector.push_back(dEdxAsymmetry);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyDepositionAsymmetryFeatureTool::GetEnergyDepositionAsymmetryForView(const CartesianVector &vertexPosition2D,
    const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList) const
{
    bool useEnergy(true);
    CartesianVector energyWeightedDirectionSum(0.f, 0.f, 0.f), hitWeightedDirectionSum(0.f, 0.f, 0.f);

    for (const VertexSelectionBaseAlgorithm::SlidingFitData &slidingFitData : slidingFitDataList)
    {
        const Cluster *const pCluster(slidingFitData.GetCluster());

        if (pCluster->GetElectromagneticEnergy() < std::numeric_limits<float>::epsilon())
            useEnergy = false;

        const CartesianVector vertexToMinLayer(slidingFitData.GetMinLayerPosition() - vertexPosition2D);
        const CartesianVector vertexToMaxLayer(slidingFitData.GetMaxLayerPosition() - vertexPosition2D);

        const bool minLayerClosest(vertexToMinLayer.GetMagnitudeSquared() < vertexToMaxLayer.GetMagnitudeSquared());
        const CartesianVector &clusterDirection((minLayerClosest) ? slidingFitData.GetMinLayerDirection() : slidingFitData.GetMaxLayerDirection());

        if ((LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster) < m_maxAsymmetryDistance))
        {
            this->IncrementAsymmetryParameters(pCluster->GetElectromagneticEnergy(), clusterDirection, energyWeightedDirectionSum);
            this->IncrementAsymmetryParameters(static_cast<float>(pCluster->GetNCaloHits()), clusterDirection, hitWeightedDirectionSum);
        }
    }

    const CartesianVector &localWeightedDirectionSum(useEnergy ? energyWeightedDirectionSum : hitWeightedDirectionSum);

    if (localWeightedDirectionSum.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon()){
      return 0.f;
    }

    return this->CalculateEnergyDepositionAsymmetry(useEnergy, vertexPosition2D, slidingFitDataList, localWeightedDirectionSum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyDepositionAsymmetryFeatureTool::IncrementAsymmetryParameters(const float weight, const CartesianVector &clusterDirection,
    CartesianVector &localWeightedDirectionSum) const
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

float EnergyDepositionAsymmetryFeatureTool::CalculateEnergyDepositionAsymmetry(const bool useEnergyMetrics, const CartesianVector &vertexPosition2D,
    const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList, const CartesianVector &localWeightedDirectionSum) const
{
    // Project every hit onto local event axis direction and record side of the projected vtx position on which it falls
    float beforeVtxHitEnergy(0.f), afterVtxHitEnergy(0.f);
    unsigned int beforeVtxHitCount(0), afterVtxHitCount(0);

    const CartesianVector localWeightedDirection(localWeightedDirectionSum.GetUnitVector());
    const float evtProjectedVtxPos(vertexPosition2D.GetDotProduct(localWeightedDirection));

    float minBeforeProjectedPos(std::numeric_limits<float>::max());
    float maxBeforeProjectedPos(-std::numeric_limits<float>::max());

    float minAfterProjectedPos(std::numeric_limits<float>::max());
    float maxAfterProjectedPos(-std::numeric_limits<float>::max());

    for (const VertexSelectionBaseAlgorithm::SlidingFitData &slidingFitData : slidingFitDataList)
    {
        const Cluster * const pCluster(slidingFitData.GetCluster());

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHitVector)
        {
            if (pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection) < evtProjectedVtxPos)
            {
	      minBeforeProjectedPos = std::min(minBeforeProjectedPos, pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection));
	      maxBeforeProjectedPos = std::max(maxBeforeProjectedPos, pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection));

	      beforeVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
	      ++beforeVtxHitCount;
            }

            else
            {
	      minAfterProjectedPos = std::min(minAfterProjectedPos, pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection));
	      maxAfterProjectedPos = std::max(maxAfterProjectedPos, pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection));
	      
	      afterVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
	      ++afterVtxHitCount;
            }
        }
    }

    // Use energy metrics if possible, otherwise fall back on hit counting.
    const unsigned int totHitCount(beforeVtxHitCount + afterVtxHitCount);

    const float beforeVtxEnergyDeposition(!std::fabs(maxBeforeProjectedPos - minBeforeProjectedPos) ? 0 : 
					  beforeVtxHitEnergy / std::fabs(maxBeforeProjectedPos - minBeforeProjectedPos));

    const float afterVtxEnergyDeposition(!std::fabs(maxAfterProjectedPos - minAfterProjectedPos) ? 0 : 
					 afterVtxHitEnergy / std::fabs(maxAfterProjectedPos - minAfterProjectedPos));

    const float totalEnergyDeposition(beforeVtxEnergyDeposition + afterVtxEnergyDeposition);

    if (useEnergyMetrics && totalEnergyDeposition)
        return std::fabs((afterVtxEnergyDeposition - beforeVtxEnergyDeposition)) / totalEnergyDeposition;

    if (0 == totHitCount)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float beforeVtxHitDeposition(!std::fabs(maxBeforeProjectedPos - minBeforeProjectedPos) ? 0 :
				       beforeVtxHitCount / std::fabs(maxBeforeProjectedPos - minBeforeProjectedPos));

    const float afterVtxHitDeposition(!std::fabs(maxAfterProjectedPos - minAfterProjectedPos) ? 0 : 
				      afterVtxHitCount / std::fabs(maxAfterProjectedPos - minAfterProjectedPos));

    const float totalHitDeposition(beforeVtxHitDeposition + afterVtxHitDeposition);

    if(totalHitDeposition)
      return std::fabs((afterVtxHitDeposition - beforeVtxHitDeposition)) / totalHitDeposition;

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyDepositionAsymmetryFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryDistance", m_maxAsymmetryDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
