/**
 *  @file   larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.cc
 *
 *  @brief  Implementation of the local asymmetry feature tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h"

using namespace pandora;

namespace lar_content
{

LocalAsymmetryFeatureTool::LocalAsymmetryFeatureTool() :
    m_maxAsymmetryDistance(5.f),
    m_minAsymmetryCosAngle(0.9962),
    m_maxAsymmetryNClusters(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LocalAsymmetryFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm, const Vertex * const pVertex,
    const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &slidingFitDataListMap, const VertexSelectionBaseAlgorithm::ClusterListMap &,
    const VertexSelectionBaseAlgorithm::KDTreeMap &, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float localAsymmetry(0.f);

    localAsymmetry += this->GetLocalAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U),
        slidingFitDataListMap.at(TPC_VIEW_U));

    localAsymmetry += this->GetLocalAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V),
        slidingFitDataListMap.at(TPC_VIEW_V));

    localAsymmetry += this->GetLocalAsymmetryForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W),
        slidingFitDataListMap.at(TPC_VIEW_W));

    featureVector.push_back(localAsymmetry);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LocalAsymmetryFeatureTool::GetLocalAsymmetryForView(const CartesianVector &vertexPosition2D,
    const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList) const
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
            useAsymmetry &= this->IncrementAsymmetryParameters(pCluster->GetElectromagneticEnergy(), clusterDirection, energyWeightedDirectionSum);
            useAsymmetry &= this->IncrementAsymmetryParameters(static_cast<float>(pCluster->GetNCaloHits()), clusterDirection, hitWeightedDirectionSum);
            asymmetryClusters.push_back(pCluster);
        }

        if (!useAsymmetry)
            return 1.f;
    }

    // Default: maximum asymmetry (i.e. not suppressed), zero for energy kick (i.e. not suppressed)
    if ((useEnergy && energyWeightedDirectionSum == CartesianVector(0.f, 0.f, 0.f)) || (!useEnergy && hitWeightedDirectionSum == CartesianVector(0.f, 0.f, 0.f)))
        return 1.f;

    const CartesianVector &localWeightedDirectionSum(useEnergy ? energyWeightedDirectionSum : hitWeightedDirectionSum);
    return useAsymmetry ? this->CalculateLocalAsymmetry(useEnergy, vertexPosition2D, asymmetryClusters, localWeightedDirectionSum) : 1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LocalAsymmetryFeatureTool::IncrementAsymmetryParameters(const float weight, const CartesianVector &clusterDirection,
    CartesianVector &localWeightedDirectionSum) const
{
    // If the new axis direction is at an angle of greater than 90 deg to the current axis direction, flip it 180 degs.
    CartesianVector newDirection(clusterDirection);

    if (localWeightedDirectionSum.GetMagnitudeSquared() > std::numeric_limits<float>::epsilon())
    {
        const float cosOpeningAngle(localWeightedDirectionSum.GetCosOpeningAngle(clusterDirection));

        if (std::fabs(cosOpeningAngle) > m_minAsymmetryCosAngle)
        {
            if (cosOpeningAngle < 0.f)
                newDirection *= -1.f;
        }

        else
            return false;
    }

    localWeightedDirectionSum += newDirection * weight;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LocalAsymmetryFeatureTool::CalculateLocalAsymmetry(const bool useEnergyMetrics, const CartesianVector &vertexPosition2D,
    const ClusterVector &asymmetryClusters, const CartesianVector &localWeightedDirectionSum) const
{
    if (asymmetryClusters.empty() || (asymmetryClusters.size() > m_maxAsymmetryNClusters))
        return 1.f;

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
    const float totHitEnergy(afterVtxHitEnergy + beforeVtxHitEnergy);
    const unsigned int totHitCount(beforeVtxHitCount + afterVtxHitCount);

    if (useEnergyMetrics && (totHitEnergy > std::numeric_limits<float>::epsilon()))
        return std::fabs((afterVtxHitEnergy - beforeVtxHitEnergy)) / totHitEnergy;

    if (0 == totHitCount)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return std::fabs((static_cast<float>(afterVtxHitCount) - static_cast<float>(beforeVtxHitCount))) / static_cast<float>(totHitCount);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LocalAsymmetryFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryDistance", m_maxAsymmetryDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAsymmetryCosAngle", m_minAsymmetryCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryNClusters", m_maxAsymmetryNClusters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
