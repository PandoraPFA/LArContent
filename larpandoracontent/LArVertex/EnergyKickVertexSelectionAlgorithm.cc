/**
 *  @file   larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of the energy kick vertex selection algorithm class.
 * 
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EnergyKickVertexSelectionAlgorithm::EnergyKickVertexSelectionAlgorithm() :
    m_minClusterCaloHits(12),
    m_slidingFitWindow(100),
    m_rOffset(10.f),
    m_xOffset(0.06),
    m_epsilon(0.06),
    m_asymmetryConstant(3.f),
    m_maxAsymmetryDistance(5.f),
    m_minAsymmetryCosAngle(0.9962),
    m_maxAsymmetryNClusters(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &/*kdTreeU*/, HitKDTree2D &/*kdTreeV*/, HitKDTree2D &/*kdTreeW*/, VertexScoreList &vertexScoreList) const
{
    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    this->CalculateClusterSlidingFits(slidingFitDataListU, slidingFitDataListV, slidingFitDataListW);

    for (const Vertex *const pVertex : vertexVector)
    {
        const float vertexMinZ(std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate()));
        const float beamDeweightingScore(this->IsBeamModeOn() ? std::exp(-(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant()) : 1.f);

        float energyKick(0.f), energyAsymmetry(0.f);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U), energyKick, energyAsymmetry, slidingFitDataListU);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V), energyKick, energyAsymmetry, slidingFitDataListV);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W), energyKick, energyAsymmetry, slidingFitDataListW);

        const float energyKickScore(std::exp(-energyKick / m_epsilon));
        const float energyAsymmetryScore(std::exp(energyAsymmetry / m_asymmetryConstant));

        vertexScoreList.push_back(VertexScore(pVertex, beamDeweightingScore * energyKickScore * energyAsymmetryScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::CalculateClusterSlidingFits(SlidingFitDataList &slidingFitDataListU, SlidingFitDataList &slidingFitDataListV,
    SlidingFitDataList &slidingFitDataListW) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
             if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                 std::cout << "EnergyKickVertexSelectionAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        SlidingFitDataList &slidingFitDataList((TPC_VIEW_U == hitType) ? slidingFitDataListU : (TPC_VIEW_V == hitType) ? slidingFitDataListV : slidingFitDataListW);

        if (!slidingFitDataListW.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterVector sortedClusters(pClusterList->begin(), pClusterList->end());
        std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster * const pCluster : sortedClusters)
        {
            if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
                continue;
            
            // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
            const int slidingFitWindow(std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch * m_slidingFitWindow)));
            slidingFitDataList.emplace_back(pCluster, slidingFitWindow, slidingFitPitch);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyScoresForView(const CartesianVector &vertexPosition2D, float &energyKick, float &energyAsymmetry, 
    const SlidingFitDataList &slidingFitDataList) const
{
    unsigned int totHits(0);
    bool useEnergy(true), useAsymmetry(true);
    float totEnergy(0.f), totEnergyKick(0.f), totHitKick(0.f);
    CartesianVector energyWeightedDirectionSum(0.f, 0.f, 0.f), hitWeightedDirectionSum(0.f, 0.f, 0.f);
    ClusterVector asymmetryClusters;

    for (const SlidingFitData &slidingFitData : slidingFitDataList)
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

        if (useAsymmetry && (LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster) < m_maxAsymmetryDistance))
        {
            useAsymmetry &= this->IncrementEnergyAsymmetryParameters(pCluster->GetElectromagneticEnergy(), clusterDirection, energyWeightedDirectionSum);
            useAsymmetry &= this->IncrementEnergyAsymmetryParameters(static_cast<float>(pCluster->GetNCaloHits()), clusterDirection, hitWeightedDirectionSum);
            asymmetryClusters.push_back(pCluster);
        }
    }

    // Default: maximum asymmetry (i.e. not suppressed), zero for energy kick (i.e. not suppressed)
    if ((0 == totHits) || (useEnergy && (totEnergy < std::numeric_limits<float>::epsilon())))
    {
        energyAsymmetry += 1.f;
        return;
    }

    energyKick += useEnergy ? (totEnergyKick / totEnergy) : (totHitKick / static_cast<float>(totHits));
    const CartesianVector &localWeightedDirectionSum(useEnergy ? energyWeightedDirectionSum : hitWeightedDirectionSum);
    energyAsymmetry += useAsymmetry ? this->CalculateEnergyAsymmetry(useEnergy, vertexPosition2D, asymmetryClusters, localWeightedDirectionSum) : 1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyKickParameters(const Cluster *const pCluster, const CartesianVector &clusterDisplacement,
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

bool EnergyKickVertexSelectionAlgorithm::IncrementEnergyAsymmetryParameters(const float weight, const CartesianVector &clusterDirection,
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
        {
            return false;
        }
    }

    localWeightedDirectionSum += newDirection * weight;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::CalculateEnergyAsymmetry(const bool useEnergyMetrics, const CartesianVector &vertexPosition2D,
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
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);

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
    const float totHitEnergy(beforeVtxHitEnergy + beforeVtxHitEnergy);
    const unsigned int totHitCount(beforeVtxHitCount + afterVtxHitCount);
    
    if (useEnergyMetrics && (totHitEnergy > std::numeric_limits<float>::max()))
        return std::fabs((afterVtxHitEnergy - beforeVtxHitEnergy)) / totHitEnergy;

    if (0 == totHitCount)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return std::fabs((static_cast<float>(afterVtxHitCount) - static_cast<float>(beforeVtxHitCount))) / static_cast<float>(totHitCount);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EnergyKickVertexSelectionAlgorithm::SlidingFitData::SlidingFitData(const pandora::Cluster *const pCluster, const int slidingFitWindow, 
        const float slidingFitPitch) :
    m_minLayerDirection(0.f, 0.f, 0.f),
    m_maxLayerDirection(0.f, 0.f, 0.f),
    m_minLayerPosition(0.f, 0.f, 0.f),
    m_maxLayerPosition(0.f, 0.f, 0.f),
    m_pCluster(pCluster)
{
    const TwoDSlidingFitResult slidingFitResult(pCluster, slidingFitWindow, slidingFitPitch);
    m_minLayerDirection = slidingFitResult.GetGlobalMinLayerDirection();
    m_maxLayerDirection = slidingFitResult.GetGlobalMaxLayerDirection();
    m_minLayerPosition = slidingFitResult.GetGlobalMinLayerPosition();
    m_maxLayerPosition = slidingFitResult.GetGlobalMaxLayerPosition();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyKickVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ROffset", m_rOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "XOffset", m_xOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Epsilon", m_epsilon));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AsymmetryConstant", m_asymmetryConstant));

    if ((m_rOffset < std::numeric_limits<float>::epsilon()) || (m_epsilon < std::numeric_limits<float>::epsilon()) || (m_asymmetryConstant < std::numeric_limits<float>::epsilon()))
    {
        std::cout << "EnergyKickVertexSelection: Invalid parameter(s), ROffset " << m_rOffset << ", Epsilon " << m_epsilon << ", AsymmetryConstant " << m_asymmetryConstant << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryDistance", m_maxAsymmetryDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAsymmetryCosAngle", m_minAsymmetryCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryNClusters", m_maxAsymmetryNClusters));

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
