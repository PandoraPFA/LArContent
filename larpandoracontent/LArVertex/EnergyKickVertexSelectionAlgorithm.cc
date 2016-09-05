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
    m_slidingFitWindow(100),
    m_rOffset(10.f),
    m_xOffset(0.06),
    m_epsilon(0.06),
    m_asymmetryConstant(3.f),
    m_minNHits(12)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &/*kdTreeU*/,
    HitKDTree2D &/*kdTreeV*/, HitKDTree2D &/*kdTreeW*/, VertexScoreList &vertexScoreList) const
{   
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    this->CalculateClusterSlidingFits(TPC_VIEW_U, slidingFitPitch, slidingFitDataListU);
    this->CalculateClusterSlidingFits(TPC_VIEW_V, slidingFitPitch, slidingFitDataListV);
    this->CalculateClusterSlidingFits(TPC_VIEW_W, slidingFitPitch, slidingFitDataListW);
   
    for (const Vertex *const pVertex : vertexVector)
    {
        const float energyScore(this->GetEnergyScore(pVertex, beamConstants, slidingFitDataListU, slidingFitDataListV, slidingFitDataListW));
        vertexScoreList.push_back(VertexScore(pVertex, energyScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::CalculateClusterSlidingFits(const HitType hitType, const float slidingFitPitch, SlidingFitDataList &slidingFitDataList) const
{
    // Get the relevant 2D cluster list - TODO, replace this, not flexible...
    std::string listName;
    switch (hitType)
    {
        case TPC_VIEW_U:
            listName = "ClustersU";
            break;
        case TPC_VIEW_V:
            listName = "ClustersV";
            break;
        case TPC_VIEW_W:
            listName = "ClustersW";
            break;
        default:
            throw;
    }
    
    const ClusterList *pClusterList(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pClusterList));

    if (!pClusterList || pClusterList->empty())
    {
         if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
             std::cout << "EnergyKickVertexSelectionAlgorithm: unable to find current cluster list " << std::endl;
    }

    for (const Cluster * const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < this->m_minNHits)
            continue;
        
        // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
        const int slidingFitWindow = std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch*m_slidingFitWindow));
        slidingFitDataList.emplace_back(pCluster, slidingFitWindow, slidingFitPitch);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::GetEnergyScore(const Vertex *const pVertex, const BeamConstants &beamConstants, 
    const SlidingFitDataList &slidingFitDataListU, const SlidingFitDataList &slidingFitDataListV, const SlidingFitDataList &slidingFitDataListW) const
{
    // Calculate the beam deweighting score.
    const float vertexMinZ = std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate());
    const float beamDeweightingScore = std::exp(-(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant());
    
    // Calculate the energy kick and energy asymmetry scores.
    float energyKick(0.f);
    float energyAsymmetry(0.f);
    
    this->IncrementEnergyScoresForView(pVertex, energyKick, energyAsymmetry, TPC_VIEW_U, slidingFitDataListU);
    this->IncrementEnergyScoresForView(pVertex, energyKick, energyAsymmetry, TPC_VIEW_V, slidingFitDataListV);
    this->IncrementEnergyScoresForView(pVertex, energyKick, energyAsymmetry, TPC_VIEW_W, slidingFitDataListW);
                                  
    const float energyKickScore = std::exp(-energyKick / this->m_epsilon);                    
    const float energyAsymmetryScore = std::exp(energyAsymmetry / this->m_asymmetryConstant);
        
    return beamDeweightingScore * energyKickScore * energyAsymmetryScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyScoresForView(const Vertex *const pVertex, float &energyKick, float &energyAsymmetry, 
    const HitType hitType, const SlidingFitDataList &slidingFitDataList) const
{
    // Project the candidate vertex to the 2D view.
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    
    // Loop over the clusters in the cluster list and add up all the energy kicks. Revert to hit-based approach if any clusters are zero-energy.
    float totEnergyKick(0.f);
    float totEnergy(0.f);
    
    float totHitKick(0.f);
    unsigned int totHits(0U);
    
    // Find the local event axis: the energy-weighted average axis direction. If any clusters have zero energy, revert to hit-weighted. 
    CartesianVector localEvtAxisDirEnergy(0.f, 0.f, 0.f);
    CartesianVector localEvtAxisDirHits(0.f, 0.f, 0.f);
    ClusterList asymmetryConsideredClusters;
    
    bool asymmetryScoreIsViable = true;
    bool useEnergyMetrics = true;

    unsigned int clusterCount(0);

    for (const SlidingFitData &slidingFitData : slidingFitDataList)
    {        
        const CartesianVector pointToFitStartVector = slidingFitData.GetMinLayerPosition() - vertexPosition2D;
        const CartesianVector pointToFitEndVector   = slidingFitData.GetMaxLayerPosition() - vertexPosition2D;
        
        const float distanceToStart = pointToFitStartVector.GetMagnitude();
        const float distanceToEnd = pointToFitEndVector.GetMagnitude();
        const CartesianVector axisDirection = distanceToStart < distanceToEnd ? slidingFitData.GetMinLayerDirection() : slidingFitData.GetMaxLayerDirection();
        
        this->IncrementEnergyKickParameters(distanceToStart, distanceToEnd, pointToFitStartVector, pointToFitEndVector, axisDirection, 
            useEnergyMetrics, totEnergyKick, totEnergy, totHitKick, totHits, slidingFitData.GetCluster());
        
        if (asymmetryScoreIsViable)
        {
            this->IncrementEnergyAsymmetryParameters(vertexPosition2D, axisDirection, useEnergyMetrics, localEvtAxisDirEnergy, 
                localEvtAxisDirHits, slidingFitData.GetCluster(), asymmetryScoreIsViable, asymmetryConsideredClusters);
        }
        
        ++clusterCount;
    }
    
    if (clusterCount == 0)
    {
        energyAsymmetry += 1.f; // default value is maximum asymmetry (i.e. not suppressed) - default value for energy kick is 0.f (i.e. not suppressed)
        return;
    }
    
    energyKick += useEnergyMetrics ? (totEnergyKick / totEnergy) : (totHitKick / totHits);

    energyAsymmetry += this->CalculateEnergyAsymmetry(asymmetryConsideredClusters, vertexPosition2D, useEnergyMetrics, localEvtAxisDirEnergy, 
        localEvtAxisDirHits, asymmetryScoreIsViable);    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyKickParameters(const float distanceToStart, const float distanceToEnd, 
    const CartesianVector &pointToFitStartVector, const CartesianVector &pointToFitEndVector, const CartesianVector &axisDirection, bool &useEnergyMetrics, 
    float &totEnergyKick, float &totEnergy, float &totHitKick, unsigned int &totHits, const Cluster * const pCluster) const
{
    // Find the norm of the cross product between this vector and the cluster's direction (unit) vector. This is the impact parameter.
    const float impactParameter = distanceToStart < distanceToEnd ? pointToFitStartVector.GetCrossProduct(axisDirection).GetMagnitude() : 
        pointToFitEndVector.GetCrossProduct(axisDirection).GetMagnitude();
    
    const float closestEndDistance = std::min(distanceToStart, distanceToEnd);
    
    // If any cluster apparently has zero energy, use a hit-count based metric instead of energy.
    if (pCluster->GetElectromagneticEnergy() == 0.f)
        useEnergyMetrics = false;
    
    if (useEnergyMetrics)
    {
        totEnergyKick += pCluster->GetElectromagneticEnergy() * (impactParameter + m_xOffset) / (closestEndDistance + m_rOffset);
        totEnergy += pCluster->GetElectromagneticEnergy();
    }
    
    totHitKick += static_cast<float>(pCluster->GetNCaloHits()) * (impactParameter + m_xOffset) / (closestEndDistance + m_rOffset);
    totHits += pCluster->GetNCaloHits();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyAsymmetryParameters(const CartesianVector &vertexPosition2D, 
    const CartesianVector &axisDirection, bool &useEnergyMetrics, CartesianVector &localEvtAxisDirEnergy, CartesianVector &localEvtAxisDirHits,
    const Cluster *const pCluster, bool &isViable, ClusterList &asymmetryConsideredClusters) const
{
    if (LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster) > 5.0)
        return;
        
    asymmetryConsideredClusters.insert(pCluster);
    
    // The position is the direction of the fit at the end nearest to the vertex.
    CartesianVector axisDirectionEnergy = axisDirection;
    CartesianVector axisDirectionHits = axisDirection;
    
    // Switch to hit-based metric if apparent cluster energy is zero.
    if (pCluster->GetElectromagneticEnergy() == 0.f)
        useEnergyMetrics = false;
    
    // If the new axis direction is at an angle of greater than 90 deg to the current axis direction, flip it 180 degs.
    if (localEvtAxisDirEnergy == CartesianVector(0.f, 0.f, 0.f)) {}
    else if (useEnergyMetrics)
    {
        const float cosOpeningAngle = localEvtAxisDirEnergy.GetCosOpeningAngle(axisDirectionEnergy);
                
        if (std::abs(cosOpeningAngle) > 0.9962) // ~cos(5 deg)
        {
            if (cosOpeningAngle < 0.f)
                axisDirectionEnergy *= -1.f;
        }
        
        else
        {
            isViable = false;
            return;
        }
    }
    
    // Do the same for the hit-based vector.
    if (localEvtAxisDirHits == CartesianVector(0.f, 0.f, 0.f)) {}
    else
    {
        const float cosOpeningAngle = localEvtAxisDirHits.GetCosOpeningAngle(axisDirectionHits);
                
        if (std::abs(cosOpeningAngle) > 0.9962) // ~cos(5 deg)
        {
            if (cosOpeningAngle < 0.f)
                axisDirectionHits *= -1.f;
        }
        
        else
        {
            isViable = false;
            return;
        }
    }
    
    if (useEnergyMetrics)
        localEvtAxisDirEnergy += axisDirectionEnergy * pCluster->GetElectromagneticEnergy();  
  
    localEvtAxisDirHits += axisDirectionHits * pCluster->GetNCaloHits(); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::CalculateEnergyAsymmetry(const ClusterList &consideredClusters, const CartesianVector &vertexPosition2D, 
    const bool useEnergyMetrics, const CartesianVector &localEvtAxisDirEnergy, const CartesianVector &localEvtAxisDirHits, bool isViable) const
{
    if (consideredClusters.empty() || consideredClusters.size() > 2 || !isViable)
        return 1.f;
    
    const CartesianVector localEvtAxisDir = useEnergyMetrics ? localEvtAxisDirEnergy.GetUnitVector() : localEvtAxisDirHits.GetUnitVector();
    
    // Now project every hit of every considered cluster onto the local event axis direction and record on what side of the projected vtx 
    // position it falls.
    const float evtProjectedVtxPos = vertexPosition2D.GetDotProduct(localEvtAxisDir);
    float beforeVtxEnergy(0.f);
    float afterVtxEnergy(0.f);
    
    unsigned int beforeVtxCount(0U);
    unsigned int afterVtxCount(0U);
    
    for (const Cluster * const pCluster : consideredClusters)
    {
        for (const std::pair<const unsigned int, std::unordered_set<const pandora::CaloHit*>*> &orderedCaloHitList : pCluster->GetOrderedCaloHitList())
        {
            for (const CaloHit * const pCaloHit : *(orderedCaloHitList.second))
            {
                const float projectedPos = pCaloHit->GetPositionVector().GetDotProduct(localEvtAxisDir);
                
                if (projectedPos < evtProjectedVtxPos)
                {
                    beforeVtxEnergy += pCaloHit->GetElectromagneticEnergy();
                    ++beforeVtxCount;
                }
                    
                else if (projectedPos > evtProjectedVtxPos)
                {
                    afterVtxEnergy += pCaloHit->GetElectromagneticEnergy();
                    ++afterVtxCount;
                }
            }
        }
    }

    const float totCaloHitEnergy = beforeVtxEnergy + afterVtxEnergy;    
    const unsigned int totCount = beforeVtxCount + afterVtxCount;    
    
    if (useEnergyMetrics && totCaloHitEnergy > 0.f)
        return std::abs((afterVtxEnergy - beforeVtxEnergy)) / totCaloHitEnergy;
    
    // Otherwise fall back on hit counting.
    return totCount > 0U ? std::abs(static_cast<float>((afterVtxCount - beforeVtxCount))) / static_cast<float>(totCount) : 1.f;
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
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNHits", m_minNHits));

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
