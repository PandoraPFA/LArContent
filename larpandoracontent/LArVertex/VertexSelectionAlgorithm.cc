/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex selection algorithm class.
 * 
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoracontent/LArVertex/VertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexSelectionAlgorithm::VertexSelectionAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_fastScoreCheck(true),
    m_fastScoreOnly(false),
    m_fullScore(false),
    m_energyScoreOnly(true),
    m_slidingFitWindow(100),
    m_rOffset(10.f),
    m_xOffset(0.06),
    m_epsilon(0.06),
    m_asymmetryConstant(3.f),
    m_minNHits(12),
    m_beamMode(false),
    m_nDecayLengthsInZSpan(2.f),
    m_kappa(0.42f),
    m_selectSingleVertex(true),
    m_maxTopScoreSelections(3),
    m_kernelEstimateSigma(0.048f),
    m_minFastScoreFraction(0.8f),
    m_fastHistogramNPhiBins(200),
    m_fastHistogramPhiMin(-1.1f * M_PI),
    m_fastHistogramPhiMax(+1.1f * M_PI),
    m_maxOnHitDisplacement(1.f),
    m_maxHitVertexDisplacement1D(100.f),
    m_minCandidateDisplacement(2.f),
    m_minCandidateScoreFraction(0.5f),
    m_enableFolding(true),
    m_useDetectorGaps(true),
    m_gapTolerance(0.f),
    m_isEmptyViewAcceptable(true),
    m_minVertexAcceptableViews(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::Run()
{
    const VertexList *pInputVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));

    if (!pInputVertexList || pInputVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexSelectionAlgorithm: unable to find current vertex list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    HitKDTree2D kdTreeU, kdTreeV, kdTreeW;
    this->InitializeKDTrees(kdTreeU, kdTreeV, kdTreeW);

    VertexList filteredVertexList;
    this->FilterVertexList(pInputVertexList, kdTreeU, kdTreeV, kdTreeW, filteredVertexList);

    if (filteredVertexList.empty())
        return STATUS_CODE_SUCCESS;

    BeamConstants beamConstants;
    this->GetBeamConstants(filteredVertexList, beamConstants);

    VertexScoreList vertexScoreList;
    this->GetVertexScoreList(filteredVertexList, beamConstants, kdTreeU, kdTreeV, kdTreeW, vertexScoreList);

    VertexList selectedVertexList;
    this->SelectTopScoreVertices(vertexScoreList, selectedVertexList);

    if (!selectedVertexList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, selectedVertexList));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    this->InitializeKDTree(m_inputCaloHitListNameU, kdTreeU);
    this->InitializeKDTree(m_inputCaloHitListNameV, kdTreeV);
    this->InitializeKDTree(m_inputCaloHitListNameW, kdTreeW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::InitializeKDTree(const std::string &caloHitListName, HitKDTree2D &kdTree) const
{
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, caloHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexSelectionAlgorithm: unable to find calo hit list " << caloHitListName << std::endl;

        return;
    }

    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(this, *pCaloHitList, hitKDNode2DList, true);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FilterVertexList(const VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
    HitKDTree2D &kdTreeW, VertexList &filteredVertexList) const
{
    for (const Vertex *const pVertex : *pInputVertexList)
    {
        unsigned int nAcceptableViews(0);

        if ((m_isEmptyViewAcceptable && kdTreeU.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeU) || this->IsVertexInGap(pVertex, TPC_VIEW_U))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeV.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeV) || this->IsVertexInGap(pVertex, TPC_VIEW_V))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeW.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeW) || this->IsVertexInGap(pVertex, TPC_VIEW_W))
            ++nAcceptableViews;

        if (nAcceptableViews >= m_minVertexAcceptableViews)
            (void) filteredVertexList.insert(pVertex);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::GetBeamConstants(const VertexList &vertexList, BeamConstants &beamConstants) const
{
    if (!m_beamMode)
        return;

    if (vertexList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    float minZCoordinate(std::numeric_limits<float>::max()), maxZCoordinate(-std::numeric_limits<float>::max());

    for (const Vertex *const pVertex : vertexList)
    {
        if (pVertex->GetPosition().GetZ() < minZCoordinate)
            minZCoordinate = pVertex->GetPosition().GetZ();

        if (pVertex->GetPosition().GetZ() > maxZCoordinate)
            maxZCoordinate = pVertex->GetPosition().GetZ();
    }

    const float zSpan(maxZCoordinate - minZCoordinate);
    const float decayConstant((zSpan < std::numeric_limits<float>::epsilon()) ? 0.f : (m_nDecayLengthsInZSpan / zSpan));
    beamConstants.SetConstants(minZCoordinate, decayConstant);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::GetVertexScoreList(const VertexList &vertexList, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
    HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    VertexVector vertexVector(vertexList.begin(), vertexList.end());
    std::sort(vertexVector.begin(), vertexVector.end(), SortByVertexZPosition);

    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    
    if (m_energyScoreOnly)
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        slidingFitDataListU = this->CalculateClusterSlidingFits(TPC_VIEW_U, slidingFitPitch);
        slidingFitDataListV = this->CalculateClusterSlidingFits(TPC_VIEW_V, slidingFitPitch);
        slidingFitDataListW = this->CalculateClusterSlidingFits(TPC_VIEW_W, slidingFitPitch);
    }

    float bestFastScore(0.f);
    
    for (const Vertex *const pVertex : vertexVector)
    {
        float finalScore(0.f);

        if (m_energyScoreOnly)
            finalScore = this->GetEnergyScore(pVertex, beamConstants, slidingFitDataListU, slidingFitDataListV, slidingFitDataListW);
        
        else
        {
            KernelEstimate kernelEstimateU(m_kernelEstimateSigma);
            KernelEstimate kernelEstimateV(m_kernelEstimateSigma);
            KernelEstimate kernelEstimateW(m_kernelEstimateSigma);

            this->FillKernelEstimate(pVertex, TPC_VIEW_U, kdTreeU, kernelEstimateU);
            this->FillKernelEstimate(pVertex, TPC_VIEW_V, kdTreeV, kernelEstimateV);
            this->FillKernelEstimate(pVertex, TPC_VIEW_W, kdTreeW, kernelEstimateW);

            const float vertexMinZ(!m_beamMode ? 0.f : std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate()));
            const float multiplier(!m_beamMode ? 1.f : std::exp(-(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant()));
            
            if (m_fastScoreCheck || m_fastScoreOnly)
            {
                const float fastScore(multiplier * this->GetFastScore(kernelEstimateU, kernelEstimateV, kernelEstimateW));

                if (m_fastScoreOnly)
                {
                    vertexScoreList.push_back(VertexScore(pVertex, fastScore));
                    continue;
                }

                if (fastScore < m_minFastScoreFraction * bestFastScore)
                    continue;

                if (fastScore > bestFastScore)
                    bestFastScore = fastScore;
            }

            finalScore = (multiplier * (m_fullScore ? this->GetFullScore(kernelEstimateU, kernelEstimateV, kernelEstimateW) :
                         this->GetMidwayScore(kernelEstimateU, kernelEstimateV, kernelEstimateW)));
        }
        
        vertexScoreList.push_back(VertexScore(pVertex, finalScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexSelectionAlgorithm::SlidingFitDataList VertexSelectionAlgorithm::CalculateClusterSlidingFits(const pandora::HitType hitType, 
                                                                                                   const float slidingFitPitch) const
{
    // Get the relevant 2D cluster list.
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
             std::cout << "VertexSelectionAlgorithm: unable to find current cluster list " << std::endl;
    }
    
    SlidingFitDataList slidingFitDataList;
    for (const Cluster * const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < this->m_minNHits)
            continue;
        
        // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
        const int slidingFitWindow = std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch*m_slidingFitWindow));
        slidingFitDataList.emplace_back(pCluster, slidingFitWindow, slidingFitPitch);
    }
        
    return slidingFitDataList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetEnergyScore(const pandora::Vertex *const pVertex, const BeamConstants &beamConstants, 
                                               const SlidingFitDataList &slidingFitDataListU, const SlidingFitDataList &slidingFitDataListV, 
                                               const SlidingFitDataList &slidingFitDataListW) const
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

void VertexSelectionAlgorithm::IncrementEnergyScoresForView(const pandora::Vertex *const pVertex, float &energyKick, float &energyAsymmetry, 
                                                            const pandora::HitType hitType, 
                                                            const SlidingFitDataList &slidingFitDataList) const
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
        const CartesianVector axisDirection = distanceToStart < distanceToEnd ? slidingFitData.GetMinLayerDirection() : 
                                                                                slidingFitData.GetMaxLayerDirection();
        
        this->IncrementEnergyKickParameters(distanceToStart, distanceToEnd, pointToFitStartVector, pointToFitEndVector, axisDirection, 
                                            useEnergyMetrics, totEnergyKick, totEnergy, totHitKick, totHits, slidingFitData.GetCluster());
        
        if (asymmetryScoreIsViable)
        {
            this->IncrementEnergyAsymmetryParameters(vertexPosition2D, axisDirection, useEnergyMetrics, localEvtAxisDirEnergy, 
                                                     localEvtAxisDirHits, slidingFitData.GetCluster(), asymmetryScoreIsViable, 
                                                     asymmetryConsideredClusters);
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

void VertexSelectionAlgorithm::IncrementEnergyKickParameters(const float distanceToStart, const float distanceToEnd, 
                                                             const CartesianVector &pointToFitStartVector, 
                                                             const CartesianVector &pointToFitEndVector, 
                                                             const CartesianVector &axisDirection, bool &useEnergyMetrics, 
                                                             float &totEnergyKick, float &totEnergy, float &totHitKick, 
                                                             unsigned int &totHits, const Cluster * const pCluster) const
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

void VertexSelectionAlgorithm::IncrementEnergyAsymmetryParameters(const CartesianVector &vertexPosition2D, 
                                                                  const CartesianVector &axisDirection, bool &useEnergyMetrics, 
                                                                  CartesianVector &localEvtAxisDirEnergy, 
                                                                  CartesianVector &localEvtAxisDirHits, const Cluster * const pCluster, 
                                                                  bool &isViable, ClusterList &asymmetryConsideredClusters) const
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

float VertexSelectionAlgorithm::CalculateEnergyAsymmetry(const ClusterList &consideredClusters, const CartesianVector &vertexPosition2D, 
                                                         const bool useEnergyMetrics, const CartesianVector &localEvtAxisDirEnergy, 
                                                         const CartesianVector &localEvtAxisDirHits, bool isViable) const
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

float VertexSelectionAlgorithm::GetFastScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    Histogram histogramU(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramV(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramW(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
        histogramU.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
        histogramV.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
        histogramW.Fill(contribution.first, contribution.second);

    // ATTN Need to renormalise histograms if ever want to directly compare fast and full scores
    float figureOfMerit(0.f);

    for (int xBin = 0; xBin < histogramU.GetNBinsX(); ++xBin)
    {
        const float binContentU(histogramU.GetBinContent(xBin));
        const float binContentV(histogramV.GetBinContent(xBin));
        const float binContentW(histogramW.GetBinContent(xBin));
        figureOfMerit += binContentU * binContentU + binContentV * binContentV + binContentW * binContentW;
    }

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetMidwayScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    Histogram histogramU(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramV(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramW(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
        histogramU.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
        histogramV.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
        histogramW.Fill(contribution.first, contribution.second);

    float figureOfMerit(0.f);

    for (int xBin = 0; xBin < histogramU.GetNBinsX(); ++xBin)
    {
        const float binCenter(histogramU.GetXLow() + (static_cast<float>(xBin) + 0.5f) * histogramU.GetXBinWidth());        
        figureOfMerit += histogramU.GetBinContent(xBin) * kernelEstimateU.Sample(binCenter);
        figureOfMerit += histogramV.GetBinContent(xBin) * kernelEstimateV.Sample(binCenter);
        figureOfMerit += histogramW.GetBinContent(xBin) * kernelEstimateW.Sample(binCenter);
    }

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFullScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    float figureOfMerit(0.f);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
        figureOfMerit += contribution.second * kernelEstimateU.Sample(contribution.first);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
        figureOfMerit += contribution.second * kernelEstimateV.Sample(contribution.first);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
        figureOfMerit += contribution.second * kernelEstimateW.Sample(contribution.first);

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::IsVertexOnHit(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxOnHitDisplacement, m_maxOnHitDisplacement);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    return (!found.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::IsVertexInGap(const Vertex *const pVertex, const HitType hitType) const
{
    if (!m_useDetectorGaps)
        return false;

    return LArGeometryHelper::IsInGap3D(this->GetPandora(), pVertex->GetPosition(), hitType, m_gapTolerance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FillKernelEstimate(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxHitVertexDisplacement1D, m_maxHitVertexDisplacement1D);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    for (const auto &hit : found)
    {
        const CartesianVector displacement(hit.data->GetPositionVector() - vertexPosition2D);
        const float magnitude(displacement.GetMagnitude());

        if (magnitude < std::numeric_limits<float>::epsilon())
            continue;

        float phi(this->atan2Fast(displacement.GetZ(), displacement.GetX()));
        float weight(1.f / (std::sqrt(magnitude + std::fabs(m_kappa))));

        if (m_enableFolding && (phi < 0.f))
        {
            phi += M_PI;
            weight *= -1.f;
        }

        kernelEstimate.AddContribution(phi, weight);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectTopScoreVertices(VertexScoreList &vertexScoreList, VertexList &selectedVertexList) const
{
    float bestScore(0.f);
    std::sort(vertexScoreList.begin(), vertexScoreList.end());

    for (const VertexScore &vertexScore : vertexScoreList)
    {
        if (selectedVertexList.size() >= m_maxTopScoreSelections)
            break;

        if (!selectedVertexList.empty() && !this->AcceptVertexLocation(vertexScore.GetVertex(), selectedVertexList))
            continue;

        if (!selectedVertexList.empty() && (vertexScore.GetScore() < m_minCandidateScoreFraction * bestScore))
            continue;

        selectedVertexList.insert(vertexScore.GetVertex());

        if (m_selectSingleVertex)
            return;

        if (vertexScore.GetScore() > bestScore)
            bestScore = vertexScore.GetScore();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexLocation(const Vertex *const pVertex, const VertexList &selectedVertexList) const
{
    const CartesianVector &position(pVertex->GetPosition());
    const float minCandidateDisplacementSquared(m_minCandidateDisplacement * m_minCandidateDisplacement);

    for (const Vertex *const pSelectedVertex : selectedVertexList)
    {
        if (pVertex == pSelectedVertex)
            return false;

        if ((position - pSelectedVertex->GetPosition()).GetMagnitudeSquared() < minCandidateDisplacementSquared)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::atan2Fast(const float y, const float x) const
{
    const float ONE_QTR_PI(0.25f * M_PI);
    const float THR_QTR_PI(0.75f * M_PI);

    const float abs_y(std::max(std::fabs(y), std::numeric_limits<float>::epsilon()));
    const float abs_x(std::fabs(x));

    const float r((x < 0.f) ? (x + abs_y) / (abs_y + abs_x) : (abs_x - abs_y) / (abs_x + abs_y));
    const float angle(((x < 0.f) ? THR_QTR_PI : ONE_QTR_PI) + (0.1963f * r * r - 0.9817f) * r);

    return ((y < 0.f) ? -angle : angle); // negate if in quad III or IV
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs)
{
    const CartesianVector deltaPosition(pRhs->GetPosition() - pLhs->GetPosition());

    if (std::fabs(deltaPosition.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    // ATTN No way to distinguish between vertices if still have a tie in y coordinate
    return (deltaPosition.GetY() > std::numeric_limits<float>::epsilon());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::KernelEstimate::Sample(const float x) const
{
    const ContributionList &contributionList(this->GetContributionList());
    ContributionList::const_iterator lowerIter(contributionList.lower_bound(x - 3.f * m_sigma));
    ContributionList::const_iterator upperIter(contributionList.upper_bound(x + 3.f * m_sigma));

    float sample(0.f);
    const float gaussConstant(1.f / std::sqrt(2.f * M_PI * m_sigma * m_sigma));

    for (ContributionList::const_iterator iter = lowerIter; iter != upperIter; ++iter)
    {
        const float deltaSigma((x - iter->first) / m_sigma);
        const float gaussian(gaussConstant * std::exp(-0.5f * deltaSigma * deltaSigma));
        sample += iter->second * gaussian;
    }

    return sample;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::KernelEstimate::AddContribution(const float x, const float weight)
{
    m_contributionList.insert(ContributionList::value_type(x, weight));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameU", m_inputCaloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameV", m_inputCaloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameW", m_inputCaloHitListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastScoreCheck", m_fastScoreCheck));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastScoreOnly", m_fastScoreOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FullScore", m_fullScore));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyScoreOnly", m_energyScoreOnly));
        
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamMode", m_beamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NDecayLengthsInZSpan", m_nDecayLengthsInZSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Kappa", m_kappa));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectSingleVertex", m_selectSingleVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreSelections", m_maxTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KernelEstimateSigma", m_kernelEstimateSigma));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinFastScoreFraction", m_minFastScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramNPhiBins", m_fastHistogramNPhiBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramPhiMin", m_fastHistogramPhiMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramPhiMax", m_fastHistogramPhiMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxOnHitDisplacement", m_maxOnHitDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitVertexDisplacement1D", m_maxHitVertexDisplacement1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateDisplacement", m_minCandidateDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateScoreFraction", m_minCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableFolding", m_enableFolding));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDetectorGaps", m_useDetectorGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GapTolerance", m_gapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IsEmptyViewAcceptable", m_isEmptyViewAcceptable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexAcceptableViews", m_minVertexAcceptableViews));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
