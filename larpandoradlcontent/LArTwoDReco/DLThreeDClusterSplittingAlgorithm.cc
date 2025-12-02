/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/DLThreeDClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArTwoDReco/DLThreeDClusterSplittingAlgorithm.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

DLThreeDClusterSplittingAlgorithm::DLThreeDClusterSplittingAlgorithm() :
    m_trainingMode(false),
    m_fileName("ClusterSplittingTraining.root"),
    m_treeName("tree"),
    m_minClusterHits(50),
    m_slidingWindow(20),
    m_lBinSize(0.5f),
    m_kalmanDelta(1.f), 
    m_kalmanProcessVarCoeff(1.f), 
    m_kalmanMeasurementVarCoeff(1.f),
    m_searchRegion1D(20.f),
    m_windowLength(48),
    m_bkgThreshold(0.1f),
    m_isContaminatedThreshold(0.5f),
    m_isSplitThreshold(0.5f),
    m_minNuVertexSep(5.f),
    m_gapSearchRegion(2.f),
    m_threeViewMatchMaxXSpan(1.f),
    m_threeViewMatchMaxChi2(1.f),
    m_twoViewMatchMaxXSpan(0.5f),
    m_twoViewMatchSearchRegion(0.5f),
    m_contFraction(0.1f),
    m_mainFraction(0.9f),
    m_minTargetMCHits(5),
    m_endpointBuffer(2.f),
    m_collinearMaxSep(3.f),
    m_collinearMaxOpeningAngle(5.f),
    m_reducedMinTargetMCHits(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DLThreeDClusterSplittingAlgorithm::~DLThreeDClusterSplittingAlgorithm()
{
    if (m_trainingMode)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDClusterSplittingAlgorithm::Run()
{
    if (this->GetLists() != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    // For cluster features, we need a KDTree
    HitKDTree2D kdTreeU, kdTreeV, kdTreeW;
    HitKDNode2DList kdNode2DListU, kdNode2DListV, kdNode2DListW;
    KDTreeBox kdTreeBoxU(fill_and_bound_2d_kd_tree(*m_pCaloHitListU, kdNode2DListU));
    KDTreeBox kdTreeBoxV(fill_and_bound_2d_kd_tree(*m_pCaloHitListV, kdNode2DListV));
    KDTreeBox kdTreeBoxW(fill_and_bound_2d_kd_tree(*m_pCaloHitListW, kdNode2DListW));
    kdTreeU.build(kdNode2DListU, kdTreeBoxU);
    kdTreeV.build(kdNode2DListV, kdTreeBoxV);
    kdTreeW.build(kdNode2DListW, kdTreeBoxW);

    // Find split positions in all three views
    CartesianPointVector splitPosU, splitPosV, splitPosW;
    std::map<const Cluster*, IntVector> clusterToSplitIndexU, clusterToSplitIndexV, clusterToSplitIndexW;
    std::map<const Cluster*, TwoDSlidingFitResult> clusterFitsU, clusterFitsV, clusterFitsW;

    // Process clusters
    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const ClusterList *pClusterList(hitType == TPC_VIEW_U ? m_pClusterListU : hitType == TPC_VIEW_V ? 
            m_pClusterListV : m_pClusterListW);
        std::map<const Cluster*, TwoDSlidingFitResult> &clusterFits(hitType == TPC_VIEW_U ? clusterFitsU : 
            hitType == TPC_VIEW_V ? clusterFitsV : clusterFitsW);
        CartesianPointVector &splitPos(hitType == TPC_VIEW_U ? splitPosU : hitType == TPC_VIEW_V ? splitPosV : splitPosW);
        std::map<const Cluster*, IntVector> &clusterToSplitIndex(hitType == TPC_VIEW_U ? clusterToSplitIndexU :
            hitType == TPC_VIEW_V ? clusterToSplitIndexV : clusterToSplitIndexW);

        // Reproducability 
        ClusterVector internalClusterVector(pClusterList->begin(), pClusterList->end());
        std::sort(internalClusterVector.begin(), internalClusterVector.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pCluster : internalClusterVector)
        {
            CaloHitList clusterHits;
            LArClusterHelper::GetAllHits(pCluster, clusterHits);

            if (clusterHits.size() < m_minClusterHits)
                continue;

            // Perform sliding fit
            try
            {
                const TwoDSlidingFitResult clusterFit(pCluster, m_slidingWindow, 
                    LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));

                CartesianPointVector foundSplitPositions;
                this->ProcessCluster(pCluster, clusterFit, (hitType == TPC_VIEW_U ? kdTreeU : 
                    hitType == TPC_VIEW_V ? kdTreeV : kdTreeW), foundSplitPositions);

                // Store split points and the clusters to which they belong
                if (!m_trainingMode)
                {
                    int posIndex(splitPos.size());
                    for (unsigned int i = 0; i < foundSplitPositions.size(); ++i)
                    {
                        splitPos.push_back(foundSplitPositions.at(i));
                        clusterToSplitIndex[pCluster].push_back(posIndex);
                        clusterFits.insert(std::make_pair(pCluster, clusterFit));
                        ++posIndex;
                    }
                }
            }
            catch(...) { continue; }
        }
    }

    // Leave if there's nothing to do
    if (m_trainingMode || (splitPosU.empty() && splitPosV.empty() && splitPosW.empty()))
        return STATUS_CODE_SUCCESS;

    // Use 3D information to refine matches
    IntVector usedU, usedV, usedW;
    this->ThreeViewMatching(splitPosU, splitPosV, splitPosW, usedU, usedV, usedW);
    this->TwoViewMatching(splitPosU, splitPosV, TPC_VIEW_U, TPC_VIEW_V, kdTreeW, usedU, usedV);
    this->TwoViewMatching(splitPosU, splitPosW, TPC_VIEW_U, TPC_VIEW_W, kdTreeV, usedU, usedW);
    this->TwoViewMatching(splitPosV, splitPosW, TPC_VIEW_V, TPC_VIEW_W, kdTreeU, usedV, usedW);

    // Split the clusters
    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        // Adjust to view
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
            PandoraContentApi::ReplaceCurrentList<Cluster>(*this, hitType == TPC_VIEW_U ? m_clusterListNameU :
            hitType == TPC_VIEW_V ? m_clusterListNameV : m_clusterListNameW));
        const IntVector &used(hitType == TPC_VIEW_U ? usedU : hitType == TPC_VIEW_V ? usedV : usedW);
        CartesianPointVector &splitPos(hitType == TPC_VIEW_U ? splitPosU : hitType == TPC_VIEW_V ? splitPosV : splitPosW);
        std::map<const Cluster*, IntVector> &clusterToSplitIndex(hitType == TPC_VIEW_U ? clusterToSplitIndexU :
            hitType == TPC_VIEW_V ? clusterToSplitIndexV : clusterToSplitIndexW);
        std::map<const Cluster*, TwoDSlidingFitResult> &clusterFits(hitType == TPC_VIEW_U ? clusterFitsU : 
            hitType == TPC_VIEW_V ? clusterFitsV : clusterFitsW);

        for (const auto& [pCluster, splitIndices] : clusterToSplitIndex) 
        {
            CartesianPointVector foundSplitPositions;

            for (const int splitIndex : splitIndices)
            {
                if (std::find(used.begin(), used.end(), splitIndex) != used.end())
                    foundSplitPositions.push_back(splitPos.at(splitIndex));
            }

            // Divide calo hits
            if (clusterFits.find(pCluster) == clusterFits.end())
                throw;

            const TwoDSlidingFitResult &clusterFit(clusterFits.at(pCluster));
            std::vector<CaloHitList> splitClusterHits(this->DivideCaloHits(clusterFit, foundSplitPositions));

            if (splitClusterHits.size() < 2)
                continue;

            // Split clusters
            this->SplitCluster(pCluster, splitClusterHits);
        }        
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDClusterSplittingAlgorithm::GetLists()
{
    // Get 2D CaloHits - must find
    m_pCaloHitListU = nullptr; m_pCaloHitListV = nullptr; m_pCaloHitListW = nullptr; 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameU, m_pCaloHitListU))
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameV, m_pCaloHitListV))
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameW, m_pCaloHitListW))

    if ((!m_pCaloHitListU) || m_pCaloHitListU->empty() ||
        (!m_pCaloHitListV) || m_pCaloHitListV->empty() ||
        (!m_pCaloHitListW) || m_pCaloHitListW->empty())
    {
        return STATUS_CODE_NOT_FOUND;
    }

    // Get 2D Clusters - must find
    m_pClusterListU = nullptr; m_pClusterListV = nullptr; m_pClusterListW = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListNameU, m_pClusterListU))
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListNameV, m_pClusterListV))
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListNameW, m_pClusterListW))

    if ((!m_pClusterListU) || m_pClusterListU->empty() ||
        (!m_pClusterListV) || m_pClusterListV->empty() ||
        (!m_pClusterListW) || m_pClusterListW->empty())
    {
        return STATUS_CODE_NOT_FOUND;
    }

    // Get MCParticle - must find if training
    m_pMCParticleList = nullptr;
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, m_pMCParticleList));

        if ((!m_pMCParticleList) || m_pMCParticleList->empty())
            return STATUS_CODE_NOT_FOUND;
    }

    // Get NeutrinoVertex - must find
    m_pNuVertexList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_nuVertexListName, m_pNuVertexList));

    if (!m_pNuVertexList || (m_pNuVertexList->size() != 1))
        return STATUS_CODE_NOT_FOUND;

    // Get secondary vertices - okay if not found
    m_pSecVertexList = nullptr;
    PandoraContentApi::GetList(*this, m_secVertexListName, m_pSecVertexList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::ProcessCluster(const Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, 
    HitKDTree2D &kdTree, CartesianPointVector &splitPositions)
{
    // Find pathway through the cluster
    ClusterPath clusterPath;
    this->FindPath(pCluster, clusterFit, clusterPath);

    if (clusterPath.size() < m_windowLength)
        return;

    // Get features
    Features features;
    this->InitialiseFeatures(features);
    this->FillFeatures(clusterPath, clusterFit, features, kdTree);

    if (m_trainingMode)
    {
#ifdef MONITORING
        // mainHitListMap - to identify targets, contHitListMap - to use for topology
        MCParticleToHitListMap mainHitListMap, contHitListMap;
        this->FillHitListMaps(pCluster, mainHitListMap, contHitListMap);

        if (mainHitListMap.empty() || contHitListMap.empty())
            return;

        // Find contaminants
        MCContaminantVector mcContaminants;
        this->FindContaminants(mainHitListMap, contHitListMap, mcContaminants);

        // Now find split positions
        this->FindTrueSplitPositions(clusterFit, mcContaminants, splitPositions);

        // Fill tree
        this->FillTree(splitPositions, mainHitListMap, clusterFit, features);
#endif
    }
    else
    {
        // Get split indices
        IntVector splitIndices;
        this->GetPredictedSplitIndices(features, splitIndices);
        
        // Remove any that are too close to the nu vertex (network picks up on hit sharing)
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        this->FilterSplitIndices(clusterPath, hitType, features, splitIndices);
        
        if (splitIndices.empty())
            return;

        // Get split positions
        for (const int splitIndex : splitIndices)
        {
            CartesianVector position(clusterPath.at(splitIndex).m_pHit->GetPositionVector());
            splitPositions.push_back(position);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::FindPath(const Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, ClusterPath &clusterPath) const
{
    // Get cluster hits
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    // Fit l decomposition
    for (const CaloHit *const pCaloHit : clusterHits)
    {
        float thisHitL(0.f), thisHitT(0.f);
        clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisHitL, thisHitT);
        const int lBinIndex(std::floor(thisHitL / m_lBinSize));

        ClusterHit clusterHit(pCaloHit, lBinIndex, thisHitT);

        const auto iter(std::find(clusterPath.begin(), clusterPath.end(), clusterHit));

        if (iter == clusterPath.end())
        {
            clusterPath.emplace_back(clusterHit);
        }
        else if (thisHitT < iter->m_t)
        {
            *iter = clusterHit;
        }
    }

    // Order low->high l
    std::sort(clusterPath.begin(), clusterPath.end(), [](const ClusterHit &lhs, const ClusterHit &rhs) { return lhs.m_l < rhs.m_l; });
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::InitialiseFeatures(Features &features) const
{
    features.insert(std::make_pair("Longitudinal", Feature((0.0), 1.0, std::vector<float>())));
    features.insert(std::make_pair("Transverse", Feature((-0.01), 2.57, std::vector<float>())));
    features.insert(std::make_pair("Energy", Feature(0.49, 0.28, std::vector<float>())));
    features.insert(std::make_pair("HitWidth", Feature(0.56, 0.16, std::vector<float>())));
    features.insert(std::make_pair("Theta", Feature((-0.01), 0.11, std::vector<float>())));
    features.insert(std::make_pair("SecVertex", Feature(93.19, 103.39, std::vector<float>())));
    features.insert(std::make_pair("GapSep", Feature(0.0, 1.0, std::vector<float>())));
    features.insert(std::make_pair("EventHitSep", Feature(5.00, 5.83, std::vector<float>())));
    features.insert(std::make_pair("ClusterHitSep", Feature(0.54, 0.13, std::vector<float>())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::FillFeatures(const ClusterPath &clusterPath, const TwoDSlidingFitResult &clusterFit, Features &features,
    HitKDTree2D &kdTree) const
{
    // Get path hits, and their total energy
    CaloHitList clusterPathHits; float totalEnergy(0.f);
    for (const auto &entry : clusterPath)
    {
        const CaloHit *const pPathCaloHit(entry.m_pHit);
        totalEnergy += pPathCaloHit->GetElectromagneticEnergy();
        clusterPathHits.push_back(pPathCaloHit);
    }

    // Initialise KalmanFilter
    KalmanFilter2D kalmanFilter2D(this->InitialiseKalmanFilter(clusterPath));

    // To reduce computations and loops, get filtered secVertex list in view
    CartesianPointVector viewSecVtx;
    this->GetFilteredViewSecVertices(clusterFit, viewSecVtx);

    // Fill features
    bool processedFirst(false); float cumulativeEnergy(0.f);
    for (unsigned int i = 0; i < clusterPath.size(); ++i)
    {
        const ClusterHit &clusterHit(clusterPath.at(i));

        if (!processedFirst)
        {
            features.at("Theta").m_sequence.push_back(-4.f);
            processedFirst = true;
        }
        else
        {
            // Make next step
            kalmanFilter2D.Predict();
            // Update feature
            Eigen::VectorXd eigenXd(2);
            const CartesianVector thisPosition(clusterHit.m_l, 0.f, clusterHit.m_t);
            eigenXd << thisPosition.GetX(), thisPosition.GetZ();
            kalmanFilter2D.Update(eigenXd);
            // Get scatter angle
            const CartesianVector newDirection(kalmanFilter2D.GetState()(2), 0.f, kalmanFilter2D.GetState()(3));
            float openingAngleL(-4.f);
            try
            {
                const float openingAngleT = CartesianVector(0.f, 0.f, 1.f).GetOpeningAngle(newDirection);
                openingAngleL = CartesianVector(1.f, 0.f, 0.f).GetOpeningAngle(newDirection);
                openingAngleL *= (openingAngleT > (M_PI * 0.5f)) ? (-1.f) : 1.f;
            }
            catch (...) {};
            features.at("Theta").m_sequence.push_back(openingAngleL);
        }

        const CaloHit *const pPathCaloHit(clusterHit.m_pHit);
        cumulativeEnergy += pPathCaloHit->GetElectromagneticEnergy();
        features.at("Longitudinal").m_sequence.push_back(clusterHit.m_l * m_lBinSize);
        features.at("Transverse").m_sequence.push_back(clusterHit.m_t);
        features.at("Energy").m_sequence.push_back(cumulativeEnergy / totalEnergy);
        features.at("HitWidth").m_sequence.push_back(pPathCaloHit->GetCellSize1());
        features.at("SecVertex").m_sequence.push_back(
            this->GetDistanceToSecVertex(pPathCaloHit, viewSecVtx));
        features.at("EventHitSep").m_sequence.push_back(
            this->GetDistanceToEventHit(kdTree, pPathCaloHit, clusterPathHits));
        features.at("ClusterHitSep").m_sequence.push_back(i == (clusterPath.size() - 1) ? -1.f :
            this->GetDistanceToClusterHit(clusterHit, clusterPath.at(i + 1)));
        features.at("GapSep").m_sequence.push_back(i == 0 ? 1.f :
            this->IsInSameVolume(clusterPath.at(i -1).m_pHit, pPathCaloHit));
    }

    // Smooth and normalise
    for (auto &entry : features)
    {
        if ((entry.first == "GapSep") || (entry.first == "Longitudinal"))
            continue;

        entry.second.Smooth();
        entry.second.Normalise();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanFilter2D DLThreeDClusterSplittingAlgorithm::InitialiseKalmanFilter(const ClusterPath &clusterPath) const
{
    // Kalman Config
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view(clusterPath.begin()->m_pHit->GetHitType());
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float processVariance{m_kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{m_kalmanMeasurementVarCoeff * pitch * pitch};

    // Initialise Kalman fit
    Eigen::VectorXd init(2);
    const float seedL(clusterPath.begin()->m_l), seedT(clusterPath.begin()->m_t);
    init << seedL, seedT;
    KalmanFilter2D kalmanFilter2D(m_kalmanDelta, processVariance, measurementVariance, init);

    return kalmanFilter2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::GetFilteredViewSecVertices(const TwoDSlidingFitResult &clusterFit, CartesianPointVector &viewSecVtx) const
{
    if (!m_pSecVertexList)
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(clusterFit.GetCluster()));
    float clusterL1(-1.f), clusterL2(-1.f), clusterT1(-1.f), clusterT2(-1.f);
    clusterFit.GetLocalPosition(clusterFit.GetGlobalMinLayerPosition(), clusterL1, clusterT1);
    clusterFit.GetLocalPosition(clusterFit.GetGlobalMaxLayerPosition(), clusterL2, clusterT2);
    float clusterMinL(std::min(clusterL1, clusterL2)), clusterMaxL(std::max(clusterL1, clusterL2));

    for (const Vertex *const pSecVertex : *m_pSecVertexList)
    {
        const CartesianVector pos(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSecVertex->GetPosition(), hitType));

        float thisL(-1.f), thisT(-1.f);
        clusterFit.GetLocalPosition(pos, thisL, thisT);

        if ((thisL < clusterMinL) || (thisL > clusterMaxL))
            continue;

        viewSecVtx.push_back(pos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLThreeDClusterSplittingAlgorithm::GetDistanceToSecVertex(const CaloHit *const pCaloHit, const CartesianPointVector &viewSecVtx) const
{
    if (viewSecVtx.empty())
        return -1.f;

    float bestSepSq(std::numeric_limits<float>::max());
    for (const CartesianVector &secVtxPos : viewSecVtx)
        bestSepSq = std::min(bestSepSq, (secVtxPos - pCaloHit->GetPositionVector()).GetMagnitudeSquared());

    return std::sqrt(bestSepSq);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLThreeDClusterSplittingAlgorithm::GetDistanceToEventHit(HitKDTree2D &kdTree, const CaloHit *const pCaloHit, const CaloHitList &clusterPathHits) const
{
    // Collect close hits
    HitKDNode2DList foundHits;
    KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));
    kdTree.search(searchRegionHits, foundHits);

    // Filter
    bool found(false);
    float minDistSq(std::numeric_limits<float>::max());
    for (const auto &hit : foundHits)
    {
        const CaloHit *const pFoundHit(hit.data);

        if (std::find(clusterPathHits.begin(), clusterPathHits.end(), pFoundHit) != clusterPathHits.end())
            continue;

        found = true;
        minDistSq = std::min(minDistSq, (pCaloHit->GetPositionVector() - pFoundHit->GetPositionVector()).GetMagnitudeSquared());
    }

    return (found ? std::sqrt(minDistSq) : -1.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLThreeDClusterSplittingAlgorithm::GetDistanceToClusterHit(const ClusterHit &currentHit, 
    const ClusterHit &nextHit) const
{
    return (currentHit.m_pHit->GetPositionVector() - nextHit.m_pHit->GetPositionVector()).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLThreeDClusterSplittingAlgorithm::IsInSameVolume(const CaloHit *const pPrevHit, const CaloHit *const pCurrentHit) const
{
    if (!pPrevHit)
        return 1.f;

    const LArCaloHit *const pPrevLArHit(dynamic_cast<const LArCaloHit *>(pPrevHit));
    const LArCaloHit *const pCurrentLArHit(dynamic_cast<const LArCaloHit *>(pCurrentHit));

    if (!pPrevLArHit || !pCurrentLArHit)
        return 1.f;

    unsigned int prevTPCID(pPrevLArHit->GetLArTPCVolumeId());
    unsigned int currentTPCID(pCurrentLArHit->GetLArTPCVolumeId());

    if (prevTPCID != currentTPCID)
        return 0.f;

    unsigned int prevChildVolID(pPrevLArHit->GetDaughterVolumeId());
    unsigned int currentChildVolID(pCurrentLArHit->GetDaughterVolumeId());

    if (prevChildVolID != currentChildVolID)
        return 0.f;

    return 1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::GetPredictedSplitIndices(const Features &features, IntVector &splitIndices)
{
    // Find the window start indices
    IntVector windowStart;
    const int sequenceLength(features.at("Transverse").m_sequence.size());
    const int nWindows(std::floor(sequenceLength) / m_windowLength);

    for (int i = 0; i < nWindows; ++i)
        windowStart.push_back(m_windowLength * i);

    if (sequenceLength % m_windowLength != 0)
        windowStart.push_back(sequenceLength - m_windowLength);

    // Now get split indices and associated scores
    FloatVector splitScores;
    this->GetPredictedSplitIndices(features, windowStart, splitIndices, splitScores);

    // The model will often identify consecutive splitting positions around the truth, so identify one point
    this->FilterModelOutput(splitScores, splitIndices);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::GetPredictedSplitIndices(const Features &features, const IntVector &windowStart, 
    IntVector &splitIndices, FloatVector &splitScores)
{
    // Create network input
    LArDLHelper::TorchInput input; // (n_windows, sequence_length, n_features)
    LArDLHelper::InitialiseInput({static_cast<int>(windowStart.size()), m_windowLength, (int(features.size()) - 1)}, input); 

    for (unsigned int i = 0; i < windowStart.size(); ++i)
    {
        for (unsigned int j = 0; j < m_windowLength; ++j)
        {
            input[i][j][0] = features.at("Transverse").m_sequence.at(windowStart.at(i) + j);
            input[i][j][1] = features.at("Energy").m_sequence.at(windowStart.at(i) + j);
            input[i][j][2] = features.at("HitWidth").m_sequence.at(windowStart.at(i) + j);
            input[i][j][3] = features.at("Theta").m_sequence.at(windowStart.at(i) + j);
            input[i][j][4] = features.at("SecVertex").m_sequence.at(windowStart.at(i) + j);
            input[i][j][5] = features.at("GapSep").m_sequence.at(windowStart.at(i) + j);
            input[i][j][6] = features.at("EventHitSep").m_sequence.at(windowStart.at(i) + j);
            input[i][j][7] = features.at("ClusterHitSep").m_sequence.at(windowStart.at(i) + j);
        }
    }

    // Get model output
    LArDLHelper::TorchOutput windowOutput, splitPosOutput;
    LArDLHelper::Forward(m_windowModel, {input}, windowOutput);
    LArDLHelper::Forward(m_splitPosModel, {input}, splitPosOutput);
    torch::TensorAccessor<float, 2> windowOutputAccessor = windowOutput.accessor<float, 2>();
    torch::TensorAccessor<float, 3> splitPosOutputAccessor = splitPosOutput.accessor<float, 3>();

    for (unsigned int i = 0; i < windowStart.size(); ++i)
    {
        // Apply softmax
        float bkgScore(exp(windowOutputAccessor[i][0])), sigScore(exp(windowOutputAccessor[i][1])), shrScore(exp(windowOutputAccessor[i][2]));
        float bkgProb = (bkgScore) / (bkgScore + sigScore + shrScore);
        float sigProb = (sigScore) / (bkgScore + sigScore + shrScore);

        if (bkgProb > m_bkgThreshold)
            continue;

        // Is contaminated?
        if (sigProb > m_isContaminatedThreshold)
        {
            // Is there a split point?
            for (unsigned int j = 0; j < m_windowLength; ++j)
            {
                const int sequenceIndex(windowStart.at(i) + j);
                const int nextWindowStart((i+1) == windowStart.size() ? std::numeric_limits<int>::max() : windowStart.at(i+1));

                // make sure we only use the last window for the overlap region...
                if (sequenceIndex >= nextWindowStart)
                    continue;

                // Apply sigmoid
                float splitProb(1.f / (1.f + exp(-splitPosOutputAccessor[i][j][0])));

                if (splitProb > m_isSplitThreshold)
                {
                    splitIndices.push_back(sequenceIndex);
                    splitScores.push_back(static_cast<float>(splitProb));
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::FilterModelOutput(const FloatVector &splitScores, IntVector &splitIndices) const
{
    if (splitIndices.empty())
        return;

    IntVector splitIndices_temp(splitIndices);
    splitIndices.clear();

    int bestIndex = splitIndices_temp.front();
    float bestScore = splitScores.front();
    int previousIndex = splitIndices_temp.front();

    for (unsigned int i = 1; i < splitIndices_temp.size(); ++i)
    {
        int thisIndex = splitIndices_temp.at(i);
        float thisScore = splitScores.at(i);

        if (thisIndex == previousIndex + 1)
        {
            if (thisScore > bestScore)
            {
                bestScore = thisScore;
                bestIndex = thisIndex;
            }
        }
        else
        {
            splitIndices.push_back(bestIndex);

            bestIndex = thisIndex;
            bestScore = thisScore;
        }

        previousIndex = thisIndex;
    }

    // Push the last run
    splitIndices.push_back(bestIndex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::FilterSplitIndices(const ClusterPath &clusterPath, const HitType hitType, const Features &features, 
    IntVector &splitIndices) const
{
    IntVector filtered;

    // Reject if near nu vertex
    const CartesianVector &nuVertex3D(m_pNuVertexList->front()->GetPosition());
    const CartesianVector nuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    for (const int splitIndex : splitIndices)
    {
        if (splitIndex >= static_cast<int>(clusterPath.size()))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const CartesianVector &position(clusterPath.at(splitIndex).m_pHit->GetPositionVector());
        
        if ((nuVertexPosition - position).GetMagnitude() > m_minNuVertexSep)
            filtered.push_back(splitIndex);
    }

    // Reject if near gap 
    if (filtered.size() != splitIndices.size())
    {
        const FloatVector &gapSep(features.at("GapSep").m_sequence);

        for (const int splitIndex : splitIndices)
        {
            if (std::find(filtered.begin(), filtered.end(), splitIndex) != filtered.end())
                continue;

            bool nearGap(false);
            for (int offset = -m_gapSearchRegion; offset <= m_gapSearchRegion; ++offset)
            {
                int newIndex(splitIndex + offset);

                // Gap indicated by 0.f
                if ((newIndex >= 0) && (newIndex < int(gapSep.size())) && (gapSep.at(newIndex) < 0.5f))
                {
                    nearGap = true;
                    break;
                }
            }

            if (!nearGap)
                filtered.push_back(splitIndex);
        }
    }

    splitIndices.swap(filtered);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::ThreeViewMatching(const CartesianPointVector &splitPosU, const CartesianPointVector &splitPosV, 
    const CartesianPointVector &splitPosW, IntVector &usedU, IntVector &usedV, IntVector &usedW)
{
    for (unsigned int iU = 0; iU < splitPosU.size(); ++iU)
    {
        if (std::find(usedU.begin(), usedU.end(), iU) != usedU.end())
            continue;

        for (unsigned int iV = 0; iV < splitPosV.size(); ++iV)
        {
            if (std::find(usedV.begin(), usedV.end(), iV) != usedV.end())
                continue;
                    
            for (unsigned int iW = 0; iW < splitPosW.size(); ++iW)
            {
                if (std::find(usedW.begin(), usedW.end(), iW) != usedW.end())
                    continue;
                        
                const float minX(std::min(splitPosU.at(iU).GetX(), std::min(splitPosV.at(iV).GetX(), splitPosW.at(iW).GetX())));
                const float maxX(std::max(splitPosU.at(iU).GetX(), std::max(splitPosV.at(iV).GetX(), splitPosW.at(iW).GetX())));

                if ((maxX - minX) > m_threeViewMatchMaxXSpan)
                    continue;

                const float predW(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, 
                    splitPosU.at(iU).GetZ(), splitPosV.at(iV).GetZ()));
                const float predV(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W,
                    splitPosU.at(iU).GetZ(), splitPosW.at(iW).GetZ()));
                const float predU(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W,
                    splitPosV.at(iV).GetZ(), splitPosW.at(iW).GetZ()));                   
                const float pseudoChi2((std::fabs(predW - splitPosW.at(iW).GetZ()) + 
                                        std::fabs(predV - splitPosV.at(iV).GetZ()) +
                                        std::fabs(predU - splitPosU.at(iU).GetZ())) / 3.f);
                        
                if (pseudoChi2 > m_threeViewMatchMaxChi2)
                    continue;
                 
                usedU.push_back(iU); usedV.push_back(iV); usedW.push_back(iW);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::TwoViewMatching(const CartesianPointVector &splitPos1, const CartesianPointVector &splitPos2, 
    const HitType hitType1, const HitType hitType2, HitKDTree2D &kdTree, IntVector &used1, IntVector &used2)
{
    for (unsigned int i1 = 0; i1 < splitPos1.size(); ++i1)
    {
        if (std::find(used1.begin(), used1.end(), i1) != used1.end())
            continue;

        for (unsigned int i2 = 0; i2 < splitPos2.size(); ++i2)
        {
            if (std::find(used2.begin(), used2.end(), i2) != used2.end())
                continue;
                    
            if (std::fabs(splitPos1.at(i1).GetX() - splitPos2.at(i2).GetX()) > m_twoViewMatchMaxXSpan)
                continue;

            const float pred3(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, splitPos1.at(i1).GetZ(), splitPos2.at(i2).GetZ()));
            const CartesianVector predPosition3((splitPos1.at(i1).GetX() + splitPos2.at(i2).GetX()) * 0.5f, 0.f, pred3);

            // Collect close hits
            HitKDNode2DList foundHits;
            KDTreeBox searchRegionHits(build_2d_kd_search_region(predPosition3, m_twoViewMatchSearchRegion, m_twoViewMatchSearchRegion));
            kdTree.search(searchRegionHits, foundHits);

            if (foundHits.empty())
                continue;

            used1.push_back(i1); used2.push_back(i2);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<CaloHitList> DLThreeDClusterSplittingAlgorithm::DivideCaloHits(const TwoDSlidingFitResult &clusterFit,
    const CartesianPointVector &splitPositions) const
{
    // Get cluster hits
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(clusterFit.GetCluster(), clusterHits);

    // Get l-coord of split index
    FloatVector lSplit;
    for (const CartesianVector &splitPosition : splitPositions)
    {
        float thisL(0.f), thisT(0.f);
        clusterFit.GetLocalPosition(splitPosition, thisL, thisT);
        lSplit.push_back(thisL);
     }
    std::sort(lSplit.begin(), lSplit.end());

    lSplit.insert(lSplit.begin(), std::numeric_limits<float>::lowest());
    lSplit.insert(lSplit.end(), std::numeric_limits<float>::max());
    const int nSplitPositions(lSplit.size() - 2);
    std::vector<CaloHitList> splitClusterHits((nSplitPositions + 1), CaloHitList());

    for (const CaloHit *const pCaloHit : clusterHits)
    {
        float thisL(0.f), thisT(0.f);
        clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

        for (int i = 0; i <= nSplitPositions; ++i)
        {
            if ((thisL > lSplit.at(i)) && (thisL < lSplit.at(i + 1)))
            {
                splitClusterHits[i].push_back(pCaloHit);
                break;
            }
        }
    }

    // Remove empty lists
    splitClusterHits.erase(std::remove_if(splitClusterHits.begin(), splitClusterHits.end(),
        [](const CaloHitList& list){ return list.empty(); }), splitClusterHits.end());

    return splitClusterHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, const std::vector<CaloHitList> &splitClusterHits) const
{
    // Begin cluster fragmentation operations
    const ClusterList clusterList(1, pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName, clusterListToSaveName));

    for (const CaloHitList &caloHitList : splitClusterHits)
    {
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList = caloHitList;

        const Cluster *pNewCluster(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
    }

    // End cluster fragmentation operations
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Training Functions
//------------------------------------------------------------------------------------------------------------------------------------------
#ifdef MONITORING

void DLThreeDClusterSplittingAlgorithm::FillHitListMaps(const Cluster *const pCluster, MCParticleToHitListMap &mainHitListMap, 
    MCParticleToHitListMap &contHitListMap) 
{
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    for (const CaloHit *const pCaloHit : clusterHits)
    {
        // Identify contributing MCParticles
        const MCParticleWeightMap hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        for (const auto &entry : hitMCParticleWeightMap)
        {
            if (entry.second > m_contFraction)
                contHitListMap[entry.first].push_back(pCaloHit);
        }
        // Identify dominant MCParticles
        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            mainHitListMap[pMCParticle].push_back(pCaloHit);
        }
        catch (...) { continue; }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::FindContaminants(const MCParticleToHitListMap &mainHitListMap, const MCParticleToHitListMap &contHitListMap,
    MCContaminantVector &mcContaminants)
{
    // Sort
    MCParticleVector mcParticleVector;
    for (const auto &entry : mainHitListMap) mcParticleVector.push_back(entry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

    // Find contaminants
    for (const MCParticle *const pMCContaminant : mcParticleVector)
    {
        // Get hits MCParticle contributed to
        if (contHitListMap.find(pMCContaminant) == contHitListMap.end())
            continue;

        const CaloHitList &mainHitList(mainHitListMap.at(pMCContaminant));
        const CaloHitList &contHitList(contHitListMap.at(pMCContaminant));

        // Can we clearly see the particle in the cluster?
        unsigned int hitCount(0);
        for (const CaloHit *const pCaloHit : mainHitList)
        {
            const MCParticleWeightMap hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
            for (const auto &[pMCParticle, weight] : hitMCParticleWeightMap)
            {
                if (weight > m_mainFraction)
                {
                    ++hitCount;
                    break;
                }
            }
        }

        if (hitCount < m_minTargetMCHits)
            continue;

        this->BuildMCContaminant(pMCContaminant, mainHitList, contHitList, mcContaminants);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::BuildMCContaminant(const MCParticle *const pMCContaminant, const CaloHitList &mainHitList, 
    const CaloHitList &contHitList, MCContaminantVector &mcContaminants)
{
    const HitType hitType(mainHitList.front()->GetHitType());
    const CartesianVector trueStart(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetVertex(), hitType));
    // This isn't great for photons, but I think that's okay, they're not the target here
    const CartesianVector trueEnd(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetEndpoint(), hitType));

    CartesianVector contStartPosition(0.f,0.f,0.f), contEndPosition(0.f,0.f,0.f);
    CartesianVector mainStartPosition(0.f,0.f,0.f), mainEndPosition(0.f,0.f,0.f);
    CartesianVector mainStartDirection(0.f,0.f,0.f), mainEndDirection(0.f,0.f,0.f);
    CartesianPointVector fitPositions;

    // Identify start/endpoints
    float contStartSepSq(std::numeric_limits<float>::max()), contEndSepSq(std::numeric_limits<float>::max());
    float mainStartSepSq(std::numeric_limits<float>::max()), mainEndSepSq(std::numeric_limits<float>::max());
        
    for (const CaloHit *const pContHit : contHitList)
    {
        float this_startSepSq((pContHit->GetPositionVector() - trueStart).GetMagnitudeSquared());
        float this_endSepSq((pContHit->GetPositionVector() - trueEnd).GetMagnitudeSquared());

        if (this_startSepSq < contStartSepSq)
        {
            contStartSepSq = this_startSepSq;
            contStartPosition = pContHit->GetPositionVector();
        }
        if (this_endSepSq < contEndSepSq)
        { 
            contEndSepSq = this_endSepSq;
            contEndPosition = pContHit->GetPositionVector();
        }

        // Now find start/end positions used for splitting
        if (std::find(mainHitList.begin(), mainHitList.end(), pContHit) != mainHitList.end())
        {
            if (this_startSepSq < mainStartSepSq)
            {
                mainStartSepSq = this_startSepSq;
                mainStartPosition = pContHit->GetPositionVector();
            }
            if (this_endSepSq < mainEndSepSq)
            { 
                mainEndSepSq = this_endSepSq;
                mainEndPosition = pContHit->GetPositionVector();
            }
        }

        // Add to fit
        fitPositions.push_back(pContHit->GetPositionVector());
    }

    // Identify start/end directions
    try
    {
        const TwoDSlidingFitResult clusterFit(&fitPositions, m_slidingWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
        float startL(-1.f), startT(-1.f), endL(-1.f), endT(-1.f);
        clusterFit.GetLocalPosition(mainStartPosition, startL, startT);
        clusterFit.GetLocalPosition(mainEndPosition, endL, endT);
        clusterFit.GetGlobalFitDirection(startL, mainStartDirection);
        clusterFit.GetGlobalFitDirection(endL, mainEndDirection);
        mcContaminants.push_back(MCContaminant(pMCContaminant, contStartPosition, contEndPosition, 
            mainStartPosition, mainEndPosition, mainStartDirection, mainEndDirection));
    }
    catch (...)
    {
        if (pMCContaminant->GetMomentum().GetMagnitude() < std::numeric_limits<float>::epsilon())
            return;

        mainStartDirection = LArGeometryHelper::ProjectDirection(this->GetPandora(), pMCContaminant->GetMomentum().GetUnitVector(), hitType);
        mainEndDirection = mainStartDirection;
        mcContaminants.push_back(MCContaminant(pMCContaminant, contStartPosition, contEndPosition, 
            mainStartPosition, mainEndPosition, mainStartDirection, mainEndDirection));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::FindTrueSplitPositions(const TwoDSlidingFitResult &clusterFit, const MCContaminantVector &mcContaminants,
    CartesianPointVector &trueSplitPositions)
{
    if (mcContaminants.size() < 2)
        return;
 
    const CartesianVector clusterMin(clusterFit.GetGlobalMinLayerPosition());
    const CartesianVector clusterMax(clusterFit.GetGlobalMaxLayerPosition());

    for (const MCContaminant &thisContaminant : mcContaminants)
    {
        // Identify rejected indices
        IntVector rejected;

        CartesianPointVector splitPositions({thisContaminant.m_startPosition, thisContaminant.m_endPosition});
        CartesianPointVector splitDirections({thisContaminant.m_startDirection, thisContaminant.m_endDirection});
        FloatVector splitL({-1.f, -1.f}), splitT({-1.f, -1.f});
        for (int i = 0; i < 2; ++i) { clusterFit.GetLocalPosition(splitPositions.at(i), splitL.at(i), splitT.at(i)); }
        float splitMinL(std::min(splitL[0], splitL[1])), splitMaxL(std::max(splitL[0], splitL[1]));

        // Too close to cluster endpoint?
        for (int i = 0; i < 2; ++i)
        {
            const float minSep((clusterMin - splitPositions.at(i)).GetMagnitude());
            const float maxSep((clusterMax - splitPositions.at(i)).GetMagnitude());

            if ((minSep < m_endpointBuffer) || (maxSep < m_endpointBuffer))
                rejected.push_back(i);
        }

        if (rejected.size() == 2)
            continue;

        // Reject if contaminant lives inside another particle
        bool contained(false);

        for (const MCContaminant &testContaminant : mcContaminants)
        {
            if (thisContaminant.m_pMCParticle == testContaminant.m_pMCParticle)
                continue;

            // ATTN: use contributing endpoints here to consider hit share regions
            const CartesianPointVector testPositions({testContaminant.m_contStartPosition, testContaminant.m_contEndPosition});
            FloatVector testL({-1.f, -1.f}), testT({-1.f, -1.f});
            for (int i = 0; i < 2; ++i) { clusterFit.GetLocalPosition(testPositions.at(i), testL.at(i), testT.at(i)); }
            float testMinL(std::min(testL[0], testL[1])), testMaxL(std::max(testL[0], testL[1]));

            if ((splitMinL > testMinL) && (splitMaxL < testMaxL))
            {
                contained = true;
                break;
            }
        }
        
        if (contained)
            continue;

        // Reject collinear
        bool collinear(false);
        for (const MCContaminant &testContaminant : mcContaminants)
        {
            if (thisContaminant.m_pMCParticle == testContaminant.m_pMCParticle)
                continue;

            const CartesianPointVector testPositions({testContaminant.m_startPosition, testContaminant.m_endPosition});
            const CartesianPointVector testDirections({testContaminant.m_startDirection, testContaminant.m_endDirection});

            for (int iCurrent = 0; iCurrent < 2; ++iCurrent)
            {
                if (std::find(rejected.begin(), rejected.end(), iCurrent) != rejected.end())
                    continue;
                    
                 for (int iTest = 0; iTest < 2; ++iTest)
                 {
                     // Evaluate collinearity at start-end
                     if (iCurrent == iTest)
                         continue;

                     const float endpointSep((splitPositions.at(iCurrent) - testPositions.at(iTest)).GetMagnitude());
                     float openingAngle(splitDirections.at(iCurrent).GetOpeningAngle(testDirections.at(iTest)));
                     openingAngle *= (180.f / 3.14);
                                            
                     if ((endpointSep < m_collinearMaxSep) && 
                         ((openingAngle < m_collinearMaxOpeningAngle) || (openingAngle > (180.f - m_collinearMaxOpeningAngle))))
                     {
                         collinear = true;
                     }
                 }

                 if (collinear)
                     rejected.push_back(iCurrent);
            }
        }

        // Finally, add split positions
        for (int i = 0; i < 2; ++i)
        {
            if (std::find(rejected.begin(), rejected.end(), i) == rejected.end())
                trueSplitPositions.push_back(splitPositions.at(i));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDClusterSplittingAlgorithm::FillTree(const CartesianPointVector &splitPositions, const MCParticleToHitListMap &mainHitListMap, 
    const TwoDSlidingFitResult &clusterFit, const Features &features)
{
    const bool isContaminated(!splitPositions.empty());

    // Work out who truly owns this cluster
    int highestNHits(-1), matchedPDG(-1);
    for (const auto &entry : mainHitListMap)
    {
        if (static_cast<int>(entry.second.size()) > highestNHits)
        {
            highestNHits = entry.second.size();
            matchedPDG = entry.first->GetParticleId();
        }
    }

    // Don't want to to train on very contaminated tracks
    int nContaminants(0);
    for (auto &entry : mainHitListMap)
        if (entry.second.size() >= m_reducedMinTargetMCHits)
            ++nContaminants;

    // Determine l of splitting positions
    FloatVector vertexL;
    if (isContaminated)
    {
        for (const CartesianVector &splitPosition : splitPositions)
        {
            float thisL(0.f), thisT(0.f);
            clusterFit.GetLocalPosition(splitPosition, thisL, thisT);
            vertexL.push_back(thisL);
        }
    }

    // Fill tree
    FloatVector energy(*features.at("Energy").GetSequence()), hitWidth(*features.at("HitWidth").GetSequence()), 
        gapSep(*features.at("GapSep").GetSequence()), eventHitSep(*features.at("EventHitSep").GetSequence()), 
        clusterHitSep(*features.at("ClusterHitSep").GetSequence()), longitudinal(*features.at("Longitudinal").GetSequence()),
        transverse(*features.at("Transverse").GetSequence()), theta(*features.at("Theta").GetSequence()),
        secVertex(*features.at("SecVertex").GetSequence());

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NContaminants", nContaminants));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "BacktrackedPDG", matchedPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsContaminated", (isContaminated ? 1 : 0)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VertexL", &vertexL));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Energy", &energy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "HitWidth", &hitWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "GapSep", &gapSep));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventHitSep", &eventHitSep));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClusterHitSep", &clusterHitSep));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Longitudinal", &longitudinal));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Transverse", &transverse));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Angle", &theta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecVertex", &secVertex));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
}
#endif
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CaloHitListNameU", m_caloHitListNameU));
    if (m_caloHitListNameU.empty())
        m_caloHitListNameU = "CaloHitListU";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CaloHitListNameV", m_caloHitListNameV));
    if (m_caloHitListNameV.empty())
        m_caloHitListNameV = "CaloHitListV";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CaloHitListNameW", m_caloHitListNameW));
    if (m_caloHitListNameW.empty())
        m_caloHitListNameW = "CaloHitListW";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    if (m_clusterListNameU.empty())
        m_clusterListNameU = "ClustersU";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    if (m_clusterListNameV.empty())
        m_clusterListNameV = "ClustersV";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));
    if (m_clusterListNameW.empty())
        m_clusterListNameW = "ClustersW";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
            XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

        if (m_mcParticleListName.empty())
            m_mcParticleListName = "Input";
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "NuVertexListName", m_nuVertexListName));
    if (m_nuVertexListName.empty())
        m_nuVertexListName = "NeutrinoVertices3D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SecVertexListName", m_secVertexListName));
    if (m_secVertexListName.empty())
        m_secVertexListName = "SecondaryVertices3D";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SlidingWindow", m_slidingWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "LBinSize", m_lBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "KalmanDelta", m_kalmanDelta));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "KalmanProcessVarCoeff", m_kalmanProcessVarCoeff));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "KalmanMeasurementVarCoeff", m_kalmanMeasurementVarCoeff));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SearchRegion1D", m_searchRegion1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "WindowLength", m_windowLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "BackgroundThreshold", m_bkgThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "IsContaminatedThreshold", m_isContaminatedThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "IsSplitThreshold", m_isSplitThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "WindowLength", m_windowLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinNuVertexSep", m_minNuVertexSep));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "GapSearchRegion", m_gapSearchRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ThreeViewMatchMaxXSpan", m_threeViewMatchMaxXSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ThreeViewMatchMaxChi2", m_threeViewMatchMaxChi2));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TwoViewMatchMaxXSpan", m_twoViewMatchMaxXSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TwoViewMatchSearchRegion", m_twoViewMatchSearchRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ContFraction", m_contFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MainFraction", m_mainFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinTargetMCHits", m_minTargetMCHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "EndpointBuffer", m_endpointBuffer));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "CollinearMaxSep", m_collinearMaxSep));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "CollinearMaxOpeningAngle", m_collinearMaxOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ReducedMinTargetMCHits", m_reducedMinTargetMCHits));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "WindowModelName", m_windowModelName));
    m_windowModelName = LArFileHelper::FindFileInPath(m_windowModelName, "FW_SEARCH_PATH");
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_windowModelName, m_windowModel));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SplitPosModelName", m_splitPosModelName));
    m_splitPosModelName = LArFileHelper::FindFileInPath(m_splitPosModelName, "FW_SEARCH_PATH");
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_splitPosModelName, m_splitPosModel));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
