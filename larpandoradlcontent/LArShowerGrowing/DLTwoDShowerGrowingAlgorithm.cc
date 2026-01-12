/**
 *  @file   larpandoradlcontent/LArShowerGrowing/DLTwoDShowerGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning shower growing algorithm
 *
 *  $Log: $
 */

#include <torch/torch.h>

#include "larpandoradlcontent/LArShowerGrowing/DLTwoDShowerGrowingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLTwoDShowerGrowingAlgorithm::HitFeatures::HitFeatures() :
    m_xRel{0.f},
    m_zRel{0.f},
    m_rRel{0.f},
    m_cosThetaRel{0.f},
    m_sinThetaRel{0.f},
    m_distToXGap{0.f},
    m_xWidth{0.f},
    m_energy{0.f}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DLTwoDShowerGrowingAlgorithm::ClusterGroup::ClusterGroup() :
    m_clusters{std::unordered_set<const Cluster *>()},
    m_representativeCluster{nullptr}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DLTwoDShowerGrowingAlgorithm::ClusterGroup::insert(const Cluster* pCluster)
{
    if (m_clusters.empty())
    {
        m_representativeCluster = pCluster;
    }
    m_clusters.insert(pCluster);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DLTwoDShowerGrowingAlgorithm::DLTwoDShowerGrowingAlgorithm() :
    m_polarRScaleFactor{1.f},
    m_cartesianXScaleFactor{1.f},
    m_cartesianZScaleFactor{1.f},
    m_detectorXGaps{std::set<double>{}},
    m_deltaRayLengthThresholdSquared{std::map<HitType, float>{}},
    m_deltaRayParentWeightThreshold{0.f},
    m_hitFeatureDim{12},
    m_trainingMode{false},
    m_trainingTreeName{""},
    m_similarityThreshold{0.5f},
    m_similarityThresholdBeta{1.f},
    m_accessoryClustersMaxHits{2},
    m_maxIterations{1},
    m_includeHitCardinalityFeatures{true},
    m_includeHitNIterationNumFeature{false}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DLTwoDShowerGrowingAlgorithm::~DLTwoDShowerGrowingAlgorithm()
{
    if (m_trainingMode)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName, m_trainingFileName, "UPDATE"));

        for ([[maybe_unused]] const HitType &view : { TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W })
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName + "_view_data", "view", static_cast<int>(view)));
            PANDORA_MONITORING_API(SetTreeVariable(
              this->GetPandora(), m_trainingTreeName + "_view_data", "pitch", LArGeometryHelper::GetWirePitch(this->GetPandora(), view)));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName + "_view_data"));
        }
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName + "_view_data", m_trainingFileName, "UPDATE"));

        const LArGeometryHelper::DetectorBoundaries detBounds{LArGeometryHelper::GetDetectorBoundaries(this->GetPandora())};
        std::cout << "x: " << detBounds.m_xBoundaries.first << " - " << detBounds.m_xBoundaries.second << "\n";
        std::cout << "z: " << detBounds.m_zBoundaries.first << " - " << detBounds.m_zBoundaries.second << "\n";
        std::cout << "y: " << detBounds.m_yBoundaries.first << " - " << detBounds.m_yBoundaries.second << "\n";
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::Run()
{
    if (m_trainingMode)
    {
        return this->PrepareTrainingSample();
    }
    else
    {
        return this->Infer();
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::PrepareTrainingSample()
{
    const std::map<HitType, CartesianVector> viewToVtxPos{this->Get2DVertices()};

    const ClusterList clusterList{this->GetAllClusters()};

    int currClusterID{0}, currMCID{0};
    std::map<const MCParticle* const, int> mcToID = { { nullptr, -1 } }, mcToPDG = { { nullptr, 0 } };
    std::vector<int> clusterView, clusterID;
    std::vector<int> mcID = { -1 }, mcPDG = { 0 };
    std::vector<int> hitClusterID;
    std::vector<float> hitXRelPos, hitZRelPos;
    std::vector<float> hitRRelPos, hitSinThetaRelPos, hitCosThetaRelPos;  // Polar coordinates
    std::vector<float> hitXWidth;
    std::vector<float> hitDistToXGap;
    std::vector<float> hitEnergy;
    std::vector<int> hitMCID;

    std::map<const MCParticle *const, const MCParticle *const> mcFoldTo;

    for (const Cluster *const pCluster : clusterList)
    {
        const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
        clusterView.emplace_back(static_cast<int>(view));
        clusterID.emplace_back(currClusterID);

        CaloHitList clusterCaloHits;
        LArClusterHelper::GetAllHits(pCluster, clusterCaloHits);
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            HitFeatures hitFeatures;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CalculateHitFeatures(pCaloHit, viewToVtxPos.at(view), hitFeatures));

            const MCParticle *const pMainMC{this->GetMainMC(pCaloHit, mcFoldTo)};
            if (mcToID.find(pMainMC) == mcToID.end())
            {
                mcToID.insert({ pMainMC, currMCID++ });
                mcToPDG.insert({ pMainMC, pMainMC->GetParticleId() });
                mcID.emplace_back(mcToID.at(pMainMC));
                mcPDG.emplace_back(mcToPDG.at(pMainMC));
            }

            hitClusterID.emplace_back(currClusterID);

            hitXRelPos.emplace_back(hitFeatures.m_xRel);
            hitZRelPos.emplace_back(hitFeatures.m_zRel);
            hitRRelPos.emplace_back(hitFeatures.m_rRel);
            hitCosThetaRelPos.emplace_back(hitFeatures.m_cosThetaRel);
            hitSinThetaRelPos.emplace_back(hitFeatures.m_sinThetaRel);
            hitXWidth.emplace_back(hitFeatures.m_xWidth);
            hitDistToXGap.emplace_back(hitFeatures.m_distToXGap);
            hitEnergy.emplace_back(hitFeatures.m_energy);

            hitMCID.emplace_back(mcToID.at(pMainMC));
        }
        currClusterID++;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "cluster_id", &clusterID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "cluster_view", &clusterView));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_id", &mcID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_pdg", &mcPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_cluster_id", &hitClusterID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_rel_pos", &hitXRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_z_rel_pos", &hitZRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_r_rel_pos", &hitRRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_ctheta_rel_pos", &hitCosThetaRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_stheta_rel_pos", &hitSinThetaRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_width", &hitXWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_gap_dist", &hitDistToXGap));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_energy", &hitEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_mc_id", &hitMCID));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* DLTwoDShowerGrowingAlgorithm::GetMainMC(
    const CaloHit *const pCaloHit, std::map<const MCParticle *const, const MCParticle *const> &mcFoldTo) const
{
    MCParticleWeightMap weightMap{pCaloHit->GetMCParticleWeightMap()};
    MCParticleWeightMap foldedWeightMap;
    for (const auto &[pMC, weight] : weightMap)
    {
        const MCParticle *pFoldedMC{nullptr};
        if (mcFoldTo.find(pMC) != mcFoldTo.end())
        {
            pFoldedMC = mcFoldTo.at(pMC);
        }
        else
        {
            pFoldedMC = this->FoldMCTo(pMC);
            mcFoldTo.insert({pMC, pFoldedMC});
        }
        foldedWeightMap[pFoldedMC] += weight;
    }
    weightMap = foldedWeightMap;

    const MCParticle *pMainMC{nullptr};
    float maxWeight{0.f};
    for (const auto &[pMC, weight] : weightMap)
    {
        if (weight > maxWeight)
        {
            pMainMC = pMC;
            maxWeight = weight;
        }
        if (weight == maxWeight) // tie-breaker (very unlikely)
        {
            if (LArMCParticleHelper::SortByMomentum(pMC, pMainMC))
            {
                pMainMC = pMC;
            }
        }
    }

    if (pMainMC)
    {
        pMainMC = this->FoldPotentialDeltaRayTo(pCaloHit, pMainMC);
    }

    return pMainMC;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* DLTwoDShowerGrowingAlgorithm::FoldMCTo(const MCParticle *const pMC) const
{
    if (!this->IsEM(pMC))
    {
        return pMC;
    }

    const MCParticle *pCurrentMC{pMC};
    const MCParticle *pLeadingMC{pMC};
    while (!pCurrentMC->IsRootParticle())
    {
        const MCParticle *const pParentMC{*(pCurrentMC->GetParentList().begin())};
        const int parentPdg{std::abs(pParentMC->GetParticleId())};
        if (parentPdg == PHOTON || parentPdg == E_MINUS)
        {
            pCurrentMC = pParentMC;
            continue;
        }
        pLeadingMC = pCurrentMC;
        break;
    }

    return pLeadingMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* DLTwoDShowerGrowingAlgorithm::FoldPotentialDeltaRayTo(const CaloHit *const pCaloHit, const MCParticle *const pMC) const
{
    // Not an electron -> not a delta ray -> do nothing
    if (pMC->IsRootParticle() || pMC->GetParticleId() != E_MINUS)
    {
        return pMC;
    }

    // Did not come from a track-like particle -> not a delta ray -> do nothing
    const MCParticle *const pParentMC{*(pMC->GetParentList().begin())};
    const int parentPdg{std::abs(pParentMC->GetParticleId())};
    if (parentPdg == PHOTON || parentPdg == E_MINUS || PdgTable::GetParticleCharge(parentPdg) == 0)
    {
        return pMC;
    }

    // Delta ray that does not start a shower and is short -> fold into parent particle
    if (!this->CausesShower(pMC, 0) &&
        (pMC->GetVertex() - pMC->GetEndpoint()).GetMagnitudeSquared() < m_deltaRayLengthThresholdSquared.at(pCaloHit->GetHitType()))
    {
        return pParentMC;
    }

    // Now have a delta ray that we would like to cluster but only the hits that are not overlapping with the parent particle
    float parentWeight{std::numeric_limits<float>::lowest()};
    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
    for (const auto &[pContributingMC, weight] : weightMap)
    {
        if (pContributingMC == pParentMC)
        {
            parentWeight = weight;
            break;
        }
    }
    if (parentWeight > m_deltaRayParentWeightThreshold)
    {
        return pParentMC;
    }
    return pMC;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline bool DLTwoDShowerGrowingAlgorithm::IsEM(const pandora::MCParticle *const pMC) const
{
    const int pdg{std::abs(pMC->GetParticleId())};
    return (pdg == E_MINUS || pdg == PHOTON);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DLTwoDShowerGrowingAlgorithm::CausesShower(const MCParticle *const pMC, int nDescendentElectrons) const
{
    if (nDescendentElectrons > 1)
    {
        return true;
    }

    if (std::abs(pMC->GetParticleId()) == E_MINUS)
    {
        nDescendentElectrons++; // Including the parent particle, ie. the first in the recursion, as a descendent
    }
    for (const MCParticle *pChildMC : pMC->GetDaughterList())
    {
        if (this->CausesShower(pChildMC, nDescendentElectrons))
        {
            return true;
        }
    }

    return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::map<HitType, CartesianVector> DLTwoDShowerGrowingAlgorithm::Get2DVertices() const
{
    const VertexList *pVertexList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    const Vertex *const pVertex{pVertexList->front()};
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, pVertex->GetVertexType() != VERTEX_3D);
    return {
        { TPC_VIEW_U, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U) },
        { TPC_VIEW_V, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V) },
        { TPC_VIEW_W, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W) }};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

ClusterList DLTwoDShowerGrowingAlgorithm::GetAllClusters() const
{
    ClusterList clusterList;
    for (const std::string &listName : m_clusterListNames)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetClusters(listName, clusterList));
    }
    return clusterList;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::GetClusters(const std::string clusterListName, ClusterList &clusterList) const
{
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));
    clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::CalculateHitFeatures(
    const CaloHit *const pCaloHit, CartesianVector vtxPos, HitFeatures &hitFeatures) const
{
    const double x{static_cast<double>(pCaloHit->GetPositionVector().GetX())};
    const double xRel{x - static_cast<double>(vtxPos.GetX())};
    const double z{static_cast<double>(pCaloHit->GetPositionVector().GetZ())};
    const double zRel{z - static_cast<double>(vtxPos.GetZ())};

    const double rRel{std::sqrt(pow(xRel, 2.) + pow(zRel, 2.))};
    const double cosThetaRel{ rRel != 0. ? xRel / rRel : 0. };
    const double sinThetaRel{ rRel != 0. ? zRel / rRel : 0. };

    double distToXGap{std::numeric_limits<double>::max()};
    for (const double xGap : m_detectorXGaps)
    {
        const double dist{x - xGap};
        if (std::abs(dist) < std::abs(distToXGap))
        {
            distToXGap = dist;
        }
    }

    hitFeatures.m_xRel = static_cast<float>(xRel);
    hitFeatures.m_zRel = static_cast<float>(zRel);
    hitFeatures.m_rRel = static_cast<float>(rRel);
    hitFeatures.m_cosThetaRel = static_cast<float>(cosThetaRel);
    hitFeatures.m_sinThetaRel = static_cast<float>(sinThetaRel);
    hitFeatures.m_distToXGap = static_cast<float>(distToXGap);
    hitFeatures.m_xWidth = pCaloHit->GetCellSize1();
    hitFeatures.m_energy = pCaloHit->GetMipEquivalentEnergy();

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::Infer()
{
    const std::map<HitType, CartesianVector> viewToVtxPos{this->Get2DVertices()};

    for (const std::string &listName : m_clusterListNames)
    {
        unsigned int nIterations{0};
        while (nIterations++ < m_maxIterations)
        {
            ClusterList clusterList;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetClusters(listName, clusterList));
            if (clusterList.empty())
            {
                break;
            }
            const HitType view{LArClusterHelper::GetClusterHitType(clusterList.front())};

            // Use the trained model to predict pairwise cluster similarities
            SimilarityMatrix clusterSimMat;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                this->PredictClusterSimilarityMatrix(clusterList, view, viewToVtxPos.at(view), nIterations, clusterSimMat));

            // Use the predicted similarities to partition the clusters into super-clusters
            std::vector<ClusterGroup> clusterGroups;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                this->ClusterFromSimilarity(clusterSimMat, m_similarityThreshold, clusterGroups));
            PANDORA_RETURN_IF(STATUS_CODE_FAILURE, this->IsInvalidPartition(clusterGroups, clusterList)); // Sanity check

            if (this->IsSingletonPartition(clusterGroups))
            {
                break;
            }

            if (m_similarityThresholdBeta >= 1.f) // Implies stage 2 is switched off
            {
                // Merge the Clusters according to their prior partition into super-clusters
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->MergeGroups(clusterGroups, listName));
                continue;
            }

            // Optionally, the predicted similarities are reused here in a second stage of clustering that uses a reduced similarity
            // threshold. The idea is to encourage merges for any clusters that underwent minimal/no merging in the first stage.

            // Get the similarities for the super-clusters
            // These similarities are the minimum predicted similarities for all pairs of the constituent clusters of the super-clusters
            SimilarityMatrix superClusterSimMat;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                this->PopulateSuperClusterSimilarityMatrix(clusterGroups, clusterSimMat, superClusterSimMat));

            // Do the first stage merges, the super-clusters are now just the clusters
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->MergeGroups(clusterGroups, listName));
            clusterGroups.clear();
            clusterSimMat = superClusterSimMat;

            // Second stage of clustering using the super-cluster similarity matrix
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                this->ClusterFromSimilarity(clusterSimMat, m_similarityThreshold * m_similarityThresholdBeta, clusterGroups));

            // Do the second stage merges
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->MergeGroups(clusterGroups, listName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::PredictClusterSimilarityMatrix(
    const ClusterList &clusterList,
    const HitType view,
    const CartesianVector &vtxPos,
    const unsigned int iterationNum,
    SimilarityMatrix &clusterSimMat)
{
    torch::InferenceMode guard;

    std::vector<torch::Tensor> tensorEncodedClusters;
    for (const Cluster *const pCluster : clusterList)
    {
        PANDORA_RETURN_IF(STATUS_CODE_NOT_INITIALIZED, !pCluster || pCluster->GetNCaloHits() == 0);
        PANDORA_RETURN_IF(STATUS_CODE_FAILURE, LArClusterHelper::GetClusterHitType(pCluster) != view);

        CaloHitList clusterCaloHits;
        LArClusterHelper::GetAllHits(pCluster, clusterCaloHits);

        std::vector<HitFeatures> clusterFeatures;
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            PANDORA_RETURN_IF(STATUS_CODE_NOT_INITIALIZED, !pCaloHit);
            HitFeatures hitFeatures;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CalculateHitFeatures(pCaloHit, vtxPos, hitFeatures));
            clusterFeatures.emplace_back(hitFeatures);
        }

        torch::Tensor tensorCluster;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            this->MakeClusterTensor(clusterFeatures, view, clusterList.size(), iterationNum, tensorCluster));
        torch::Tensor tensorEncodedCluster{m_modelEncoder.forward({tensorCluster}).toTensor()};
        tensorEncodedClusters.emplace_back(tensorEncodedCluster);
    }
    torch::Tensor tensorEncodedEvent{torch::cat(tensorEncodedClusters, 1)};
    tensorEncodedClusters.clear(); // Free memory

    torch::Tensor tensorAttnEvent{m_modelAttn.forward({tensorEncodedEvent}).toTensor()};
    tensorEncodedEvent = torch::Tensor(); // Free memory

    torch::Tensor tensorSimMat{m_modelSim.forward({tensorAttnEvent}).toTensor()};
    tensorAttnEvent = torch::Tensor();

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PopulateClusterSimilarityMatrix(tensorSimMat, clusterList, clusterSimMat));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::MakeClusterTensor(
    const std::vector<HitFeatures> &clusterFeatures,
    const HitType view,
    const size_t nClusters,
    const unsigned int iterationNum,
    torch::Tensor &tensorCluster) const
{
    const int nHits{static_cast<int>(clusterFeatures.size())};

    const float nHitsFeat{std::log(static_cast<float>(nHits))};
    const float nClustersFeat{std::log(static_cast<float>(nClusters))};
    const float iterationNumFeat{std::log(static_cast<float>(iterationNum))};
    const float zWidthFeat{LArGeometryHelper::GetWirePitch(this->GetPandora(), view) * m_cartesianZScaleFactor};

    tensorCluster = torch::zeros({1, nHits, m_hitFeatureDim});
    auto accessor = tensorCluster.accessor<float, 3>();

    for (int i = 0; i < nHits; i++)
    {
        const HitFeatures hitFeatures{clusterFeatures.at(i)};
        accessor[0][i][0] = hitFeatures.m_rRel * m_polarRScaleFactor;
        accessor[0][i][1] = hitFeatures.m_cosThetaRel;
        accessor[0][i][2] = hitFeatures.m_sinThetaRel;
        accessor[0][i][3] = hitFeatures.m_xRel * m_cartesianXScaleFactor;
        accessor[0][i][4] = hitFeatures.m_zRel * m_cartesianZScaleFactor;
        accessor[0][i][5] = hitFeatures.m_xWidth * m_cartesianXScaleFactor;
        accessor[0][i][6] = zWidthFeat;
        accessor[0][i][7] = hitFeatures.m_distToXGap * m_cartesianXScaleFactor;
        accessor[0][i][8] = std::log(hitFeatures.m_energy);
        accessor[0][i][9] = 0.f;
        accessor[0][i][10] = 0.f;
        accessor[0][i][11] = 0.f;
        if (view == TPC_VIEW_U)
        {
            accessor[0][i][9] = 1.f;
        }
        else if (view == TPC_VIEW_V)
        {
            accessor[0][i][10] = 1.f;
        }
        else
        {
            accessor[0][i][11] = 1.f;
        }
        if (m_includeHitCardinalityFeatures)
        {
            accessor[0][i][m_hitFeaturesNHitsIdx] = nHitsFeat;
            accessor[0][i][m_hitFeaturesNClustersIdx] = nClustersFeat;
        }
        if (m_includeHitNIterationNumFeature)
        {
            accessor[0][i][m_hitFeaturesIterationNumIdx] = iterationNumFeat;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::PopulateClusterSimilarityMatrix(
    const torch::Tensor &tensorSimMat, const ClusterList &clusterList, SimilarityMatrix &clusterSimMat) const
{
    PANDORA_RETURN_IF(STATUS_CODE_NOT_ALLOWED, !clusterSimMat.empty());
    const int64_t nClusters{static_cast<int64_t>(clusterList.size())};
    PANDORA_RETURN_IF(STATUS_CODE_NOT_ALLOWED,
        tensorSimMat.dim() != 3 || nClusters != tensorSimMat.size(-1) || nClusters != tensorSimMat.size(-2));

    auto accessor = tensorSimMat.accessor<float, 3>();

    auto iterI{clusterList.begin()};
    for (int i = 0; i < nClusters; i++, iterI++)
    {
        auto iterJ{clusterList.begin()};
        for (int j = 0; j < nClusters; j++, iterJ++)
        {
            clusterSimMat[*iterI][*iterJ] = accessor[0][i][j];
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::ClusterFromSimilarity(
    const SimilarityMatrix &clusterSimMat, const float similarityThreshold, std::vector<ClusterGroup> &clusterGroups) const
{
    AdjacencyLists coreClusterAdjLists, accClusterAdjLists;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->PopulateAdjacencyLists(clusterSimMat, similarityThreshold, coreClusterAdjLists, accClusterAdjLists));

    // Connected grouping of core clusters
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CalculateConnectedGroups(coreClusterAdjLists, clusterGroups));

    // Add accessory clusters into the core groups
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        this->AddClustersToGroups(clusterSimMat, accClusterAdjLists, similarityThreshold, clusterGroups));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CalculateConnectedGroups(accClusterAdjLists, clusterGroups));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::PopulateAdjacencyLists(
    const SimilarityMatrix &simMat, const float similarityThreshold, AdjacencyLists &coreAdjLists, AdjacencyLists &accAdjLists) const
{
    PANDORA_RETURN_IF(STATUS_CODE_NOT_ALLOWED, !coreAdjLists.empty() || !accAdjLists.empty());

    for (const auto &[pClusterI, simRow] : simMat)
    {
        const bool iIsCore{pClusterI->GetNCaloHits() > m_accessoryClustersMaxHits};
        // I rely on all clusters appearing in an adjList, this ensures singletons appear
        if (iIsCore && coreAdjLists.find(pClusterI) == coreAdjLists.end())
        {
            coreAdjLists[pClusterI] = {};
        }
        else if (!iIsCore && accAdjLists.find(pClusterI) == accAdjLists.end())
        {
            accAdjLists[pClusterI] = {};
        }

        for (const auto &[pClusterJ, sim] : simRow)
        {
            if (pClusterI == pClusterJ || sim <= similarityThreshold)
            {
                continue;
            }
            const bool jIsCore{pClusterJ->GetNCaloHits() > m_accessoryClustersMaxHits};
            if (iIsCore && jIsCore)
            {
                coreAdjLists.at(pClusterI).emplace_back(pClusterJ);
            }
            else if (!iIsCore && !jIsCore)
            {
                accAdjLists.at(pClusterI).emplace_back(pClusterJ);
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::CalculateConnectedGroups(
    const AdjacencyLists &clusterAdjLists, std::vector<ClusterGroup> &clusterGroups) const
{
    std::unordered_set<const Cluster *> visitedClusters;
    for (const auto &[pClusterRoot, _] : clusterAdjLists)
    {
        if (visitedClusters.find(pClusterRoot) != visitedClusters.end())
        {
            continue;
        }

        std::vector<const Cluster *> toSearch = { pClusterRoot };
        ClusterGroup group;
        while (!toSearch.empty())
        {
            const Cluster *const pClusterI{toSearch.back()};
            toSearch.pop_back();
            group.insert(pClusterI);
            visitedClusters.insert(pClusterI);

            for (const Cluster *const pClusterJ : clusterAdjLists.at(pClusterI))
            {
                if (visitedClusters.find(pClusterJ) != visitedClusters.end())
                {
                    continue;
                }
                toSearch.push_back(pClusterJ);
            }
        }

        clusterGroups.emplace_back(group);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::AddClustersToGroups(
    const SimilarityMatrix &clusterSimMat,
    AdjacencyLists &ungroupedClusterAdjLists,
    const float similarityThreshold,
    std::vector<ClusterGroup> &clusterGroups) const
{
    std::map<const Cluster *const, int> plannedMerges;
    do
    {
        plannedMerges.clear();
        for (const auto &[pClusterUngrouped, _] : ungroupedClusterAdjLists)
        {
            int bestGroup{-1};
            float bestSim{std::numeric_limits<float>::lowest()};

            auto iter{clusterGroups.begin()};
            for (int i = 0; i < static_cast<int>(clusterGroups.size()); i++, iter++)
            {
                for (const Cluster *const pClusterGrouped : *iter)
                {
                    const float sim{clusterSimMat.at(pClusterUngrouped).at(pClusterGrouped)};
                    if (sim > similarityThreshold && sim > bestSim)
                    {
                        bestSim = sim;
                        bestGroup = i;
                    }
                }
            }

            if (bestGroup >= 0)
            {
                plannedMerges[pClusterUngrouped] = bestGroup;
            }
        }

        for (const auto &[pClusterUngrouped, bestGroup] : plannedMerges)
        {
            clusterGroups[bestGroup].insert(pClusterUngrouped);
            ungroupedClusterAdjLists.erase(pClusterUngrouped);
        }
    } while (!plannedMerges.empty());

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::PopulateSuperClusterSimilarityMatrix(
    const std::vector<ClusterGroup> &clusterGroups, const SimilarityMatrix &clusterSimMat, SimilarityMatrix &superClusterSimMat) const
{
    for (const ClusterGroup &clusterGroupI : clusterGroups)
    {
        const Cluster *const pClusterRepI{clusterGroupI.GetRepresentativeCluster()};
        for (const ClusterGroup &clusterGroupJ : clusterGroups)
        {
            const Cluster *const pClusterRepJ{clusterGroupJ.GetRepresentativeCluster()};
            if (pClusterRepI == pClusterRepJ)
            {
                superClusterSimMat[pClusterRepI][pClusterRepJ] = 1.f;
                continue;
            }
            float minSim{1.f};
            for (const Cluster *const pClusterI : clusterGroupI)
            {
                for (const Cluster *const pClusterJ : clusterGroupJ)
                {
                    minSim = std::min(minSim, clusterSimMat.at(pClusterI).at(pClusterJ));
                }
            }
            superClusterSimMat[pClusterRepI][pClusterRepJ] = minSim;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::MergeGroups(const std::vector<ClusterGroup> &clusterGroups, const std::string &listName) const
{
    for (const ClusterGroup &clusterGroup : clusterGroups)
    {
        if (clusterGroup.empty())
        {
            continue;
        }

        auto iter{clusterGroup.begin()};
        const Cluster *const pClusterToEnlarge{clusterGroup.GetRepresentativeCluster()};
        for (; iter != clusterGroup.end(); ++iter)
        {
            if (*iter == pClusterToEnlarge)
            {
                continue;
            }
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, *iter, listName, listName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DLTwoDShowerGrowingAlgorithm::IsInvalidPartition(const std::vector<ClusterGroup> &clusterGroups, const ClusterList &clusterList) const
{
    std::unordered_set<const Cluster *> uniqueGroupedClusters;
    size_t totalGroupedClusters{0};
    for (const ClusterGroup &group : clusterGroups)
    {
        uniqueGroupedClusters.insert(group.begin(), group.end());
        totalGroupedClusters += group.size();
    }
    return totalGroupedClusters != clusterList.size() || uniqueGroupedClusters.size() != clusterList.size();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DLTwoDShowerGrowingAlgorithm::IsSingletonPartition(const std::vector<ClusterGroup> &clusterGroups) const
{
    for (const ClusterGroup &group : clusterGroups)
    {
        if (group.size() > 1)
        {
            return false;
        }
    }
    return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoDShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingTreeName", m_trainingTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
    }
    else
    {
        std::string modelEncoderName, modelAttnName, modelSimName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelEncoderFileName", modelEncoderName));
        modelEncoderName = LArFileHelper::FindFileInPath(modelEncoderName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelEncoderName, m_modelEncoder);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelAttnFileName", modelAttnName));
        modelAttnName = LArFileHelper::FindFileInPath(modelAttnName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelAttnName, m_modelAttn);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelSimFileName", modelSimName));
        modelSimName = LArFileHelper::FindFileInPath(modelSimName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelSimName, m_modelSim);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF( STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SimilarityThreshold", m_similarityThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SimilarityThresholdBeta", m_similarityThresholdBeta));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "AccessoryClusterMaxHits", m_accessoryClustersMaxHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxIterations", m_maxIterations));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "IncludeHitCardinalityFeatures", m_includeHitCardinalityFeatures));
    if (m_includeHitCardinalityFeatures)
    {
        m_hitFeaturesNHitsIdx = m_hitFeatureDim++;
        m_hitFeaturesNClustersIdx = m_hitFeatureDim++;
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "IncludeHitNIterationNumFeature", m_includeHitNIterationNumFeature));
    if (m_includeHitNIterationNumFeature)
    {
        m_hitFeaturesIterationNumIdx = m_hitFeatureDim++;
    }

    m_deltaRayLengthThresholdSquared = {
        { TPC_VIEW_U, static_cast<float>(std::pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U), 2.)) },
        { TPC_VIEW_V, static_cast<float>(std::pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V), 2.)) },
        { TPC_VIEW_W, static_cast<float>(std::pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W), 2.)) }
    };

    const LArGeometryHelper::DetectorBoundaries detBounds{LArGeometryHelper::GetDetectorBoundaries(this->GetPandora())};
    double xLow{static_cast<double>(detBounds.m_xBoundaries.first)}, xHigh{static_cast<double>(detBounds.m_xBoundaries.second)};
    double yLow{static_cast<double>(detBounds.m_yBoundaries.first)}, yHigh{static_cast<double>(detBounds.m_yBoundaries.second)};
    double zLow{static_cast<double>(detBounds.m_zBoundaries.first)}, zHigh{static_cast<double>(detBounds.m_zBoundaries.second)};
    m_polarRScaleFactor = static_cast<float>(
        1. / std::sqrt(std::pow(xHigh - xLow, 2.) + std::pow(yHigh - yLow, 2.) + std::pow(zHigh - zLow, 2.)));
    m_cartesianXScaleFactor = static_cast<float>(1. / (xHigh - xLow));
    m_cartesianZScaleFactor = static_cast<float>(1. / (zHigh - zLow));

    for (const DetectorGap *const pDetectorGap : this->GetPandora().GetGeometry()->GetDetectorGapList())
    {
        const LineGap *const pLineGap = dynamic_cast<const LineGap *>(pDetectorGap);
        if (pLineGap->GetLineGapType() == TPC_DRIFT_GAP)
        {
            m_detectorXGaps.insert(static_cast<double>(pLineGap->GetLineStartX()));
            m_detectorXGaps.insert(static_cast<double>(pLineGap->GetLineEndX()));
        }
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
