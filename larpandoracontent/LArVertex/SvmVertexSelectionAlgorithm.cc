/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.cc
 *
 *  @brief  Implementation of the Svm vertex selection algorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArSvmHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArVertex/EnergyKickFeatureTool.h"
#include "larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/GlobalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/ShowerAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/RPhiFeatureTool.h"

#include "larpandoracontent/LArVertex/SvmVertexSelectionAlgorithm.h"

#include <random>

using namespace pandora;

namespace lar_content
{
SvmVertexSelectionAlgorithm::SvmVertexSelectionAlgorithm() : VertexSelectionBaseAlgorithm(),
    m_trainingSetMode(false),
    m_allowClassifyDuringTraining(false),
    m_selectInputHits(true),
    m_mcVertexXCorrection(0.f),
    m_minHitSharingFraction(0.9f),
    m_maxPhotonPropagation(2.5f),
    m_minHitsForGoodView(5),
    m_minPrimaryGoodHits(15),
    m_minPrimaryGoodViews(2),
    m_minClusterCaloHits(12),
    m_slidingFitWindow(100),
    m_minShowerSpineLength(15.f),
    m_beamDeweightingConstant(0.4),
    m_localAsymmetryConstant(3.f),
    m_globalAsymmetryConstant(1.f),
    m_showerAsymmetryConstant(1.f),
    m_energyKickConstant(0.06),
    m_showerClusteringDistance(3.f),
    m_minShowerClusterHits(1),
    m_useShowerClusteringApproximation(false),
    m_regionRadius(10.f),
    m_rPhiFineTuningRadius(2.f),
    m_maxTrueVertexRadius(1.f),
    m_useRPhiFeatureForRegion(false),
    m_dropFailedRPhiFastScoreCandidates(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    ClusterList clustersU, clustersV, clustersW;
    this->GetClusterLists(m_inputClusterListNames, clustersU, clustersV, clustersW);

    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    this->CalculateClusterSlidingFits(clustersU, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListU);
    this->CalculateClusterSlidingFits(clustersV, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListV);
    this->CalculateClusterSlidingFits(clustersW, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListW);

    ShowerClusterList showerClusterListU, showerClusterListV, showerClusterListW;
    this->CalculateShowerClusterList(clustersU, showerClusterListU);
    this->CalculateShowerClusterList(clustersV, showerClusterListV);
    this->CalculateShowerClusterList(clustersW, showerClusterListW);

    // Create maps from hit types to objects for passing to feature tools.
    const ClusterListMap clusterListMap{{TPC_VIEW_U, clustersU},
                                        {TPC_VIEW_V, clustersV},
                                        {TPC_VIEW_W, clustersV}};

    const SlidingFitDataListMap slidingFitDataListMap{{TPC_VIEW_U, slidingFitDataListU},
                                                      {TPC_VIEW_V, slidingFitDataListV},
                                                      {TPC_VIEW_W, slidingFitDataListW}};

    const ShowerClusterListMap showerClusterListMap{{TPC_VIEW_U, showerClusterListU},
                                                    {TPC_VIEW_V, showerClusterListV},
                                                    {TPC_VIEW_W, showerClusterListW}};

    const KDTreeMap kdTreeMap{{TPC_VIEW_U, kdTreeU},
                              {TPC_VIEW_V, kdTreeV},
                              {TPC_VIEW_W, kdTreeW}};

    // Calculate the event feature list and the vertex feature map.
    EventFeatureInfo eventFeatureInfo(this->CalculateEventFeatures(clustersU, clustersV, clustersW, vertexVector));

    SupportVectorMachine::DoubleVector eventFeatureList;
    this->AddEventFeaturesToVector(eventFeatureInfo, eventFeatureList);

    VertexFeatureInfoMap vertexFeatureInfoMap;
    for (const Vertex *const pVertex : vertexVector)
    {
        this->PopulateVertexFeatureInfoMap(beamConstants, clusterListMap, slidingFitDataListMap, showerClusterListMap, kdTreeMap, pVertex,
            vertexFeatureInfoMap);
    }

    // Use a simple score to get the list of vertices representing good regions.
    VertexScoreList initialScoreList;
    for (const Vertex *const pVertex : vertexVector)
        PopulateInitialScoreList(vertexFeatureInfoMap, pVertex, initialScoreList);

    VertexVector bestRegionVertices;
    this->GetBestRegionVertices(initialScoreList, bestRegionVertices);

    if (m_trainingSetMode)
        this->ProduceTrainingSets(vertexVector, bestRegionVertices, vertexFeatureInfoMap, eventFeatureList, kdTreeMap);

    if ((!m_trainingSetMode || m_allowClassifyDuringTraining) && !bestRegionVertices.empty())
    {
        // Use svm to choose the region.
        const Vertex *const pBestRegionVertex(this->CompareVertices(bestRegionVertices, vertexFeatureInfoMap, eventFeatureList, m_svMachineRegion,
            m_useRPhiFeatureForRegion));

        // Get all the vertices in the best region.
        VertexVector regionalVertices{pBestRegionVertex};
        for (const Vertex *const pVertex : vertexVector)
        {
            if (pVertex == pBestRegionVertex)
                continue;

            if ((pBestRegionVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude() < m_regionRadius)
                regionalVertices.push_back(pVertex);
        }

        this->CalculateRPhiScores(regionalVertices, vertexFeatureInfoMap, kdTreeMap);

        if (!regionalVertices.empty())
        {
            // Use svm to choose the vertex and then fine-tune using the RPhi score.
            const Vertex *const pBestVertex(this->CompareVertices(regionalVertices, vertexFeatureInfoMap, eventFeatureList, m_svMachineVertex, true));
            this->PopulateFinalVertexScoreList(vertexFeatureInfoMap, pBestVertex, vertexVector, vertexScoreList);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::CalculateShowerClusterList(const ClusterList &inputClusterList, ShowerClusterList &showerClusterList) const
{
    ClusterEndPointsMap clusterEndPointsMap;
    ClusterList showerLikeClusters;
    this->GetShowerLikeClusterEndPoints(inputClusterList, showerLikeClusters, clusterEndPointsMap);

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    ClusterList availableShowerLikeClusters(showerLikeClusters.begin(), showerLikeClusters.end());

    while (!availableShowerLikeClusters.empty())
    {
        ClusterList showerCluster;
        showerCluster.push_back(availableShowerLikeClusters.back());
        availableShowerLikeClusters.pop_back();

        bool addedCluster(true);
        while (addedCluster && !availableShowerLikeClusters.empty())
        {
            addedCluster = false;
            for (const Cluster *const pCluster : showerCluster)
            {
                addedCluster = this->AddClusterToShower(clusterEndPointsMap, availableShowerLikeClusters, pCluster, showerCluster);

                if (addedCluster)
                    break;
            }
        }

        unsigned int totHits(0);
        for (const Cluster *const pCluster : showerCluster)
            totHits += pCluster->GetNCaloHits();

        if (totHits < m_minClusterCaloHits)
            continue;

        showerClusterList.emplace_back(showerCluster, slidingFitPitch, m_slidingFitWindow);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetShowerLikeClusterEndPoints(const ClusterList &clusterList, ClusterList &showerLikeClusters,
    ClusterEndPointsMap &clusterEndPointsMap) const
{
    for (const Cluster *const pCluster : clusterList)
    {
        if (pCluster->GetNCaloHits() < m_minShowerClusterHits)
            continue;

        if (this->IsClusterShowerLike(pCluster))
            showerLikeClusters.push_back(pCluster);

        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        CaloHitVector clusterCaloHitVector(clusterCaloHitList.begin(), clusterCaloHitList.end());
        std::sort(clusterCaloHitVector.begin(), clusterCaloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        if (clusterCaloHitVector.empty())
            continue;

        ClusterEndPoints clusterEndPoints(clusterCaloHitVector.front()->GetPositionVector(), clusterCaloHitVector.back()->GetPositionVector());
        clusterEndPointsMap.emplace(pCluster, clusterEndPoints);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SvmVertexSelectionAlgorithm::AddClusterToShower(const ClusterEndPointsMap &clusterEndPointsMap, ClusterList &availableShowerLikeClusters,
    const Cluster *const pCluster, ClusterList &showerCluster) const
{
    const auto existingEndPointsIter(clusterEndPointsMap.find(pCluster));
    if (existingEndPointsIter == clusterEndPointsMap.end())
        return false;

    const ClusterEndPoints &existingClusterEndPoints(existingEndPointsIter->second);

    for (auto iter = availableShowerLikeClusters.begin(); iter != availableShowerLikeClusters.end(); /* no increment */)
    {
        const auto &newEndPointsIter(clusterEndPointsMap.find(*iter));
        if (newEndPointsIter == clusterEndPointsMap.end())
            continue;

        const ClusterEndPoints &newClusterEndPoints(newEndPointsIter->second);
        bool satisfiesClosestDistance(true);

        if (m_useShowerClusteringApproximation)
        {
            const float startStartDistance((newClusterEndPoints.first - existingClusterEndPoints.first).GetMagnitude());
            const float startEndDistance((newClusterEndPoints.first - existingClusterEndPoints.second).GetMagnitude());
            const float endStartDistance((newClusterEndPoints.second - existingClusterEndPoints.first).GetMagnitude());
            const float endEndDistance((newClusterEndPoints.second - existingClusterEndPoints.second).GetMagnitude());

            const float smallestDistance(std::min(std::min(startStartDistance, startEndDistance), std::min(endStartDistance, endEndDistance)));
            satisfiesClosestDistance = (smallestDistance < m_showerClusteringDistance);
        }

        else
            satisfiesClosestDistance = (LArClusterHelper::GetClosestDistance(pCluster, *iter) < m_showerClusteringDistance);

        if (satisfiesClosestDistance)
        {
            showerCluster.push_back(*iter);
            availableShowerLikeClusters.erase(iter);
            return true;
        }

        ++iter;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

SvmVertexSelectionAlgorithm::EventFeatureInfo SvmVertexSelectionAlgorithm::CalculateEventFeatures(const ClusterList &clusterListU,
    const ClusterList &clusterListV, const ClusterList &clusterListW, const VertexVector &vertexVector) const
{
    float eventEnergy(0.f);
    unsigned int nShoweryHits(0), nHits(0);

    this->IncrementShoweryParameters(clusterListU, nShoweryHits, nHits, eventEnergy);
    this->IncrementShoweryParameters(clusterListV, nShoweryHits, nHits, eventEnergy);
    this->IncrementShoweryParameters(clusterListW, nShoweryHits, nHits, eventEnergy);

    const unsigned int nClusters(clusterListU.size() + clusterListV.size() + clusterListW.size());
    const float eventShoweryness((nHits > 0) ? static_cast<float>(nShoweryHits) / static_cast<float>(nHits) : 0.f);

    ClusterList allClusters(clusterListU);
    allClusters.insert(allClusters.end(), clusterListV.begin(), clusterListV.end());
    allClusters.insert(allClusters.end(), clusterListW.begin(), clusterListW.end());

    float eventVolume(0.f), longitudinality(0.f);
    this->GetEventShapeFeatures(allClusters, eventVolume, longitudinality);

    return EventFeatureInfo(eventShoweryness, eventEnergy, eventVolume, longitudinality, nHits, nClusters, vertexVector.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::IncrementShoweryParameters(const ClusterList &clusterList, unsigned int &nShoweryHits, unsigned int &nHits,
    float &eventEnergy) const
{
    for (const Cluster *const pCluster : clusterList)
    {
        if (this->IsClusterShowerLike(pCluster))
            nShoweryHits += pCluster->GetNCaloHits();

        eventEnergy += pCluster->GetElectromagneticEnergy();
        nHits += pCluster->GetNCaloHits();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool SvmVertexSelectionAlgorithm::IsClusterShowerLike(const Cluster *const pCluster) const
{
    return (pCluster->GetParticleId() == E_MINUS && LArClusterHelper::GetLength(pCluster) < m_minShowerSpineLength);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetEventShapeFeatures(const ClusterList &clusterList, float &eventVolume, float &longitudinality) const
{
    InputFloat xMin, yMin, zMin, xMax, yMax, zMax;

    for (const Cluster *const pCluster : clusterList)
    {
        CartesianVector minPosition(0.f, 0.f, 0.f), maxPosition(0.f, 0.f, 0.f);
        LArClusterHelper::GetClusterBoundingBox(pCluster, minPosition, maxPosition);

        this->UpdateSpanCoordinate(minPosition.GetX(), maxPosition.GetX(), xMin, xMax);
        this->UpdateSpanCoordinate(minPosition.GetY(), maxPosition.GetY(), yMin, yMax);
        this->UpdateSpanCoordinate(minPosition.GetZ(), maxPosition.GetZ(), zMin, zMax);
    }

    const float xSpan(this->GetCoordinateSpan(xMax, xMin));
    const float ySpan(this->GetCoordinateSpan(yMax, zMin));
    const float zSpan(this->GetCoordinateSpan(yMax, zMin));

    // Calculate the volume and longitudinality of the event (ySpan often 0 - to be investigated).
    if ((xSpan > std::numeric_limits<float>::epsilon()) && (ySpan > std::numeric_limits<float>::epsilon()))
    {
        eventVolume     = xSpan * ySpan * zSpan;
        longitudinality = zSpan / (xSpan + ySpan);
    }

    else if (ySpan > std::numeric_limits<float>::epsilon())
    {
        eventVolume     = ySpan * ySpan * zSpan;
        longitudinality = zSpan / (ySpan + ySpan);
    }

    else if (xSpan > std::numeric_limits<float>::epsilon())
    {
        eventVolume     = xSpan * xSpan * zSpan;
        longitudinality = zSpan / (xSpan + xSpan);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void SvmVertexSelectionAlgorithm::UpdateSpanCoordinate(const float minPositionCoord, const float maxPositionCoord, InputFloat &minCoord,
    InputFloat &maxCoord) const
{
    if (!minCoord.IsInitialized() || minPositionCoord < minCoord.Get())
        minCoord = minPositionCoord;

    if (!maxCoord.IsInitialized() || maxPositionCoord > maxCoord.Get())
        maxCoord = maxPositionCoord;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SvmVertexSelectionAlgorithm::GetCoordinateSpan(const InputFloat &minCoord, const InputFloat &maxCoord) const
{
   if (minCoord.IsInitialized() && maxCoord.IsInitialized())
        return std::fabs(maxCoord.Get() - minCoord.Get());

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::AddEventFeaturesToVector(const EventFeatureInfo &eventFeatureInfo,
    SupportVectorMachine::DoubleVector &featureVector) const
{
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventShoweryness));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventEnergy));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventVolume));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_longitudinality));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nHits));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nClusters));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nCandidates));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::PopulateVertexFeatureInfoMap(const BeamConstants &beamConstants, const ClusterListMap &clusterListMap,
    const SlidingFitDataListMap &slidingFitDataListMap, const ShowerClusterListMap &showerClusterListMap, const KDTreeMap &kdTreeMap,
    const Vertex *const pVertex, VertexFeatureInfoMap &vertexFeatureInfoMap) const
{
    float bestFastScore(-std::numeric_limits<float>::max()); // not actually used - artefact of toolizing RPhi score and still using performance trick

    const double beamDeweighting(this->GetBeamDeweightingScore(beamConstants, pVertex));

    const double energyKick(LArSvmHelper::CalculateFeaturesOfType<EnergyKickFeatureTool>(m_featureToolVector, this, pVertex,
        slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));

    const double localAsymmetry(LArSvmHelper::CalculateFeaturesOfType<LocalAsymmetryFeatureTool>(m_featureToolVector, this, pVertex,
        slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));

    const double globalAsymmetry(LArSvmHelper::CalculateFeaturesOfType<GlobalAsymmetryFeatureTool>(m_featureToolVector, this, pVertex,
        slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));

    const double showerAsymmetry(LArSvmHelper::CalculateFeaturesOfType<ShowerAsymmetryFeatureTool>(m_featureToolVector, this, pVertex,
        slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));

    //const double rPhiFeature(LArSvmHelper::CalculateFeaturesOfType<RPhiFeatureTool>(m_featureToolVector, this, pVertex,
    //    slidingFitDataListMap, clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));

    VertexFeatureInfo vertexFeatureInfo(beamDeweighting, 0.f, energyKick, localAsymmetry, globalAsymmetry, showerAsymmetry);
    vertexFeatureInfoMap.emplace(pVertex, vertexFeatureInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::PopulateInitialScoreList(VertexFeatureInfoMap &vertexFeatureInfoMap, const Vertex *const pVertex,
                                                           VertexScoreList &initialScoreList) const
{
    VertexFeatureInfo vertexFeatureInfo = vertexFeatureInfoMap.at(pVertex);

    const float beamDeweightingScore(vertexFeatureInfo.m_beamDeweighting / m_beamDeweightingConstant);
    const float energyKickScore(-vertexFeatureInfo.m_energyKick / m_energyKickConstant);
    const float localAsymmetryScore(vertexFeatureInfo.m_localAsymmetry / m_localAsymmetryConstant);
    const float globalAsymmetryScore(vertexFeatureInfo.m_globalAsymmetry / m_globalAsymmetryConstant);
    const float showerAsymmetryScore(vertexFeatureInfo.m_showerAsymmetry / m_showerAsymmetryConstant);

    initialScoreList.emplace_back(pVertex, beamDeweightingScore + energyKickScore + localAsymmetryScore + globalAsymmetryScore + showerAsymmetryScore);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetBestRegionVertices(VertexScoreList &initialScoreList, VertexVector &bestRegionVertices) const
{
    std::sort(initialScoreList.begin(), initialScoreList.end());

    for (const VertexScore &vertexScore : initialScoreList)
    {
        const Vertex *const pVertex(vertexScore.GetVertex());
        bool farEnoughAway(true);

        for (const Vertex *const pRegionVertex: bestRegionVertices)
        {
            if (pRegionVertex == pVertex)
            {
                farEnoughAway = false;
                break;
            }

            const float distance = (pRegionVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude();

            if (distance <= m_regionRadius)
            {
                farEnoughAway = false;
                break;
            }
        }

        if (farEnoughAway)
            bestRegionVertices.push_back(pVertex);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::ProduceTrainingSets(const VertexVector &vertexVector, const VertexVector &bestRegionVertices,
    VertexFeatureInfoMap &vertexFeatureInfoMap, const SupportVectorMachine::DoubleVector &eventFeatureList, const KDTreeMap &kdTreeMap) const
{
    if (vertexVector.empty())
        return;

    // Create a distribution for random coin flips.
    std::random_device device;
    std::mt19937 generator(device());
    std::bernoulli_distribution coinFlip(0.5);

    std::string interactionType(this->GetInteractionType(vertexVector));

    // Produce training examples for the vertices representing regions.
    const Vertex *const pBestRegionVertex(this->ProduceTrainingExamples(bestRegionVertices, vertexFeatureInfoMap, coinFlip, generator,
        interactionType, m_trainingOutputFileRegion, eventFeatureList, m_regionRadius, m_useRPhiFeatureForRegion));

    // Get all the vertices in the best region.
    VertexVector regionalVertices{pBestRegionVertex};
    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestRegionVertex)
            continue;

        if ((pBestRegionVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude() < m_regionRadius)
            regionalVertices.push_back(pVertex);
    }

    this->CalculateRPhiScores(regionalVertices, vertexFeatureInfoMap, kdTreeMap);

    // Produce training examples for the final vertices within the best region.
    if (!regionalVertices.empty())
    {
        this->ProduceTrainingExamples(regionalVertices, vertexFeatureInfoMap, coinFlip, generator, interactionType, m_trainingOutputFileVertex,
            eventFeatureList, m_maxTrueVertexRadius, true);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::CalculateRPhiScores(VertexVector &vertexVector, VertexFeatureInfoMap &vertexFeatureInfoMap,
    const KDTreeMap &kdTreeMap) const
{
    float bestFastScore(-std::numeric_limits<float>::max());

    for (auto iter = vertexVector.begin(); iter != vertexVector.end(); /* no increment */)
    {
        VertexFeatureInfo &vertexFeatureInfo = vertexFeatureInfoMap.at(*iter);
        vertexFeatureInfo.m_rPhiFeature = static_cast<float>(LArSvmHelper::CalculateFeaturesOfType<RPhiFeatureTool>(m_featureToolVector, this, *iter,
            SlidingFitDataListMap(), ClusterListMap(), kdTreeMap, ShowerClusterListMap(), vertexFeatureInfo.m_beamDeweighting,
            bestFastScore).at(0));

        if (m_dropFailedRPhiFastScoreCandidates && (vertexFeatureInfo.m_rPhiFeature <= std::numeric_limits<float>::epsilon()))
            iter = vertexVector.erase(iter);

        else
            ++iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string SvmVertexSelectionAlgorithm::GetInteractionType(const VertexVector &vertexVector) const
{
    // Extract input collections
    const MCParticleList *pMCParticleList(nullptr);
    PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);

    const CaloHitList *pCaloHitList(nullptr);
    PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList);

    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::SelectTrueNeutrinos(pMCParticleList, mcNeutrinoVector);

    std::string interactionType("UNKNOWN");

    for (const Vertex *const pVertex : vertexVector)
    {
        float mcVertexDr(std::numeric_limits<float>::max());

        for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
        {
            const CartesianVector mcNeutrinoPosition(pMCNeutrino->GetEndpoint().GetX() + m_mcVertexXCorrection, pMCNeutrino->GetEndpoint().GetY(),
                pMCNeutrino->GetEndpoint().GetZ());

            const float dr((mcNeutrinoPosition - pVertex->GetPosition()).GetMagnitude());

            if (dr < mcVertexDr)
            {
                if (const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino))
                {
                    interactionType = LArMCParticleHelper::ToString(LArMCParticleHelper::GetInteractionType(pLArMCNeutrino, pMCParticleList,
                        pCaloHitList, m_minPrimaryGoodHits, m_minHitsForGoodView, m_minPrimaryGoodViews, m_selectInputHits, m_maxPhotonPropagation,
                        m_minHitSharingFraction));
                }

                mcVertexDr = dr;
            }
        }
    }

    return interactionType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Vertex * SvmVertexSelectionAlgorithm::ProduceTrainingExamples(const VertexVector &vertexVector,
    const VertexFeatureInfoMap &vertexFeatureInfoMap, std::bernoulli_distribution &coinFlip, std::mt19937 &generator,
    const std::string &interactionType, const std::string &trainingOutputFile, const SupportVectorMachine::DoubleVector &eventFeatureList,
    const float maxRadius, const bool useRPhi) const
{
    const Vertex *pBestVertex(nullptr);
    float bestVertexDr(std::numeric_limits<float>::max());

    SupportVectorMachine::DoubleVector bestVertexFeatureList;
    this->GetBestVertex(vertexVector, pBestVertex, bestVertexDr);

    VertexFeatureInfo bestVertexFeatureInfo(vertexFeatureInfoMap.at(pBestVertex));
    this->AddVertexFeaturesToVector(bestVertexFeatureInfo, bestVertexFeatureList, useRPhi);

    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestVertex)
            continue;

        SupportVectorMachine::DoubleVector featureList;
        VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
        this->AddVertexFeaturesToVector(vertexFeatureInfo, featureList, useRPhi);

        if (pBestVertex && (bestVertexDr < maxRadius))
        {
            if (coinFlip(generator))
            {
                LArSvmHelper::ProduceTrainingExample(trainingOutputFile + "_" + interactionType + ".txt", true, eventFeatureList,
                    bestVertexFeatureList, featureList);
            }

            else
            {
                LArSvmHelper::ProduceTrainingExample(trainingOutputFile + "_" + interactionType + ".txt", false, eventFeatureList, featureList,
                    bestVertexFeatureList);
            }
        }
    }

    return pBestVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::GetBestVertex(const VertexVector &vertexVector, const Vertex *&pBestVertex, float &bestVertexDr) const
{
    // Extract input collections
    const MCParticleList *pMCParticleList(nullptr);
    PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);

    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);

    for (const Vertex *const pVertex : vertexVector)
    {
        float mcVertexDr(std::numeric_limits<float>::max());
        for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
        {
            const CartesianVector mcNeutrinoPosition(pMCNeutrino->GetEndpoint().GetX() + m_mcVertexXCorrection, pMCNeutrino->GetEndpoint().GetY(),
                pMCNeutrino->GetEndpoint().GetZ());

            const float dr = (mcNeutrinoPosition - pVertex->GetPosition()).GetMagnitude();
            if (dr < mcVertexDr)
                mcVertexDr = dr;
        }

        if (mcVertexDr < bestVertexDr)
        {
            bestVertexDr = mcVertexDr;
            pBestVertex = pVertex;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::AddVertexFeaturesToVector(const VertexFeatureInfo &vertexFeatureInfo,
    SupportVectorMachine::DoubleVector &featureVector, const bool useRPhi) const
{
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_beamDeweighting));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_energyKick));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_globalAsymmetry));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_localAsymmetry));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_showerAsymmetry));

    if (useRPhi)
        featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_rPhiFeature));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Vertex * SvmVertexSelectionAlgorithm::CompareVertices(const VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
    const SupportVectorMachine::DoubleVector &eventFeatureList, const SupportVectorMachine &supportVectorMachine, const bool useRPhi) const
{
    const Vertex *pBestVertex(vertexVector.front());
    SupportVectorMachine::DoubleVector chosenFeatureList;

    VertexFeatureInfo chosenVertexFeatureInfo(vertexFeatureInfoMap.at(pBestVertex));
    this->AddVertexFeaturesToVector(chosenVertexFeatureInfo, chosenFeatureList, useRPhi);

    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestVertex)
            continue;

        SupportVectorMachine::DoubleVector featureList;
        VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
        this->AddVertexFeaturesToVector(vertexFeatureInfo, featureList, useRPhi);

        if (LArSvmHelper::Classify(supportVectorMachine, eventFeatureList, featureList, chosenFeatureList))
        {
            pBestVertex = pVertex;
            chosenFeatureList = featureList;
        }
    }

    return pBestVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SvmVertexSelectionAlgorithm::PopulateFinalVertexScoreList(const VertexFeatureInfoMap &vertexFeatureInfoMap, const Vertex *const pFavouriteVertex,
    const VertexVector &vertexVector, VertexScoreList &finalVertexScoreList) const
{
    if (pFavouriteVertex)
    {
        const CartesianVector vertexPosition(pFavouriteVertex->GetPosition());

        for (const Vertex *const pVertex : vertexVector)
        {
            if ((pVertex->GetPosition() - vertexPosition).GetMagnitude() < m_rPhiFineTuningRadius)
            {
                const float rPhiScore(vertexFeatureInfoMap.at(pVertex).m_rPhiFeature);
                finalVertexScoreList.emplace_back(pVertex, rPhiScore);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SvmVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArSvmHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SvmFileName", m_svmFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RegionSvmName", m_regionSvmName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexSvmName", m_vertexSvmName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AllowClassifyDuringTraining", m_allowClassifyDuringTraining));

    if ((!m_trainingSetMode || m_allowClassifyDuringTraining))
    {
        if (m_svmFileName.empty() || m_regionSvmName.empty() || m_vertexSvmName.empty())
        {
            std::cout << "SvmVertexSelectionAlgorithm: SvmFileName, RegionSvmName and VertexSvmName must be set if training set mode is" <<
                         "off or we allow classification during training" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        m_svMachineRegion.Initialize(m_svmFileName, m_regionSvmName);
        m_svMachineVertex.Initialize(m_svmFileName, m_vertexSvmName);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCVertexXCorrection", m_mcVertexXCorrection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForGoodView", m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodHits", m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodViews", m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingOutputFileRegion", m_trainingOutputFileRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingOutputFileVertex", m_trainingOutputFileVertex));

    if (m_trainingSetMode && (m_trainingOutputFileRegion.empty() || m_trainingOutputFileVertex.empty()))
    {
        std::cout << "SvmVertexSelectionAlgorithm: TrainingOutputFileRegion and TrainingOutputFileVertex are required for training set " <<
                     "mode" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListName", m_caloHitListName));

    if (m_trainingSetMode && (m_mcParticleListName.empty() || m_caloHitListName.empty()))
    {
        std::cout << "SvmVertexSelectionAlgorithm: MCParticleListName and CaloHitListName are required for training set mode" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerSpineLength", m_minShowerSpineLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamDeweightingConstant", m_beamDeweightingConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LocalAsymmetryConstant", m_localAsymmetryConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GlobalAsymmetryConstant", m_globalAsymmetryConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerAsymmetryConstant", m_showerAsymmetryConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyKickConstant", m_energyKickConstant));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerClusteringDistance", m_showerClusteringDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerClusterHits", m_minShowerClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseShowerClusteringApproximation", m_useShowerClusteringApproximation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RegionRadius", m_regionRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RPhiFineTuningRadius", m_rPhiFineTuningRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrueVertexRadius", m_maxTrueVertexRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseRPhiFeatureForRegion", m_useRPhiFeatureForRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DropFailedRPhiFastScoreCandidates", m_dropFailedRPhiFastScoreCandidates));

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}
} // namespace lar_content
