/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPNeutrinoHierarchyAlgorithm.cc
 *
 *  @brief  Implementation of the MLP neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoradlcontent/LArCheating/MLPCheatHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPLaterTierHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPNeutrinoHierarchyAlgorithm.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MLPNeutrinoHierarchyAlgorithm::MLPNeutrinoHierarchyAlgorithm() :
    m_eventID(0),
    m_trainingMode(false),
    m_trainingFileName("MLPHierarchyTrainingFile.root"),
    m_eventTreeName("EventTree"),
    m_primaryTrackTreeName("PrimaryTrackTree"),
    m_primaryShowerTreeName("PrimaryShowerTree"),
    m_laterTierTrackTrackTreeName("LaterTierTrackTrackTree"),
    m_laterTierTrackShowerTreeName("LaterTierTrackShowerTree"),
    m_mcParticleListName("Input"),
    m_trainingVertexAccuracy(5.f),    
    m_bogusFloat(-999.f),
    m_minClusterSize(5),
    m_slidingFitWindow(20),
    m_regionForDirFit(25.f),
    m_nAngularBins(180),
    m_primaryRegion(15.f),
    m_primaryThresholdTrackPass1(0.8f),
    m_primaryThresholdShowerPass1(0.45f),
    m_laterTierThresholdTrackPass1(0.8f),
    m_laterTierThresholdShowerPass1(0.8f),
    m_primaryThresholdTrackPass2(0.9f),
    m_primaryThresholdShowerPass2(0.9f),
    m_laterTierThresholdTrackPass2(0.0f),
    m_laterTierThresholdShowerPass2(0.0f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MLPNeutrinoHierarchyAlgorithm::~MLPNeutrinoHierarchyAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_eventTreeName.c_str(), m_trainingFileName.c_str(), "UPDATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_primaryTrackTreeName.c_str(), m_trainingFileName.c_str(), "UPDATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_primaryShowerTreeName.c_str(), m_trainingFileName.c_str(), "UPDATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_laterTierTrackTrackTreeName.c_str(), m_trainingFileName.c_str(), "UPDATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_laterTierTrackShowerTreeName.c_str(), m_trainingFileName.c_str(), "UPDATE"));
        }
        catch (...) {}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPNeutrinoHierarchyAlgorithm::Run()
{
    // Get neutrino pfo
    const ParticleFlowObject *pNeutrinoPfo(nullptr);
    if (!this->GetNeutrinoPfo(pNeutrinoPfo))
        return STATUS_CODE_SUCCESS;

    // Ensure neutrino has a vertex
    if (pNeutrinoPfo->GetVertexList().empty())
        return STATUS_CODE_NOT_INITIALIZED;

    // Fill the track/shower vectors
    HierarchyPfoMap trackPfos, showerPfos;
    this->FillTrackShowerVectors(pNeutrinoPfo, trackPfos, showerPfos);

    #ifdef MONITORING
    
    if (m_trainingMode)
    {
        // Do ID stuff
        ++m_eventID;

        std::map<const pandora::ParticleFlowObject *, int> particleIDMap;
        this->GetParticleIDMap(trackPfos, showerPfos, particleIDMap);
        
        // Get truth
        PfoToMCParticleMap matchingMap; ChildToParentPfoMap childToParentPfoMap;
        m_cheatHierarchyTool->FillHierarchyMap(this, matchingMap, childToParentPfoMap);

        // Fill trees
        int nPrimaryTrackLinks(0), nPrimaryShowerLinks(0);
        int nLaterTierTrackTrackLinks(0), nLaterTierTrackShowerLinks(0);
        
        this->FillPrimaryTrees(matchingMap, childToParentPfoMap, pNeutrinoPfo, trackPfos, showerPfos, particleIDMap,
            nPrimaryTrackLinks, nPrimaryShowerLinks); 
        this->FillLaterTierTrees(matchingMap, childToParentPfoMap, pNeutrinoPfo, trackPfos, showerPfos, particleIDMap,
            nLaterTierTrackTrackLinks, nLaterTierTrackShowerLinks);
        this->FillEventTree(trackPfos, showerPfos, nPrimaryTrackLinks, nPrimaryShowerLinks,
            nLaterTierTrackTrackLinks, nLaterTierTrackShowerLinks);        

        return STATUS_CODE_SUCCESS;
    }

    #endif

    // Calculate primary scores
    this->SetPrimaryScores(pNeutrinoPfo, trackPfos, showerPfos);

    // Build initial primary tier
    Hierarchy hierarchy({PfoVector()});
    this->UpdateHierarchy(pNeutrinoPfo, true, true, m_primaryThresholdTrackPass1, m_primaryThresholdShowerPass1, 
        true, trackPfos, showerPfos, hierarchy);

    // If we didn't find any primaries
    if (hierarchy.at(0).empty())
    {
        // Set everything as primary and leave.
        this->BuildPandoraHierarchy(pNeutrinoPfo, trackPfos, showerPfos);

        return STATUS_CODE_SUCCESS;
    }

    // Set later tier scores
    this->SetLaterTierScores(pNeutrinoPfo, trackPfos, showerPfos);

    // Build the later tier
    this->UpdateHierarchy(pNeutrinoPfo, false, false, m_laterTierThresholdTrackPass1, m_laterTierThresholdShowerPass1, 
        true, trackPfos, showerPfos, hierarchy);

    // Try to recover primaries using laterTierScore
    this->UpdateHierarchy(pNeutrinoPfo, true, false, m_primaryThresholdTrackPass2, m_primaryThresholdShowerPass2, 
        false, trackPfos, showerPfos, hierarchy);

    // Try to recover any children using laterTierScore
    this->UpdateHierarchy(pNeutrinoPfo, false, false, m_laterTierThresholdTrackPass2, m_laterTierThresholdShowerPass2, 
        true, trackPfos, showerPfos, hierarchy);

    this->BuildPandoraHierarchy(pNeutrinoPfo, trackPfos, showerPfos);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPNeutrinoHierarchyAlgorithm::GetNeutrinoPfo(const ParticleFlowObject *&pNeutrinoPfo) const
{
    const PfoList *pPfoList = nullptr;

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, 
        PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
        return false;

    // ATTN Enforces that only one pfo, of neutrino-type, be in the specified input list
    pNeutrinoPfo = (1 == pPfoList->size()) ? *(pPfoList->begin()) : nullptr;

    if (!pNeutrinoPfo || !LArPfoHelper::IsNeutrino(pNeutrinoPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::FillTrackShowerVectors(const ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, 
    HierarchyPfoMap &showerPfos) const
{
    // Sliding fit shenanigans
    const float pitchU{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U)};
    const float pitchV{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V)};
    const float pitchW{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W)};
    const float slidingFitPitch(std::max({pitchU, pitchV, pitchW}));

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject * pPfo : *pPfoList)
        {
            // Apply hit cut
            if (this->GetNSpacepoints(pPfo) < m_minClusterSize)
                continue;

            // Attempt sliding linear fit
            try
            {
                ClusterList clusters3D;
                LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

                if (clusters3D.size() != 1)
                    continue;

                const ThreeDSlidingFitResult slidingFitResult(*clusters3D.begin(), m_slidingFitWindow, slidingFitPitch);

                // We need directions...
                CartesianVector upstreamVertex(m_bogusFloat, m_bogusFloat, m_bogusFloat), upstreamDirection(m_bogusFloat, m_bogusFloat, m_bogusFloat), 
                    downstreamVertex(m_bogusFloat, m_bogusFloat, m_bogusFloat), downstreamDirection(m_bogusFloat, m_bogusFloat, m_bogusFloat);

                if (!this->GetExtremalVerticesAndDirections(pNeutrinoPfo, pPfo, slidingFitResult, upstreamVertex, upstreamDirection, downstreamVertex, downstreamDirection))
                    continue;

                // Create track/shower objects
                if (pPfo->GetParticleId() == 13)
                {
                    trackPfos.insert(std::make_pair(pPfo, HierarchyPfo(true, pPfo, slidingFitResult, upstreamVertex, upstreamDirection, downstreamVertex, downstreamDirection)));
                }
                else if (pPfo->GetParticleId() == 11) 
                {
                    showerPfos.insert(std::make_pair(pPfo, HierarchyPfo(false, pPfo, slidingFitResult, upstreamVertex, upstreamDirection, downstreamVertex, downstreamDirection)));
                }
            }
            catch (...) { continue; }
        }
    }

    std::cout << "We have " << trackPfos.size() << " track(s) and " << showerPfos.size() << " shower(s)" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetNSpacepoints(const ParticleFlowObject *const pPfo) const
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    int total3DHits(0);

    for (const Cluster *const pCluster3D : clusterList3D)
        total3DHits += pCluster3D->GetNCaloHits();

    return total3DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetSeparation(const ParticleFlowObject *const pParentPfo, const ParticleFlowObject *const pChildPfo) const
{
    CartesianPointVector parentPositions3D, childPositions3D;
    LArPfoHelper::GetCoordinateVector(pParentPfo, TPC_3D, parentPositions3D);
    LArPfoHelper::GetCoordinateVector(pChildPfo, TPC_3D, childPositions3D);

    float sepSq(std::numeric_limits<float>::max());

    for (const CartesianVector &parentPosition3D : parentPositions3D)
    {
        for (const CartesianVector &childPosition3D : childPositions3D)
        {
            const float thisSepSq((parentPosition3D - childPosition3D).GetMagnitudeSquared());

            sepSq = std::min(thisSepSq, sepSq);
        }
    }

    return std::sqrt(sepSq);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPNeutrinoHierarchyAlgorithm::GetExtremalVerticesAndDirections(const ParticleFlowObject *const pNeutrinoPfo, const ParticleFlowObject *const pPfo, 
    const ThreeDSlidingFitResult &slidingFitResult, CartesianVector &upstreamVertex, CartesianVector &upstreamDirection, CartesianVector &downstreamVertex, 
    CartesianVector &downstreamDirection) const
{
    // First, get the neutrino vertex
    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector nuVertex(pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

    if (pPfo->GetParticleId() == 13)
    {
        try
        {
            const LArPointingCluster pointingCluster(slidingFitResult); // this could throw an exception
            const CartesianVector innerVertex(pointingCluster.GetInnerVertex().GetPosition().GetX(), pointingCluster.GetInnerVertex().GetPosition().GetY(), pointingCluster.GetInnerVertex().GetPosition().GetZ()); 
            CartesianVector innerDirection(pointingCluster.GetInnerVertex().GetDirection().GetX(), pointingCluster.GetInnerVertex().GetDirection().GetY(), pointingCluster.GetInnerVertex().GetDirection().GetZ()); 
            const CartesianVector outerVertex(pointingCluster.GetOuterVertex().GetPosition().GetX(), pointingCluster.GetOuterVertex().GetPosition().GetY(), pointingCluster.GetOuterVertex().GetPosition().GetZ());
            CartesianVector outerDirection(pointingCluster.GetOuterVertex().GetDirection().GetX(), pointingCluster.GetOuterVertex().GetDirection().GetY(), pointingCluster.GetOuterVertex().GetDirection().GetZ());

            // want them to both point into the particle
            if (innerDirection.GetDotProduct(outerVertex - innerVertex) < 0.f)
                innerDirection *= (-1.f);

            if (outerDirection.GetDotProduct(innerVertex - outerVertex) < 0.f)
                outerDirection *= (-1.f);

            // determine upstream/downstream
            if ((innerVertex - nuVertex).GetMagnitudeSquared() < (outerVertex - nuVertex).GetMagnitudeSquared())
            {
                upstreamVertex = innerVertex;
                upstreamDirection = innerDirection;
                downstreamVertex = outerVertex;
                downstreamDirection = outerDirection;
            }
            else
            {
                upstreamVertex = outerVertex;
                upstreamDirection = outerDirection;
                downstreamVertex = innerVertex;
                downstreamDirection = innerDirection;
            }
        }
        catch(...) { return false; }        
    }
    else
    {
        const float minSepSq((slidingFitResult.GetGlobalMinLayerPosition() - nuVertex).GetMagnitudeSquared());
        const float maxSepSq((slidingFitResult.GetGlobalMaxLayerPosition() - nuVertex).GetMagnitudeSquared());

        upstreamVertex = (minSepSq < maxSepSq) ? slidingFitResult.GetGlobalMinLayerPosition() : slidingFitResult.GetGlobalMaxLayerPosition();
        downstreamVertex = (minSepSq < maxSepSq) ? slidingFitResult.GetGlobalMaxLayerPosition() : slidingFitResult.GetGlobalMinLayerPosition();

        // find directions, demand that we have at least one
        const bool upstreamDirectionSet(this->GetShowerDirection(pPfo, upstreamVertex, upstreamDirection));
        const bool downstreamDirectionSet(this->GetShowerDirection(pPfo, downstreamVertex, downstreamDirection));

        if (!upstreamDirectionSet && !downstreamDirectionSet)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPNeutrinoHierarchyAlgorithm::GetShowerDirection(const ParticleFlowObject *const pPfo, const CartesianVector &vertex,
    CartesianVector &direction) const
{
    CartesianPointVector pointVector;
    LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, pointVector);

    if (pointVector.empty())
        return false;

    // set direction by angular distribution
    const float angleMin(0.f), angleMax(2.f * M_PI);
    const float binWidth((angleMax - angleMin) / static_cast<float>(m_nAngularBins));

    std::vector<std::vector<int>> spatialDist(m_nAngularBins, std::vector<int>(m_nAngularBins, 0));
    std::vector<std::vector<float>> tieBreakerDist(m_nAngularBins, std::vector<float>(m_nAngularBins, 0.f));

    int highestSP(0), bestTheta0YZBin(-1), bestTheta0XZBin(-1);
    float tieBreaker(std::numeric_limits<float>::max());

    for (const CartesianVector &position3D : pointVector)
    {
        const CartesianVector displacement(position3D - vertex);
        const float mag = displacement.GetMagnitude();

        if (mag > m_regionForDirFit)
            continue;

        const float magXZ = sqrt((displacement.GetX() * displacement.GetX()) + (displacement.GetZ() * displacement.GetZ()));

        float theta0YZ = (mag < std::numeric_limits<float>::epsilon()) ? 0.f : 
            (std::fabs(std::fabs(displacement.GetY() / mag) - 1.f) < std::numeric_limits<float>::epsilon()) ? 0.f : 
            std::acos(displacement.GetY() / mag);

        float theta0XZ = (magXZ < std::numeric_limits<float>::epsilon()) ? 0.f : 
            (std::fabs(std::fabs(displacement.GetX() / magXZ) - 1.f) < std::numeric_limits<float>::epsilon()) ? 0.f :
            std::acos(displacement.GetX() / magXZ);

        // try do signed-ness
        if (displacement.GetZ() < 0.f)
            theta0XZ = (2.0 * M_PI) - theta0XZ;

        const int bin0YZ = std::floor(theta0YZ / binWidth);
        const int bin0XZ = std::floor(theta0XZ / binWidth);

        spatialDist[bin0YZ][bin0XZ] += 1;
        tieBreakerDist[bin0YZ][bin0XZ] += (1.f / mag); // tie-breaker

        if (((spatialDist[bin0YZ][bin0XZ] == highestSP) && (tieBreakerDist[bin0YZ][bin0XZ] < tieBreaker)) ||
            (spatialDist[bin0YZ][bin0XZ] > highestSP))
        {
            highestSP = spatialDist[bin0YZ][bin0XZ];
            tieBreaker = tieBreakerDist[bin0YZ][bin0XZ];
            bestTheta0YZBin = bin0YZ;
            bestTheta0XZBin = bin0XZ;
        }
    }

    if ((bestTheta0YZBin < 0) || (bestTheta0XZBin < 0))
        return false;

    const float bestTheta0YZ = angleMin + ((static_cast<float>(bestTheta0YZBin) + 0.5f) * binWidth);
    const float bestTheta0XZ = angleMin + ((static_cast<float>(bestTheta0XZBin) + 0.5f) * binWidth);

    direction = CartesianVector(std::sin(bestTheta0YZ) * std::cos(bestTheta0XZ), std::cos(bestTheta0YZ), std::sin(bestTheta0YZ) * std::sin(bestTheta0XZ));

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::SetPrimaryScores(const ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, 
    HierarchyPfoMap &showerPfos) const
{
    for (auto& [pPfo, hierarchyPfo] : trackPfos)
        hierarchyPfo.SetPrimaryScore(this->GetPrimaryScore(pNeutrinoPfo, trackPfos, hierarchyPfo));

    for (auto& [pPfo, hierarchyPfo] : showerPfos)
        hierarchyPfo.SetPrimaryScore(this->GetPrimaryScore(pNeutrinoPfo, trackPfos, hierarchyPfo));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetPrimaryScore(const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, 
    const HierarchyPfo &hierarchyPfo) const
{
    std::vector<MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams> networkParamVector;
    float primaryScore(m_bogusFloat);

    if (m_primaryHierarchyTool->Run(this, pNeutrinoPfo, trackPfos, hierarchyPfo, networkParamVector, primaryScore) != STATUS_CODE_SUCCESS)
        return m_bogusFloat;

    return primaryScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::UpdateHierarchy(const ParticleFlowObject *const pNeutrinoPfo, const bool buildPrimaryTier, 
    const bool usePrimaryScore, const float trackThreshold, const float showerThreshold, const bool isLowerThreshold, 
    HierarchyPfoMap &trackPfos, HierarchyPfoMap &showerPfos, Hierarchy &hierarchy) const
{
    bool found(true);
    unsigned int tierIndex(buildPrimaryTier ? 0 : 1);

    while (found)
    {
        found = false;

        std::vector<const ParticleFlowObject*> &thisTier(hierarchy.at(tierIndex));

        for (bool isTrack : {true, false})
        {
            std::map<const ParticleFlowObject*, HierarchyPfo> &hierarchyPfoMap(isTrack ? trackPfos : showerPfos);

            for (auto& [pChildPfo, childPfo] : hierarchyPfoMap)
            {
                // Continue if parent already found
                if (childPfo.GetIsInHierarchy())
                    continue;

                const ParticleFlowObject *pPredictedParent(pNeutrinoPfo);

                // If not buildPrimaryTier correct pPredictedParent
                if (!buildPrimaryTier)
                {
                    // Do we have a predicted parent
                    if (!childPfo.GetPredictedParentPfo())
                        continue;

                    // Is predicted parent in preceeding tier?
                    std::vector<const ParticleFlowObject*> &preceedingTier(hierarchy.at(tierIndex - 1));
                    pPredictedParent = childPfo.GetPredictedParentPfo();
                    if (std::find(preceedingTier.begin(), preceedingTier.end(), pPredictedParent) == preceedingTier.end())
                        continue;
                }

                // Does it pass tier cut?
                const float networkScore(usePrimaryScore ? childPfo.GetPrimaryScore() : childPfo.GetLaterTierScore());
                const float thresholdScore(isTrack ? trackThreshold : showerThreshold);

                if ((isLowerThreshold && (networkScore > thresholdScore)) || (!isLowerThreshold && (networkScore < thresholdScore)))
                {
                    found = true;

                    // Add child to hierarchy tier
                    thisTier.emplace_back(pChildPfo);

                    // Add info to childPfo
                    childPfo.SetIsInHierarchy(true);
                    childPfo.SetParentPfo(pPredictedParent);

                    // Add info to parentPfo
                    if (!buildPrimaryTier)
                    {
                        if (trackPfos.find(pPredictedParent) == trackPfos.end())
                            throw StatusCodeException(STATUS_CODE_FAILURE);

                        trackPfos.at(pPredictedParent).AddChildPfo(pChildPfo);
                    }
                }
            }
        }

        // Make sure to add the next tier
        if (hierarchy.size() == (tierIndex + 1))
            hierarchy.push_back(PfoVector());

        // If buildPrimaryTier, our work is done!
        if (buildPrimaryTier)
            break;

        ++tierIndex;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::SetLaterTierScores(const ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoMap &trackPfos, 
    HierarchyPfoMap &showerPfos) const
{
    for (bool isTrack : {true, false})
    {
        std::map<const ParticleFlowObject*, HierarchyPfo> &hierarchyPfoMap(isTrack ? trackPfos : showerPfos);

        for (auto& [pChildPfo, childPfo] : hierarchyPfoMap)
        {
            // Continue if parent already found
            if (childPfo.GetIsInHierarchy())
                continue;

            // Only consider for later tier if far away from nu vertex
            const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
            const CartesianVector nuVertex(pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

            if ((childPfo.GetUpstreamVertex() - nuVertex).GetMagnitudeSquared() < (m_primaryRegion * m_primaryRegion))
                continue;

            float highestScore(0.f);
            int highestNHits(0);

            for (auto& [pParentPfo, parentPfo] : trackPfos)
            {
                if (pChildPfo == pParentPfo)
                    continue;

                const float thisScore(this->GetLaterTierScore(pNeutrinoPfo, parentPfo, childPfo));

                const bool best((std::fabs(thisScore - highestScore) < std::numeric_limits<float>::epsilon()) ?
                                (this->GetSeparation(pParentPfo, pChildPfo) > highestNHits) : (thisScore > highestScore));
                
                if (best)
                {
                    highestScore = thisScore;
                    highestNHits = this->GetSeparation(pParentPfo, pChildPfo);
                    childPfo.SetLaterTierScore(thisScore);
                    childPfo.SetPredictedParentPfo(pParentPfo);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetLaterTierScore(const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfo &parentPfo, 
    const HierarchyPfo &childPfo) const
{
    std::vector<MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams> networkParamVector;
    float laterTierScore(m_bogusFloat);

    if (m_laterTierHierarchyTool->Run(this, pNeutrinoPfo, parentPfo, childPfo, networkParamVector, laterTierScore) != STATUS_CODE_SUCCESS)
        return m_bogusFloat;

    return laterTierScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::BuildPandoraHierarchy(const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, 
    const HierarchyPfoMap &showerPfos) const
{
    PfoVector pfoVector;

    // Get All pfos
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
            pfoVector.emplace_back(pPfo);
    }

    // Sort to maintain reproducibility
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);

    // Build known hierarchy
    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        const HierarchyPfo* hierarchyPfo(nullptr);

        if (trackPfos.find(pPfo) != trackPfos.end())
        {
            hierarchyPfo = &trackPfos.at(pPfo);
        }
        else if (showerPfos.find(pPfo) != showerPfos.end())
        {
            hierarchyPfo = &showerPfos.at(pPfo);
        }
        else
        {
            // If alg never handled the pfo, assign as primary
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pPfo));
            continue;
        }

        // If pfo was never assigned to hierarchy, add as primary
        if (!hierarchyPfo->GetIsInHierarchy())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pPfo));
            continue;
        }

        // If parent is a neutrino we have to assign its parent
        if (hierarchyPfo->GetParentPfo() == pNeutrinoPfo)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pPfo));

        // Assign parent-child links of its children
        PfoVector childPfos(hierarchyPfo->GetChildPfoVector());
        std::sort(childPfos.begin(), childPfos.end(), LArPfoHelper::SortByNHits);

        if ((!childPfos.empty()) && (showerPfos.find(pPfo) != showerPfos.end()))
        {
            std::cout << "ISOBEL YOU MORON, SHOWER WITH CHILDREN" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        for (const ParticleFlowObject *const pChildPfo : childPfos)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfo, pChildPfo));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Training functions
//------------------------------------------------------------------------------------------------------------------------------------------

#ifdef MONITORING

bool MLPNeutrinoHierarchyAlgorithm::ShouldTrainOnEvent(const ParticleFlowObject *const pNeutrinoPfo) const
{
    const MCParticleList *pMCParticleList(nullptr);
    if (PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList) != STATUS_CODE_SUCCESS)
        return false;

    if (!pMCParticleList || pMCParticleList->empty())
        return false;

    // Apply nu vertex accuracy requirement
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);

    if (mcNeutrinoVector.size() != 1)
        return false;

    const CartesianVector &trueNuVertex(mcNeutrinoVector.front()->GetVertex());
    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector recoNuVertex(pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

    if ((trueNuVertex - recoNuVertex).GetMagnitudeSquared() > (m_trainingVertexAccuracy * m_trainingVertexAccuracy))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::GetParticleIDMap(const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos, 
    std::map<const ParticleFlowObject *, int> &particleIDMap)
{
    PfoVector pfoVector;

    for (const auto& [pPfo, hierarchyPfo] : trackPfos)
        pfoVector.push_back(pPfo);

    for (const auto& [pPfo, hierarchyPfo] : showerPfos)
        pfoVector.push_back(pPfo);

    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject * pPfo : pfoVector)
        particleIDMap[pPfo] = particleIDMap.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::FillEventTree(const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos,
    const int nPrimaryTrackLinks, const int nPrimaryShowerLinks, const int nLaterTierTrackTrackLinks, const int nLaterTierTrackShowerLinks) const
{
    const int nParticles(trackPfos.size() + showerPfos.size());

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "NRecoParticles", nParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "EventID", m_eventID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "NPrimaryTrackLinks", nPrimaryTrackLinks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "NPrimaryShowerLinks", nPrimaryShowerLinks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "NLaterTierTrackTrackLinks", nLaterTierTrackTrackLinks));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "NLaterTierTrackShowerLinks", nLaterTierTrackShowerLinks));    
    
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_eventTreeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::FillPrimaryTrees(const PfoToMCParticleMap &matchingMap, const ChildToParentPfoMap &childToParentPfoMap,
    const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos, 
    const std::map<const ParticleFlowObject *, int> &particleIDMap, int &nPrimaryTrackLinks, int &nPrimaryShowerLinks) const
{
    const bool isTrainingEvent(this->ShouldTrainOnEvent(pNeutrinoPfo));
  
    for (const bool isTrack : {true, false})
    {
        const std::string treeName(isTrack ? m_primaryTrackTreeName : m_primaryShowerTreeName);      
        const HierarchyPfoMap &hierarchyPfoMap(isTrack ? trackPfos : showerPfos);

        for (const auto& [pPfo, hierarchyPfo] : hierarchyPfoMap)
        {
            // Get particle ID
            if (particleIDMap.find(pPfo) == particleIDMap.end())
                continue;

            const int particleID(particleIDMap.at(pPfo));

	    // Is training link?
            bool isNeutronChild(false);
            bool isTrainingLink(isTrainingEvent && !((m_cheatHierarchyTool->IsNeutronChild(pPfo, matchingMap, isNeutronChild) != STATUS_CODE_SUCCESS) || isNeutronChild));

            // Get the true info
            // Orientation == true if true start point is at upstream position
            bool isTrueLink(false), trueChildOrientation(false);
            if (m_cheatHierarchyTool->Run(matchingMap, childToParentPfoMap, pNeutrinoPfo, hierarchyPfo, 
                isTrueLink, trueChildOrientation) != STATUS_CODE_SUCCESS)
            {
		continue;
            }

            // Get true hierarchy info
            const ParticleFlowObject *const pTrueParentPfo(childToParentPfoMap.at(pPfo).first);
            const int trueVisibleGen(childToParentPfoMap.at(pPfo).second);

            int trueParentID(-1); // default neutrino
            if (pTrueParentPfo != pNeutrinoPfo)
            {
                if (particleIDMap.find(pTrueParentPfo) == particleIDMap.end())
                    continue;

                trueParentID = particleIDMap.at(pTrueParentPfo);
            }

	    // Get reco info
            std::vector<MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams> networkParamVector;
            float primaryScore(m_bogusFloat);

            if (m_primaryHierarchyTool->Run(this, pNeutrinoPfo, trackPfos, hierarchyPfo, networkParamVector, primaryScore) != STATUS_CODE_SUCCESS)
		continue;

            for (const MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams &networkParams : networkParamVector)
            {
                const bool useChildUpstream(std::fabs(networkParams.m_isPOIClosestToNu - 1.f) < std::numeric_limits<float>::epsilon());
                this->FillPrimaryTree(treeName, isTrainingLink, isTrueLink, (useChildUpstream == trueChildOrientation), trueVisibleGen,
                    trueParentID, particleID, hierarchyPfo.GetUpstreamVertex(), hierarchyPfo.GetDownstreamVertex(), networkParams);
                
                isTrack ? ++nPrimaryTrackLinks : ++nPrimaryShowerLinks;                
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::FillPrimaryTree(const std::string &treeName, const bool isTrainingLink, const bool isTrueLink,
    const bool isOrientationCorrect, const int trueVisibleGen, const int trueParentID, const int particleID,
    const CartesianVector &upstreamVertex, const CartesianVector &downstreamVertex,
    const MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams &primaryNetworkParams) const
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EventID", m_eventID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParticleID", particleID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "IsTrainingLink", (isTrainingLink ? 1 : 0)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "IsTrueLink", (isTrueLink ? 1 : 0)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "IsOrientationCorrect", (isOrientationCorrect ? 1 : 0)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TrueVisibleGeneration", trueVisibleGen));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TrueVisibleParentID", trueParentID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpstreamVertexX", upstreamVertex.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpstreamVertexY", upstreamVertex.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpstreamVertexZ", upstreamVertex.GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DownstreamVertexX", downstreamVertex.GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DownstreamVertexY", downstreamVertex.GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DownstreamVertexZ", downstreamVertex.GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "NSpacepoints", primaryNetworkParams.m_nSpacepoints));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "NuSeparation", primaryNetworkParams.m_nuSeparation));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "VertexRegionNHits", primaryNetworkParams.m_vertexRegionNHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "VertexRegionNParticles", primaryNetworkParams.m_vertexRegionNParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DCA", primaryNetworkParams.m_dca));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ConnectionExtrapDistance", primaryNetworkParams.m_connectionExtrapDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "IsPOIClosestToNu", primaryNetworkParams.m_isPOIClosestToNu));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentConnectionDistance", primaryNetworkParams.m_parentConnectionDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildConnectionDistance", primaryNetworkParams.m_childConnectionDistance));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::FillLaterTierTrees(const PfoToMCParticleMap &matchingMap, const ChildToParentPfoMap &childToParentPfoMap,
    const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const HierarchyPfoMap &showerPfos, 
    const std::map<const ParticleFlowObject *, int> &particleIDMap, int &nTrackLinks, int &nShowerLinks) const
{
    const bool isTrainingEvent(this->ShouldTrainOnEvent(pNeutrinoPfo));
  
    for (const auto& [pParentTrackPfo, hierarchyParentTrackPfo] : trackPfos)
    {
        for (const bool isTrack : {true, false})
        {
  	    const std::string treeName(isTrack ? m_laterTierTrackTrackTreeName : m_laterTierTrackShowerTreeName);	  
            const HierarchyPfoMap &hierarchyPfoMap(isTrack ? trackPfos : showerPfos);
 
            for (const auto& [pChildPfo, hierarchyChildPfo] : hierarchyPfoMap)
            {
                if (pParentTrackPfo == pChildPfo)
                    continue;

                // Get particle IDs
                if ((particleIDMap.find(pParentTrackPfo) == particleIDMap.end()) ||
                    (particleIDMap.find(pChildPfo) == particleIDMap.end()))
                {
                    continue;
                }

                const int parentID(particleIDMap.at(pParentTrackPfo));
                const int childID(particleIDMap.at(pChildPfo));

		// Get true primary info
                bool isTruePrimaryLink(false), tempBool(false);
                if (m_cheatHierarchyTool->Run(matchingMap, childToParentPfoMap, pNeutrinoPfo, hierarchyChildPfo, 
                    isTruePrimaryLink, tempBool) != STATUS_CODE_SUCCESS)
                {
                    continue;
                }

		// Do not train on child primaries
		bool isTrainingLink(isTrainingEvent && !isTruePrimaryLink);

                // Get the true info
                // Orientation == true if true start point is at upstream position
                bool isTrueLink(false), trueParentOrientation(false), trueChildOrientation(false);
                if (m_cheatHierarchyTool->Run(matchingMap, childToParentPfoMap, hierarchyParentTrackPfo, hierarchyChildPfo, 
                    isTrueLink, trueParentOrientation, trueChildOrientation) != STATUS_CODE_SUCCESS)
                {
                    continue;
                }

                // Get true hierarchy info
                const int trueVisibleGen(childToParentPfoMap.at(pChildPfo).second);

                // For parent we care about the endpoint orientation
                trueParentOrientation = !trueParentOrientation;

                // Now get MLP info for ALL orientations
                float laterTierScore(m_bogusFloat);
                std::vector<MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams> networkParamVector;

                if (m_laterTierHierarchyTool->Run(this, pNeutrinoPfo, hierarchyParentTrackPfo, hierarchyChildPfo, networkParamVector, laterTierScore) != STATUS_CODE_SUCCESS)
                    continue;

                // Check that we have an edge with correctLinkOrientation
                int nCorrectLinkOrientation(0);
                for (const MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams &networkParams : networkParamVector)
                {
                    const bool useParentUpstream(std::fabs(networkParams.m_parentIsPOIClosestToNu - 1.f) < std::numeric_limits<float>::epsilon());
                    const bool useChildUpstream(std::fabs(networkParams.m_childIsPOIClosestToNu - 1.f) < std::numeric_limits<float>::epsilon());
                    const bool isOrientationCorrect((useParentUpstream == trueParentOrientation) && (useChildUpstream == trueChildOrientation));

                    if (isOrientationCorrect)
                        ++nCorrectLinkOrientation;
                }

		// For training, make sure we have one correct orientation link
                isTrainingLink = (isTrainingLink && (nCorrectLinkOrientation == 1));

                // Get training cut info
		std::pair<float, float> trainingCuts({-999.f, -999.f});
		
		if (isTrainingLink)
		  trainingCuts = this->GetTrainingCuts(hierarchyParentTrackPfo, hierarchyChildPfo, trueParentOrientation, trueChildOrientation);

                for (const MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams &networkParams : networkParamVector)
                {
                    const bool useParentUpstream(std::fabs(networkParams.m_parentIsPOIClosestToNu - 1.f) < std::numeric_limits<float>::epsilon());
                    const bool useChildUpstream(std::fabs(networkParams.m_childIsPOIClosestToNu - 1.f) < std::numeric_limits<float>::epsilon());
                    const bool isOrientationCorrect((useParentUpstream == trueParentOrientation) && (useChildUpstream == trueChildOrientation));
                    
                    this->FillLaterTierTree(treeName, isTrainingLink, isTrueLink, isOrientationCorrect, trueVisibleGen, trainingCuts, parentID, childID,
                        networkParams);
                    
                    isTrack ? ++nTrackLinks : ++nShowerLinks;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<float, float> MLPNeutrinoHierarchyAlgorithm::GetTrainingCuts(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo,
    const bool trueParentOrientation, const bool trueChildOrientation) const
{
    const CartesianVector &parentEnd(trueParentOrientation ? parentHierarchyPfo.GetUpstreamVertex() : 
        parentHierarchyPfo.GetDownstreamVertex());
    const CartesianVector &parentEndDirection(trueParentOrientation ? parentHierarchyPfo.GetUpstreamDirection() : 
        parentHierarchyPfo.GetDownstreamDirection());
    const CartesianVector &childStart(trueChildOrientation ? childHierarchyPfo.GetUpstreamVertex() : 
        childHierarchyPfo.GetDownstreamVertex());
    const CartesianVector &childStartDirection(trueChildOrientation ? childHierarchyPfo.GetUpstreamDirection() : 
        childHierarchyPfo.GetDownstreamDirection());

    std::pair<float, float> trainingCuts(std::make_pair(m_bogusFloat, m_bogusFloat));

    m_primaryHierarchyTool->CalculateConnectionDistances(parentEnd, parentEndDirection, childStart, childStartDirection, 
        trainingCuts.first, trainingCuts.second);

    return trainingCuts;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::FillLaterTierTree(const std::string &treeName, const bool isTrainingLink, const bool isTrueLink,
    const bool isOrientationCorrect, const int childTrueGen, const std::pair<float, float> &trainingCuts, 
    const int parentID, const int childID, const MLPLaterTierHierarchyTool::MLPLaterTierNetworkParams &networkParams) const
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EventID", m_eventID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentParticleID", parentID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildParticleID", childID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "IsTrainingLink", (isTrainingLink ? 1 : 0)));    
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "IsTrueLink", (isTrueLink ? 1 : 0)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "IsOrientationCorrect", (isOrientationCorrect ? 1 : 0)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildTrueVisibleGeneration", childTrueGen));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TrainingCutL", trainingCuts.first));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TrainingCutT", trainingCuts.second));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentTrackScore", networkParams.m_parentTrackScore));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildTrackScore", networkParams.m_childTrackScore));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentNSpacepoints", networkParams.m_parentNSpacepoints));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildNSpacepoints", networkParams.m_childNSpacepoints));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "Separation3D", networkParams.m_separation3D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentNuVertexSep", networkParams.m_parentNuVertexSep));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildNuVertexSep", networkParams.m_childNuVertexSep));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentEndRegionNHits", networkParams.m_parentEndRegionNHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentEndRegionNParticles", networkParams.m_parentEndRegionNParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentEndRegionRToWall", networkParams.m_parentEndRegionRToWall));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "VertexSeparation", networkParams.m_vertexSeparation));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DoesChildConnect", networkParams.m_doesChildConnect));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "OvershootStartDCA", networkParams.m_overshootStartDCA));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "OvershootStartL", networkParams.m_overshootStartL));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "OvershootEndDCA", networkParams.m_overshootEndDCA));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "OvershootEndL", networkParams.m_overshootEndL));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildCPDCA", networkParams.m_childCPDCA));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildCPExtrapDistance", networkParams.m_childCPExtrapDistance));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildCPLRatio", networkParams.m_childCPLRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentCPNUpstreamHits", networkParams.m_parentCPNUpstreamHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentCPNDownstreamHits", networkParams.m_parentCPNDownstreamHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentCPNHitRatio", networkParams.m_parentCPNHitRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentCPEigenvalueRatio", networkParams.m_parentCPEigenvalueRatio));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentCPOpeningAngle", networkParams.m_parentCPOpeningAngle));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ParentIsPOIClosestToNu", networkParams.m_parentIsPOIClosestToNu));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChildIsPOIClosestToNu", networkParams.m_childIsPOIClosestToNu));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName.c_str()));
}

#endif

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPNeutrinoHierarchyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventID", m_eventID));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingVertexAccuracy", m_trainingVertexAccuracy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventTreeName", m_eventTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryTrackTreeName", m_primaryTrackTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryShowerTreeName", m_primaryShowerTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierTrackTrackTreeName", m_laterTierTrackTrackTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierTrackShowerTreeName", m_laterTierTrackShowerTreeName));
   
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BogusFloat", m_bogusFloat));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterSize", m_minClusterSize));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RegionForDirectionFit", m_regionForDirFit));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NAngularBins", m_nAngularBins));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryRegion", m_primaryRegion));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryThresholdTrackPass1", m_primaryThresholdTrackPass1));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryThresholdShowerPass1", m_primaryThresholdShowerPass1));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierThresholdTrackPass1", m_laterTierThresholdTrackPass1));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierThresholdShowerPass1", m_laterTierThresholdShowerPass1));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryThresholdTrackPass2", m_primaryThresholdTrackPass2));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryThresholdShowerPass2", m_primaryThresholdShowerPass2));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierThresholdTrackPass2", m_laterTierThresholdTrackPass2));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierThresholdShowerPass2", m_laterTierThresholdShowerPass2));

    AlgorithmTool *pAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "MLPPrimaryHierarchyTool", pAlgorithmTool));
    m_primaryHierarchyTool = dynamic_cast<MLPPrimaryHierarchyTool *>(pAlgorithmTool);

    if (!m_primaryHierarchyTool)
        return STATUS_CODE_INVALID_PARAMETER;

    pAlgorithmTool = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "MLPLaterTierHierarchyTool", pAlgorithmTool));
    m_laterTierHierarchyTool = dynamic_cast<MLPLaterTierHierarchyTool *>(pAlgorithmTool);

    if (!m_laterTierHierarchyTool)
        return STATUS_CODE_INVALID_PARAMETER;

    if (m_trainingMode)
    {
        pAlgorithmTool = nullptr;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "MLPCheatHierarchyTool", pAlgorithmTool));
        m_cheatHierarchyTool = dynamic_cast<MLPCheatHierarchyTool *>(pAlgorithmTool);

        if (!m_cheatHierarchyTool)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
