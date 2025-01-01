/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPNeutrinoHierarchyAlgorithm.cc
 *
 *  @brief  Implementation of the MLP neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

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

StatusCode MLPNeutrinoHierarchyAlgorithm::Run()
{
    // Get neutrino pfo
    const ParticleFlowObject *pNeutrinoPfo(nullptr);
    if (!this->GetNeutrinoPfo(pNeutrinoPfo))
        return STATUS_CODE_SUCCESS;

    // Ensure neutrino has a vertex
    if (pNeutrinoPfo->GetVertexList().empty())
        return STATUS_CODE_NOT_INITIALIZED;

    // Give PFPs IDs to keep track of them
    this->DetermineIsobelID();

    // Fill the track/shower vectors
    HierarchyPfoMap trackPfos, showerPfos;
    this->FillTrackShowerVectors(pNeutrinoPfo, trackPfos, showerPfos);

    // Calculate primary scores
    this->SetPrimaryScores(pNeutrinoPfo, trackPfos, showerPfos);

    //////////////////////////////////////////////////////////////////
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "PRIMARY SCORES" << std::endl;
    for (const auto& [pPfo, hierarchyPfo] : trackPfos)
        std::cout << m_isobelID[pPfo] << ": " << hierarchyPfo.GetPrimaryScore() << std::endl;

    for (const auto& [pPfo, hierarchyPfo] : showerPfos)
        std::cout << m_isobelID[pPfo] << ": " << hierarchyPfo.GetPrimaryScore() << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    //////////////////////////////////////////////////////////////////

    // Build initial primary tier
    Hierarchy hierarchy({PfoVector()});
    this->UpdateHierarchy(pNeutrinoPfo, true, true, m_primaryThresholdTrackPass1, m_primaryThresholdShowerPass1, 
        true, trackPfos, showerPfos, hierarchy);

    if (hierarchy.empty())
    {
        // Set everything as primary and leave.
        this->BuildPandoraHierarchy(pNeutrinoPfo, trackPfos, showerPfos);

        return STATUS_CODE_SUCCESS;
    }

    //////////////////////////////////////////////////////////////////    
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "PRIMARY TIER:" << std::endl;
    for (const ParticleFlowObject* pPfo : hierarchy.at(0))
        std::cout << m_isobelID[pPfo] << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    //////////////////////////////////////////////////////////////////    

    // Set later tier scores
    this->SetLaterTierScores(pNeutrinoPfo, trackPfos, showerPfos);

    //////////////////////////////////////////////////////////////////
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "LATER SCORES" << std::endl;
    for (const auto& [pPfo, hierarchyPfo] : trackPfos)
        if (hierarchyPfo.GetPredictedParentPfo())
            std::cout << m_isobelID[pPfo] << ": " << m_isobelID[hierarchyPfo.GetPredictedParentPfo()] << ", " << hierarchyPfo.GetLaterTierScore() << std::endl;

    for (const auto& [pPfo, hierarchyPfo] : showerPfos)
        if (hierarchyPfo.GetPredictedParentPfo())
            std::cout << m_isobelID[pPfo] << ": " << m_isobelID[hierarchyPfo.GetPredictedParentPfo()] << ", " << hierarchyPfo.GetLaterTierScore() << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    //////////////////////////////////////////////////////////////////

    // Build the later tier
    this->UpdateHierarchy(pNeutrinoPfo, false, false, m_laterTierThresholdTrackPass1, m_laterTierThresholdShowerPass1, 
        true, trackPfos, showerPfos, hierarchy);

    // Try to recover primaries using laterTierScore
    this->UpdateHierarchy(pNeutrinoPfo, true, false, m_primaryThresholdTrackPass2, m_primaryThresholdShowerPass2, 
        false, trackPfos, showerPfos, hierarchy);

    // Try to recover any children using laterTierScore
    this->UpdateHierarchy(pNeutrinoPfo, false, false, m_laterTierThresholdTrackPass2, m_laterTierThresholdShowerPass2, 
        true, trackPfos, showerPfos, hierarchy);

    //////////////////////////////////////////////////////////////////
    this->PrintHierarchy(hierarchy);
    //////////////////////////////////////////////////////////////////

    this->BuildPandoraHierarchy(pNeutrinoPfo, trackPfos, showerPfos);

    //////////////////////////////////////////////////////////////////
    this->PrintPandoraHierarchy(pNeutrinoPfo);
    //////////////////////////////////////////////////////////////////

    this->CheckForOrphans();

    ///////////////////////////
    // Reset member variables
    m_isobelID.clear();
    ///////////////////////////

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
    float primaryScore(m_bogusFloat);

    m_primaryHierarchyTool->Run(this, pNeutrinoPfo, trackPfos, hierarchyPfo, primaryScore);

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

                if ((thisScore > highestScore) || ((std::fabs(thisScore - highestScore) < std::numeric_limits<float>::epsilon()) &&
                    (this->GetNSpacepoints(pParentPfo) > highestNHits)))
                {
                    highestScore = thisScore;
                    highestNHits = this->GetNSpacepoints(pParentPfo);
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
    float laterTierScore(m_bogusFloat);

    m_laterTierHierarchyTool->Run(this, pNeutrinoPfo, parentPfo, childPfo, laterTierScore);

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

StatusCode MLPNeutrinoHierarchyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
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

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::DetermineIsobelID()
{
    PfoVector pfoVector;

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject * pPfo : *pPfoList)
            pfoVector.push_back(pPfo);
    }

    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject * pPfo : pfoVector)
        m_isobelID[pPfo] = m_isobelID.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::PrintHierarchy(const Hierarchy &hierarchy) const
{
    int count = 2;

    for (const std::vector<const ParticleFlowObject*> &hierarchyTier : hierarchy)
    {
        std::cout << "-------" << std::endl;
        std::cout << "Hierarchy Tier: " << count << std::endl;
        std::cout << "-------" << std::endl;

        for (const ParticleFlowObject * pPfo : hierarchyTier)
            std::cout << "PFP: " << m_isobelID.at(pPfo) << std::endl;

        ++count;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::PrintPandoraHierarchy(const ParticleFlowObject *const pNeutrinoPfo) const
{
    PfoVector parentPfos;
    parentPfos.push_back(pNeutrinoPfo);

    int count = 2;

    while (!parentPfos.empty())
    {
        std::cout << "-------" << std::endl;
        std::cout << "Pandora Tier: " << count << std::endl;
        std::cout << "-------" << std::endl;

        PfoVector tierPfos;

        for (const ParticleFlowObject *const pPfo : parentPfos)
        {
            const PfoList childList(pPfo->GetDaughterPfoList());

            PfoVector childVector;
            for (const ParticleFlowObject *const pChildPfo : childList)
                childVector.emplace_back(pChildPfo);

            std::sort(childVector.begin(), childVector.end(), LArPfoHelper::SortByNHits);

            for (const ParticleFlowObject *const pChildPfo : childVector)
            {
                std::cout << "PFP: " << m_isobelID.at(pChildPfo) << std::endl;
                tierPfos.push_back(pChildPfo);
            }
        }

        parentPfos = tierPfos;
        ++count;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::CheckForOrphans() const
{
    PfoVector pfoVector;

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject * pPfo : *pPfoList)
        {
            if (pPfo->GetParentPfoList().empty())
            {
                std::cout << "ISOBEL THERE IS A SAD LITTLE ORPHAN" << std::endl;
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------


} // namespace lar_dl_content
