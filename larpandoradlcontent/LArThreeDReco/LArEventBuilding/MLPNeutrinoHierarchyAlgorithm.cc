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
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPLaterTierHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPNeutrinoHierarchyAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MLPNeutrinoHierarchyAlgorithm::MLPNeutrinoHierarchyAlgorithm() :
    m_primaryThresholdTrackPass1(0.8f),
    m_primaryThresholdShowerPass1(0.45f),
    m_laterTierThresholdTrackPass1(0.8f),
    m_laterTierThresholdShowerPass1(0.8f),
    m_pNeutrinoPfo(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPNeutrinoHierarchyAlgorithm::Run()
{
    srand(0.1);
    std::cout << "RUNNING THE NEW HIERARCHY ALG!" << std::endl;

    if (!this->GetNeutrinoPfo())
        return STATUS_CODE_SUCCESS;

    // Give PFPs IDs to keep track of them
    this->DetermineIsobelID();

    // Fill the track/shower vectors
    this->FillTrackShowerVectors();

    // Calculate primary scores
    this->SetPrimaryScores();

    //////////////////////////////////////////////////////////////////
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "PRIMARY SCORES" << std::endl;
    for (const auto& [pPfo, hierarchyPfo] : m_trackPfos)
        std::cout << m_isobelID[pPfo] << ": " << hierarchyPfo.GetPrimaryScore() << std::endl;

    for (const auto& [pPfo, hierarchyPfo] : m_showerPfos)
        std::cout << m_isobelID[pPfo] << ": " << hierarchyPfo.GetPrimaryScore() << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    //////////////////////////////////////////////////////////////////

    // Set what primary things we can
    this->BuildPrimaryTierPass1();

    if (m_hierarchy.size() == 0)
    {
        std::cout << "NO PRIMARIES FOUND!" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    //////////////////////////////////////////////////////////////////    
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "PRIMARY TIER:" << std::endl;
    for (const ParticleFlowObject* pPfo : m_hierarchy.at(0))
        std::cout << m_isobelID[pPfo] << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    //////////////////////////////////////////////////////////////////    

    // Set later tier scores
    this->SetLaterTierScores();

    //////////////////////////////////////////////////////////////////
    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "LATER SCORES" << std::endl;
    for (const auto& [pPfo, hierarchyPfo] : m_trackPfos)
        if (hierarchyPfo.GetPredictedParentPfo())
            std::cout << m_isobelID[pPfo] << ": " << m_isobelID[hierarchyPfo.GetPredictedParentPfo()] << ", " << hierarchyPfo.GetLaterTierScore() << std::endl;

    for (const auto& [pPfo, hierarchyPfo] : m_showerPfos)
        if (hierarchyPfo.GetPredictedParentPfo())
            std::cout << m_isobelID[pPfo] << ": " << m_isobelID[hierarchyPfo.GetPredictedParentPfo()] << ", " << hierarchyPfo.GetLaterTierScore() << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
    //////////////////////////////////////////////////////////////////

    // Build the later tier
    this->BuildLaterTierPass1();

    //////////////////////////////////////////////////////////////////
    this->PrintHierarchy();
    //////////////////////////////////////////////////////////////////

    this->BuildPandoraHierarchy();

    //////////////////////////////////////////////////////////////////
    this->PrintPandoraHierarchy();
    //////////////////////////////////////////////////////////////////

    this->CheckForOrphans();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPNeutrinoHierarchyAlgorithm::GetNeutrinoPfo()
{
    const PfoList *pPfoList = nullptr;

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, 
        PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
        return false;

    // ATTN Enforces that only one pfo, of neutrino-type, be in the specified input list
    m_pNeutrinoPfo = (1 == pPfoList->size()) ? *(pPfoList->begin()) : nullptr;

    if (!m_pNeutrinoPfo || !LArPfoHelper::IsNeutrino(m_pNeutrinoPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::FillTrackShowerVectors()
{
    const int MIN_3D_CLUSTER_SIZE = 5;

    // Sliding fit shenanigans
    const int HALF_WINDOW_LAYERS(20);
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
            if (this->GetNSpacepoints(pPfo) < MIN_3D_CLUSTER_SIZE)
                continue;

            // Attempt sliding linear fit
            try
            {
                CartesianPointVector pointVector;
                LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, pointVector);

                const ThreeDSlidingFitResult slidingFitResult(&pointVector, HALF_WINDOW_LAYERS, slidingFitPitch);

                // We need directions...
                CartesianVector upstreamVertex(-999.f, -999.f, -999.f), upstreamDirection(-999.f, -999.f, -999.f), 
                    downstreamVertex(-999.f, -999.f, -999.f), downstreamDirection(-999.f, -999.f, -999.f);

                if (!this->GetExtremalVerticesAndDirections(pPfo, slidingFitResult, upstreamVertex, upstreamDirection, downstreamVertex, downstreamDirection))
                    continue;

                // Create track/shower objects
                if (pPfo->GetParticleId() == 13)
                {
                    m_trackPfos.insert(HierarchyPfoMapEntry({pPfo, HierarchyPfo(true, pPfo, slidingFitResult, upstreamVertex, upstreamDirection, downstreamVertex, downstreamDirection)}));
                }
                else if (pPfo->GetParticleId() == 11) 
                {
                    m_showerPfos.insert(HierarchyPfoMapEntry({pPfo, HierarchyPfo(false, pPfo, slidingFitResult, upstreamVertex, upstreamDirection, downstreamVertex, downstreamDirection)}));
                }
            }
            catch (...) { continue; }
        }
    }

    std::cout << "We have " << m_trackPfos.size() << " track(s) and " << m_showerPfos.size() << " shower(s)" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetNSpacepoints(const ParticleFlowObject *const pPfo)
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    int total3DHits(0);

    for (const Cluster *const pCluster3D : clusterList3D)
        total3DHits += pCluster3D->GetNCaloHits();

    return total3DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPNeutrinoHierarchyAlgorithm::GetExtremalVerticesAndDirections(const ParticleFlowObject *const pPfo, const ThreeDSlidingFitResult &slidingFitResult,
    CartesianVector &upstreamVertex, CartesianVector &upstreamDirection, CartesianVector &downstreamVertex, CartesianVector &downstreamDirection)
{
    if (m_pNeutrinoPfo->GetVertexList().empty())
        return false;
 
    // First, get the neutrino vertex
    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(m_pNeutrinoPfo));
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
        const bool upstreamDirectionSet(this->GetShowerDirection(pPfo, upstreamVertex, 25.f, upstreamDirection));
        const bool downstreamDirectionSet(this->GetShowerDirection(pPfo, downstreamVertex, 25.f, downstreamDirection));

        if (!upstreamDirectionSet && !downstreamDirectionSet)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPNeutrinoHierarchyAlgorithm::GetShowerDirection(const ParticleFlowObject *const pPfo, const CartesianVector &vertex, const float searchRegion, 
    CartesianVector &direction)
{
    CartesianPointVector pointVector;
    LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, pointVector);

    if (pointVector.empty())
        return false;

    // set direction by angular distribution
    const int nBins(180);
    const float angleMin(0.f), angleMax(2.f * M_PI);
    const float binWidth((angleMax - angleMin) / static_cast<float>(nBins));

    std::vector<std::vector<int>> spatialDist(nBins, std::vector<int>(nBins, 0));
    std::vector<std::vector<float>> tieBreakerDist(nBins, std::vector<float>(nBins, 0.f));

    int highestSP(0), bestTheta0YZBin(-1), bestTheta0XZBin(-1);
    float tieBreaker(std::numeric_limits<float>::max());

    for (const CartesianVector &position3D : pointVector)
    {
        const CartesianVector displacement(position3D - vertex);
        const float mag = displacement.GetMagnitude();

        if (mag > searchRegion)
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

void MLPNeutrinoHierarchyAlgorithm::SetPrimaryScores()
{
    // For tracks
    for (auto& [pPfo, hierarchyPfo] : m_trackPfos)
        m_primaryHierarchyTool->Run(this, m_pNeutrinoPfo, m_trackPfos, hierarchyPfo);

    // For showers
    for (auto& [pPfo, hierarchyPfo] : m_showerPfos)
        m_primaryHierarchyTool->Run(this, m_pNeutrinoPfo, m_trackPfos, hierarchyPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::BuildPrimaryTierPass1()
{
    std::vector<const ParticleFlowObject*> primaryTier;

    for (bool isTrack : {true, false})
    {
        std::map<const ParticleFlowObject*, HierarchyPfo> &hierarchyPfoMap(isTrack ? m_trackPfos : m_showerPfos);

        for (auto& [pPfo, hierarchyPfo] : hierarchyPfoMap)
        {
            const float thisPrimaryScore(hierarchyPfo.GetPrimaryScore());

            if ((thisPrimaryScore > 0.f) && ((isTrack && (thisPrimaryScore > m_primaryThresholdTrackPass1)) || 
                (!isTrack && (thisPrimaryScore > m_primaryThresholdShowerPass1))))
            {
                // Add to primary tier
                primaryTier.emplace_back(pPfo);

                // Add to hierarchyPfo info
                hierarchyPfo.SetParentPfo(m_pNeutrinoPfo);
                hierarchyPfo.SetIsInHierarchy(true);
            }
        }
    }

    m_hierarchy.emplace_back(primaryTier);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::SetLaterTierScores()
{
    for (bool isTrack : {true, false})
    {
        std::map<const ParticleFlowObject*, HierarchyPfo> &hierarchyPfoMap(isTrack ? m_trackPfos : m_showerPfos);

        for (auto& [pChildPfo, childPfo] : hierarchyPfoMap)
        {
            // Don't bother if we have already found a parent
            if (childPfo.GetIsInHierarchy())
                continue;

            float highestScore(-std::numeric_limits<float>::max());

            for (auto& [pParentPfo, parentPfo] : m_trackPfos)
            {
                if (pChildPfo == pParentPfo)
                    continue;

                int parentOrientation(-1), childOrientation(-1);
                const float thisScore(isTrack ? this->GetLaterTierScoreTrackToTrack(parentPfo, childPfo, parentOrientation, childOrientation) :
                    this->GetLaterTierScoreTrackToShower(parentPfo, childPfo, parentOrientation, childOrientation));

                if (thisScore > highestScore)
                {
                    highestScore = thisScore;
                    childPfo.SetLaterTierScore(thisScore);
                    childPfo.SetParentOrientation(parentOrientation);
                    childPfo.SetChildOrientation(childOrientation);
                    childPfo.SetPredictedParentPfo(pParentPfo);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::BuildLaterTierPass1()
{
    while (m_hierarchy.at(m_hierarchy.size() - 1).size() != 0)
    {
        std::vector<const ParticleFlowObject*> &preceedingTier(m_hierarchy.at(m_hierarchy.size() - 1));
        std::vector<const ParticleFlowObject*> thisTier;

        for (bool isTrack : {true, false})
        {
            std::map<const ParticleFlowObject*, HierarchyPfo> &hierarchyPfoMap(isTrack ? m_trackPfos : m_showerPfos);

            for (auto& [pChildPfo, childPfo] : hierarchyPfoMap)
            {
                // Don't bother if we have already found a parent
                if (childPfo.GetIsInHierarchy())
                    continue;

                // Do we have a predicted parent
                if (!childPfo.GetPredictedParentPfo())
                    continue;

                // Is predicted parent in preceeding tier?
                const ParticleFlowObject *const pPredictedParent = childPfo.GetPredictedParentPfo();
                if (std::find(preceedingTier.begin(), preceedingTier.end(), pPredictedParent) == preceedingTier.end())
                    continue;

                // Does it pass tier cut?
                if ((isTrack && (childPfo.GetLaterTierScore() > m_laterTierThresholdTrackPass1)) ||
                    (!isTrack && (childPfo.GetLaterTierScore() > m_laterTierThresholdShowerPass1)))
                {
                    // Add child to hierarchy tier
                    thisTier.emplace_back(pChildPfo);

                    // Add info to childPfo
                    childPfo.SetIsInHierarchy(true);
                    childPfo.SetParentPfo(pPredictedParent);

                    // Add info to parentPfo
                    if (m_trackPfos.find(pPredictedParent) == m_trackPfos.end())
                        throw StatusCodeException(STATUS_CODE_FAILURE);

                    m_trackPfos.at(pPredictedParent).AddChildPfo(pChildPfo);
                }
            }
        }

        m_hierarchy.emplace_back(thisTier);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetLaterTierScoreTrackToTrack(HierarchyPfo &parentPfo, HierarchyPfo &childPfo, 
    int &parentOrientation, int &childOrientation) const
{
    //m_laterTierHierarchyTool->Run(this, m_pNeutrinoPfo, parentPfo, childPfo);

    return this->GetRandomNumber();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetLaterTierScoreTrackToShower(HierarchyPfo &parentPfo, HierarchyPfo &childPfo,
    int &parentOrientation, int &childOrientation) const
{
    //m_laterTierHierarchyTool->Run(this, m_pNeutrinoPfo, parentPfo, childPfo);

    return this->GetRandomNumber();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::BuildPandoraHierarchy()
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
        HierarchyPfo* hierarchyPfo(nullptr);

        if (m_trackPfos.find(pPfo) != m_trackPfos.end())
        {
            hierarchyPfo = &m_trackPfos.at(pPfo);
        }
        else if (m_showerPfos.find(pPfo) != m_showerPfos.end())
        {
            hierarchyPfo = &m_showerPfos.at(pPfo);
        }
        else
        {
            // If alg never handled the pfo, assign as primary
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, m_pNeutrinoPfo, pPfo));
            continue;
        }

        // If pfo was never assigned to hierarchy, add as primary
        if (!hierarchyPfo->GetIsInHierarchy())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, m_pNeutrinoPfo, pPfo));
            continue;
        }

        // If parent is a neutrino we have to assign its parent
        if (hierarchyPfo->GetParentPfo() == m_pNeutrinoPfo)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, m_pNeutrinoPfo, pPfo));

        // Assign parent-child links of its children
        const PfoVector &childPfos(hierarchyPfo->GetSortedChildPfoVector());

        if ((!childPfos.empty()) && (m_showerPfos.find(pPfo) != m_showerPfos.end()))
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

    AlgorithmTool *pAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "MLPPrimaryHierarchyTool", pAlgorithmTool));
    m_primaryHierarchyTool = dynamic_cast<MLPPrimaryHierarchyTool *>(pAlgorithmTool);

    if (!m_primaryHierarchyTool)
        return STATUS_CODE_INVALID_PARAMETER;

    // pAlgorithmTool = nullptr;
    // PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "MLPLaterTierHierarchyTool", pAlgorithmTool));
    // m_laterTierHierarchyTool = dynamic_cast<MLPLaterTierHierarchyTool *>(pAlgorithmTool);

    // if (!m_laterTierHierarchyTool)
    //     return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryThresholdTrackPass1", m_primaryThresholdTrackPass1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryThresholdShowerPass1", m_primaryThresholdShowerPass1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierThresholdTrackPass1", m_laterTierThresholdTrackPass1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierThresholdShowerPass1", m_laterTierThresholdShowerPass1));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetRandomNumber() const
{
    int randomNumber1 = rand() % 10 + 1;
    int randomNumber2 = rand() % 10 + 1;

    return ((static_cast<float>(randomNumber1) / 10.f) + (static_cast<float>(randomNumber2) / 100.f));
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

void MLPNeutrinoHierarchyAlgorithm::PrintHierarchy()
{
    int count = 2;

    for (const std::vector<const ParticleFlowObject*> &hierarchyTier : m_hierarchy)
    {
        std::cout << "-------" << std::endl;
        std::cout << "m_hierarchy Tier: " << count << std::endl;
        std::cout << "-------" << std::endl;

        for (const ParticleFlowObject * pPfo : hierarchyTier)
            std::cout << "PFP: " << m_isobelID[pPfo] << std::endl;

        ++count;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::PrintPandoraHierarchy()
{
    PfoVector parentPfos;
    parentPfos.push_back(m_pNeutrinoPfo);

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
                std::cout << "PFP: " << m_isobelID[pChildPfo] << std::endl;
                tierPfos.push_back(pChildPfo);
            }
        }

        parentPfos = tierPfos;
        ++count;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::CheckForOrphans()
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
