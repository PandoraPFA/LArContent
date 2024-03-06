/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronInitialRegionRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the electron initial region refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArShowerRefinement/ElectronInitialRegionRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h"

using namespace pandora;

namespace lar_content
{

ElectronInitialRegionRefinementAlgorithm::ElectronInitialRegionRefinementAlgorithm() :
    m_minShowerHits3D(50),
    m_showerSlidingFitWindow(1000),
    m_maxCoincidenceTransverseSeparation(5.f),
    m_minSpinePurity(0.7f),
    m_trainingMode(false),
    m_trainingFileName("ConnectionPathwayTrain.txt"),
    m_unambiguousThreshold(0.5f),
    m_maxConnectionDistance(1.f),
    m_minNConnectedHits(2),
    m_minElectronCompleteness(0.33f),
    m_minElectronPurity(0.5f),
    m_maxSeparationFromHit(3.f),
    m_maxProjectionSeparation(5.f),
    m_maxXSeparation(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronInitialRegionRefinementAlgorithm::Run()
{
    PfoVector showerPfoVector;
    this->FillShowerPfoVector(showerPfoVector);

    if (showerPfoVector.empty())
        return STATUS_CODE_SUCCESS;

    for (const ParticleFlowObject *const pShowerPfo : showerPfoVector)
    {
        // Only consider significant showers
        CaloHitList caloHits3D;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

        if (caloHits3D.size() < m_minShowerHits3D)
            continue;

        this->RefineShower(pShowerPfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::FillShowerPfoVector(PfoVector &showerPfoVector) const
{
    const PfoList *pPfoList(nullptr);

    if (PandoraContentApi::GetList(*this, m_showerPfoListName, pPfoList) != STATUS_CODE_SUCCESS)
        return;

    if (!pPfoList || pPfoList->empty())
    {
        std::cout << "ElectronInitialRegionRefinementAlgorithm: unable to find shower pfo list " << m_showerPfoListName << std::endl;
        return;
    }

    showerPfoVector.insert(showerPfoVector.begin(), pPfoList->begin(), pPfoList->end());

    std::sort(showerPfoVector.begin(), showerPfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::RefineShower(const ParticleFlowObject *const pShowerPfo) const
{
    CartesianVector nuVertex3D(0.f, 0.f, 0.f);

    if (this->GetNeutrinoVertex(nuVertex3D) != STATUS_CODE_SUCCESS)
        return;

    // Create the 2D connetion pathways
    ProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_U, protoShowerVectorU);
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_V, protoShowerVectorV);
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_W, protoShowerVectorW);

    // 2D->3D connection pathway matching
    ProtoShowerMatchVector protoShowerMatchVector;
    m_pProtoShowerMatchingTool->Run(protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, protoShowerMatchVector);

    if (protoShowerMatchVector.empty())
        return;

    for (ProtoShowerMatch &protoShowerMatch : protoShowerMatchVector)
    {
        // Remove ambiguous hits from hits to add list
        ConnectionPathwayVector viewPathwaysU, viewPathwaysV, viewPathwaysW;

        this->BuildViewPathways(pShowerPfo, protoShowerMatch.GetProtoShowerU().GetSpineHitList(), nuVertex3D, TPC_VIEW_U, viewPathwaysU);
        this->BuildViewPathways(pShowerPfo, protoShowerMatch.GetProtoShowerV().GetSpineHitList(), nuVertex3D, TPC_VIEW_V, viewPathwaysV);
        this->BuildViewPathways(pShowerPfo, protoShowerMatch.GetProtoShowerW().GetSpineHitList(), nuVertex3D, TPC_VIEW_W, viewPathwaysW);

        this->RefineHitsToAdd(nuVertex3D, TPC_VIEW_U, viewPathwaysU, protoShowerMatch.GetProtoShowerToModify(TPC_VIEW_U));
        this->RefineHitsToAdd(nuVertex3D, TPC_VIEW_V, viewPathwaysV, protoShowerMatch.GetProtoShowerToModify(TPC_VIEW_V));
        this->RefineHitsToAdd(nuVertex3D, TPC_VIEW_W, viewPathwaysW, protoShowerMatch.GetProtoShowerToModify(TPC_VIEW_W));

        // Determine the 3D shower vertex
        CartesianPointVector showerStarts3D;
        if (!LArConnectionPathwayHelper::FindShowerStarts3D(this, pShowerPfo, protoShowerMatch, nuVertex3D, m_maxSeparationFromHit,
                m_maxProjectionSeparation, m_maxXSeparation, showerStarts3D))
        {
            continue;
        }

        // Fill BDT information
        StringVector featureOrder;
        const LArMvaHelper::MvaFeatureMap featureMap(LArMvaHelper::CalculateFeatures(
            m_algorithmToolNames, m_featureToolMap, featureOrder, this, pShowerPfo, nuVertex3D, protoShowerMatch, showerStarts3D));

        this->SetMetadata(pShowerPfo, featureMap);

        if (m_trainingMode)
        {
            // To truth match electrons
            HitOwnershipMap electronHitMap;
            this->FillElectronHitMap(electronHitMap);

            if (this->IsElectron(pShowerPfo, electronHitMap))
            {
                LArMvaHelper::ProduceTrainingExample(m_trainingFileName, true, featureOrder, featureMap);
                break;
            }

            LArMvaHelper::ProduceTrainingExample(m_trainingFileName, false, featureOrder, featureMap);

            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronInitialRegionRefinementAlgorithm::GetNeutrinoVertex(CartesianVector &nuVertex3D) const
{
    const VertexList *pNuVertexList(nullptr);
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_neutrinoVertexListName, pNuVertexList));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ElectronInitialRegionRefinementAlgorithm: unable to find vertex list " << m_neutrinoVertexListName
                      << " if it does exist, it may have more than one nu vertex" << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    nuVertex3D = pNuVertexList->front()->GetPosition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::BuildViewProtoShowers(const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertex3D, HitType hitType, ProtoShowerVector &protoShowerVector) const
{
    const CaloHitList *pViewHitList(nullptr);

    if (this->GetHitListOfType(hitType, pViewHitList) != STATUS_CODE_SUCCESS)
        return;

    CartesianVector showerVertexPosition(0.f, 0.f, 0.f);
    try
    {
        showerVertexPosition = this->GetShowerVertex(pShowerPfo, hitType, nuVertex3D);
    }
    catch (...)
    {
        return;
    }

    // Determine directions of pathways out of neutrino vertex
    CartesianPointVector peakDirectionVector;
    if (m_pShowerPeakDirectionFinderTool->Run(pShowerPfo, nuVertex3D, pViewHitList, hitType, peakDirectionVector) != STATUS_CODE_SUCCESS)
        return;

    // Investigate each direction
    CaloHitList unavailableHitList;
    for (CartesianVector &peakDirection : peakDirectionVector)
    {
        // Collect the hits associated with the pathway (the shower spine)
        CaloHitList showerSpineHitList;
        if (m_pShowerSpineFinderTool->Run(nuVertex3D, pViewHitList, hitType, peakDirection, unavailableHitList, showerSpineHitList) != STATUS_CODE_SUCCESS)
            continue;

        this->RefineShowerVertex(pShowerPfo, hitType, nuVertex3D, peakDirection, showerVertexPosition);

        // Quality check: If the spine passes the shower vertex, does it live inside the shower?
        if (!this->IsSpineCoincident(pShowerPfo, nuVertex3D, hitType, showerVertexPosition, showerSpineHitList))
            continue;

        // Find the 2D shower start position
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        CartesianVector showerStartDirection(0.f, 0.f, 0.f);

        if (m_pShowerStartFinderTool->Run(pShowerPfo, peakDirection, hitType, showerSpineHitList, showerStartPosition, showerStartDirection) != STATUS_CODE_SUCCESS)
            continue;

        // Create the ProtoShower object
        const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

        ProtoShower protoShower(ShowerCore(showerStartPosition, showerStartDirection), ConnectionPathway(nuVertex2D, peakDirection),
            showerSpineHitList, CaloHitList(), CartesianPointVector(), CaloHitList());

        // Now determine the hits to be added to the shower from the found spine
        CaloHitList viewShowerHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

        const float showerVertexL(std::max(
            peakDirection.GetDotProduct(showerStartPosition - nuVertex2D), peakDirection.GetDotProduct(showerVertexPosition - nuVertex2D)));

        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            if (std::find(viewShowerHitList.begin(), viewShowerHitList.end(), pCaloHit) != viewShowerHitList.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const float showerL(peakDirection.GetDotProduct(hitPosition - nuVertex2D));

            if ((showerL > 0.f) && (showerL < showerVertexL))
                protoShower.AddHitToAdd(pCaloHit);
        }

        protoShowerVector.push_back(protoShower);
        unavailableHitList.insert(unavailableHitList.begin(), showerSpineHitList.begin(), showerSpineHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronInitialRegionRefinementAlgorithm::GetHitListOfType(const HitType hitType, const CaloHitList *&pCaloHitList) const
{
    const std::string typeHitListName(hitType == TPC_VIEW_U   ? m_caloHitListNameU
                                      : hitType == TPC_VIEW_V ? m_caloHitListNameV
                                                              : m_caloHitListNameW);

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, typeHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementBaseTool: unable to find calo hit list " << typeHitListName << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector ElectronInitialRegionRefinementAlgorithm::GetShowerVertex(
    const ParticleFlowObject *const pShowerPfo, const HitType hitType, const CartesianVector &nuVertex3D) const
{
    ClusterList viewCusterList;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, viewCusterList);

    if (viewCusterList.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const TwoDSlidingShowerFitResult twoDShowerSlidingFit(
        viewCusterList.front(), m_showerSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
    const CartesianVector &minLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMinLayerPosition());
    const CartesianVector &maxLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMaxLayerPosition());
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    if ((nuVertex2D.GetZ() > minLayerPosition.GetZ()) && (nuVertex2D.GetZ() < maxLayerPosition.GetZ()))
        return nuVertex2D;

    const float minSeparation((nuVertex2D - minLayerPosition).GetMagnitudeSquared());
    const float maxSeparation((nuVertex2D - maxLayerPosition).GetMagnitudeSquared());
    CartesianVector showerVertexPosition(minSeparation < maxSeparation ? minLayerPosition : maxLayerPosition);

    return (minSeparation < maxSeparation ? minLayerPosition : maxLayerPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::RefineShowerVertex(const ParticleFlowObject *const pShowerPfo, const HitType hitType,
    const CartesianVector &nuVertex3D, const CartesianVector &peakDirection, CartesianVector &showerVertexPosition) const
{
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    if (!this->IsShowerConnected(showerVertexPosition, nuVertex2D, peakDirection))
    {
        CaloHitList viewShowerHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

        float minL(std::numeric_limits<float>::max());

        for (const CaloHit *const pCaloHit : viewShowerHitList)
        {
            const CartesianVector displacement(pCaloHit->GetPositionVector() - nuVertex2D);
            const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

            if ((longitudinalSeparation < (-1.f)) || (longitudinalSeparation > minL))
                continue;

            const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

            if (transverseSeparation < m_maxCoincidenceTransverseSeparation)
            {
                showerVertexPosition = pCaloHit->GetPositionVector();
                minL = longitudinalSeparation;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronInitialRegionRefinementAlgorithm::IsShowerConnected(
    const CartesianVector &showerVertexPosition, const CartesianVector &nuVertex2D, const CartesianVector &peakDirection) const
{
    CartesianVector displacement(showerVertexPosition - nuVertex2D);

    const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

    if (longitudinalSeparation < (-1.f))
        return false;

    const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

    if (transverseSeparation > m_maxCoincidenceTransverseSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronInitialRegionRefinementAlgorithm::IsSpineCoincident(const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertex3D, const HitType hitType, const CartesianVector &showerVertex, const CaloHitList &showerSpineHitList) const
{
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

    CaloHitList postShowerVertexSpineHits;
    const float showerDistanceFromNuSquared((showerVertex - nuVertex2D).GetMagnitudeSquared());

    for (const CaloHit *const pSpineHit : showerSpineHitList)
    {
        const CartesianVector &hitPosition(pSpineHit->GetPositionVector());
        const float separationSquared((hitPosition - nuVertex2D).GetMagnitudeSquared());

        if (separationSquared > showerDistanceFromNuSquared)
            postShowerVertexSpineHits.push_back(pSpineHit);
    }

    if (postShowerVertexSpineHits.size() == 0)
        return true;

    // Check whether shower spine is pure
    const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(postShowerVertexSpineHits, viewShowerHitList));
    const float spinePurity(static_cast<float>(sharedHitList.size()) / static_cast<float>(postShowerVertexSpineHits.size()));

    if (spinePurity < m_minSpinePurity)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::BuildViewPathways(const ParticleFlowObject *const pShowerPfo,
    const CaloHitList &protectedHits, const CartesianVector &nuVertex3D, HitType hitType, ConnectionPathwayVector &viewPathways) const
{
    const CaloHitList *pViewHitList(nullptr);

    if (this->GetHitListOfType(hitType, pViewHitList) != STATUS_CODE_SUCCESS)
        return;

    // Get the peak direction vector
    CartesianPointVector eventPeakDirectionVector;
    m_pEventPeakDirectionFinderTool->Run(pShowerPfo, nuVertex3D, pViewHitList, hitType, eventPeakDirectionVector);

    CaloHitList unavailableHitList(protectedHits);
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, unavailableHitList);

    for (CartesianVector &eventPeakDirection : eventPeakDirectionVector)
    {
        CaloHitList pathwayHitList;
        if (m_pEventPathwayFinderTool->Run(nuVertex3D, pViewHitList, hitType, eventPeakDirection, unavailableHitList, pathwayHitList) != STATUS_CODE_SUCCESS)
            continue;

        const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));
        const ConnectionPathway connectionPathway(nuVertex2D, eventPeakDirection);

        viewPathways.push_back(connectionPathway);
        unavailableHitList.insert(unavailableHitList.begin(), pathwayHitList.begin(), pathwayHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::RefineHitsToAdd(
    const CartesianVector &nuVertex3D, const HitType hitType, const ConnectionPathwayVector &viewPathways, ProtoShower &protoShower) const
{
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));
    const CartesianVector &peakDirection(protoShower.GetConnectionPathway().GetStartDirection());

    CaloHitList refinedHitList;

    for (const CaloHit *const pHitToAdd : protoShower.GetHitsToAddList())
    {
        bool found(false);

        const CartesianVector &hitPosition(pHitToAdd->GetPositionVector());
        const CartesianVector displacement(hitPosition - nuVertex2D);
        const float thisT((peakDirection.GetCrossProduct(displacement)).GetMagnitudeSquared());

        for (const ConnectionPathway &connectionPathway : viewPathways)
        {
            const CartesianVector &eventPeakDirection(connectionPathway.GetStartDirection());

            if (peakDirection.GetOpeningAngle(eventPeakDirection) > (M_PI / 2))
                continue;

            const float otherT((eventPeakDirection.GetCrossProduct(displacement)).GetMagnitudeSquared());

            if ((otherT < thisT) || (otherT < m_unambiguousThreshold))
            {
                found = true;

                if (std::find(protoShower.GetAmbiguousDirectionVector().begin(), protoShower.GetAmbiguousDirectionVector().end(),
                        eventPeakDirection) == protoShower.GetAmbiguousDirectionVector().end())
                {
                    protoShower.AddAmbiguousDirection(connectionPathway.GetStartDirection());
                }
            }
        }

        found ? protoShower.AddAmbiguousHit(pHitToAdd) : refinedHitList.push_back(pHitToAdd);
    }

    CaloHitList continuousHitList;
    this->FindContinuousPath(refinedHitList, nuVertex2D, continuousHitList);

    protoShower.SetHitsToAddList(continuousHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::FindContinuousPath(
    const CaloHitList &refinedHitList, const CartesianVector &nuVertex2D, CaloHitList &continuousHitList) const
{
    continuousHitList.clear();

    CaloHitVector refinedHitVector(refinedHitList.begin(), refinedHitList.end());
    std::sort(refinedHitVector.begin(), refinedHitVector.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(nuVertex2D));

    unsigned int startIndex(refinedHitVector.size());

    for (unsigned int i = 0; i < refinedHitVector.size(); ++i)
    {
        CaloHitList connectedHitList;
        connectedHitList.push_back(refinedHitVector[i]);

        bool found(true);

        while (found)
        {
            found = false;

            for (unsigned int j = (i + 1); j < refinedHitVector.size(); ++j)
            {
                const CaloHit *const pCaloHit(refinedHitVector[j]);

                if (std::find(connectedHitList.begin(), connectedHitList.end(), pCaloHit) != connectedHitList.end())
                    continue;

                if (LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), connectedHitList) < m_maxConnectionDistance)
                {
                    found = true;
                    connectedHitList.push_back(pCaloHit);
                    break;
                }
            }
        }

        if ((connectedHitList.size() >= m_minNConnectedHits) || (connectedHitList.size() == refinedHitVector.size()))
        {
            startIndex = i;
            break;
        }
    }

    for (unsigned int i = startIndex; i < refinedHitVector.size(); ++i)
        continuousHitList.push_back(refinedHitVector[i]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::SetMetadata(const ParticleFlowObject *const pShowerPfo, const LArMvaHelper::MvaFeatureMap &featureMap) const
{
    object_creation::ParticleFlowObject::Metadata metadata;

    metadata.m_propertiesToAdd["FoundConnectionPathway"] = true;

    if (featureMap.find("LArInitialRegionFeatureTool_initialGapSize") != featureMap.end())
        metadata.m_propertiesToAdd["MaxInitialGapSize"] = featureMap.at("LArInitialRegionFeatureTool_initialGapSize").Get();

    if (featureMap.find("LArInitialRegionFeatureTool_largestGapSize") != featureMap.end())
        metadata.m_propertiesToAdd["MinLargestProjectedGapSize"] = featureMap.at("LArInitialRegionFeatureTool_largestGapSize").Get();

    if (featureMap.find("LArConnectionRegionFeatureTool_pathwayLength") != featureMap.end())
        metadata.m_propertiesToAdd["PathwayLengthMin"] = featureMap.at("LArConnectionRegionFeatureTool_pathwayLength").Get();

    if (featureMap.find("LArConnectionRegionFeatureTool_pathwayScatteringAngle2D") != featureMap.end())
        metadata.m_propertiesToAdd["MaxShowerStartPathwayScatteringAngle2D"] =
            featureMap.at("LArConnectionRegionFeatureTool_pathwayScatteringAngle2D").Get();

    if (featureMap.find("LArShowerRegionFeatureTool_nShowerHits") != featureMap.end())
        metadata.m_propertiesToAdd["MaxNPostShowerStartHits"] = featureMap.at("LArShowerRegionFeatureTool_nShowerHits").Get();

    if (featureMap.find("LArShowerRegionFeatureTool_foundHitRatio") != featureMap.end())
        metadata.m_propertiesToAdd["MaxFoundHitRatio"] = featureMap.at("LArShowerRegionFeatureTool_foundHitRatio").Get();

    if (featureMap.find("LArShowerRegionFeatureTool_scatterAngle") != featureMap.end())
        metadata.m_propertiesToAdd["MaxPostShowerStartScatterAngle"] = featureMap.at("LArShowerRegionFeatureTool_scatterAngle").Get();

    if (featureMap.find("LArShowerRegionFeatureTool_openingAngle") != featureMap.end())
        metadata.m_propertiesToAdd["MaxPostShowerStartOpeningAngle"] = featureMap.at("LArShowerRegionFeatureTool_openingAngle").Get();

    if (featureMap.find("LArShowerRegionFeatureTool_nuVertexEnergyAsymmetry") != featureMap.end())
        metadata.m_propertiesToAdd["MaxPostShowerStartNuVertexEnergyAsymmetry"] =
            featureMap.at("LArShowerRegionFeatureTool_nuVertexEnergyAsymmetry").Get();

    if (featureMap.find("LArShowerRegionFeatureTool_nuVertexEnergyWeightedMeanRadialDistance") != featureMap.end())
        metadata.m_propertiesToAdd["MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance"] =
            featureMap.at("LArShowerRegionFeatureTool_nuVertexEnergyWeightedMeanRadialDistance").Get();

    if (featureMap.find("LArShowerRegionFeatureTool_showerStartEnergyAsymmetry") != featureMap.end())
        metadata.m_propertiesToAdd["MaxPostShowerStartShowerStartEnergyAsymmetry"] =
            featureMap.at("LArShowerRegionFeatureTool_showerStartEnergyAsymmetry").Get();

    if (featureMap.find("LArShowerRegionFeatureTool_showerStartMoliereRadius") != featureMap.end())
        metadata.m_propertiesToAdd["MinPostShowerStartShowerStartMoliereRadius"] =
            featureMap.at("LArShowerRegionFeatureTool_showerStartMoliereRadius").Get();

    if (featureMap.find("LArAmbiguousRegionFeatureTool_nAmbiguousViews") != featureMap.end())
        metadata.m_propertiesToAdd["NViewsWithAmbiguousHits"] = featureMap.at("LArAmbiguousRegionFeatureTool_nAmbiguousViews").Get();

    if (featureMap.find("LArAmbiguousRegionFeatureTool_maxUnaccountedEnergy") != featureMap.end())
        metadata.m_propertiesToAdd["AmbiguousHitMaxUnaccountedEnergy"] = featureMap.at("LArAmbiguousRegionFeatureTool_maxUnaccountedEnergy").Get();

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pShowerPfo, metadata));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::FillElectronHitMap(HitOwnershipMap &electronHitMap) const
{
    electronHitMap.clear();

    for (HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const CaloHitList *pViewHitList(nullptr);

        if (this->GetHitListOfType(hitType, pViewHitList) != STATUS_CODE_SUCCESS)
            continue;

        for (const CaloHit *const pCaloHit : *pViewHitList)
        {
            MCParticleVector contributingMCParticleVector;
            const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());

            for (const auto &mapEntry : weightMap)
                contributingMCParticleVector.push_back(mapEntry.first);

            std::sort(contributingMCParticleVector.begin(), contributingMCParticleVector.end(), PointerLessThan<MCParticle>());

            float highestWeight(0.f);
            const MCParticle *highestElectronContributor(nullptr);

            for (const MCParticle *const pMCParticle : contributingMCParticleVector)
            {
                const bool isLeadingElectron((std::abs(pMCParticle->GetParticleId()) == 11) &&
                                             (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle));

                if (isLeadingElectron)
                {
                    const float weight(weightMap.at(pMCParticle));

                    if (weight > highestWeight)
                    {
                        highestWeight = weight;
                        highestElectronContributor = pMCParticle;
                    }
                }
            }

            if (highestElectronContributor)
                electronHitMap[highestElectronContributor].push_back(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronInitialRegionRefinementAlgorithm::IsElectron(const ParticleFlowObject *const pPfo, const HitOwnershipMap &electronHitMap) const
{
    MCParticleVector mcElectronVector;

    for (auto &entry : electronHitMap)
        mcElectronVector.push_back(entry.first);

    CaloHitList pfoHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, pfoHitList);

    for (const MCParticle *const pMCElectron : mcElectronVector)
    {
        const CaloHitList &mcElectronHitList(electronHitMap.at(pMCElectron));
        const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(pfoHitList, mcElectronHitList));

        const float completeness(static_cast<float>(sharedHitList.size()) / static_cast<float>(mcElectronHitList.size()));
        const float purity(static_cast<float>(sharedHitList.size()) / static_cast<float>(pfoHitList.size()));

        if (completeness < m_minElectronCompleteness)
            continue;

        if (purity < m_minElectronPurity)
            continue;

        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronInitialRegionRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameU", m_caloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameV", m_caloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameW", m_caloHitListNameW));

    AlgorithmTool *pAlgorithmTool1(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerPeakDirectionFinder", pAlgorithmTool1));
    m_pShowerPeakDirectionFinderTool = dynamic_cast<PeakDirectionFinderTool *>(pAlgorithmTool1);

    if (!m_pShowerPeakDirectionFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool2(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "EventPeakDirectionFinder", pAlgorithmTool2));
    m_pEventPeakDirectionFinderTool = dynamic_cast<PeakDirectionFinderTool *>(pAlgorithmTool2);

    if (!m_pEventPeakDirectionFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool3(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerSpineFinder", pAlgorithmTool3));
    m_pShowerSpineFinderTool = dynamic_cast<ShowerSpineFinderTool *>(pAlgorithmTool3);

    if (!m_pShowerSpineFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool4(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "EventPathwayFinder", pAlgorithmTool4));
    m_pEventPathwayFinderTool = dynamic_cast<ShowerSpineFinderTool *>(pAlgorithmTool4);

    if (!m_pEventPathwayFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool5(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerStartFinder", pAlgorithmTool5));
    m_pShowerStartFinderTool = dynamic_cast<ShowerStartFinderTool *>(pAlgorithmTool5);

    if (!m_pShowerSpineFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool6(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ProtoShowerMatching", pAlgorithmTool6));
    m_pProtoShowerMatchingTool = dynamic_cast<ProtoShowerMatchingTool *>(pAlgorithmTool6);

    if (!m_pProtoShowerMatchingTool)
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinShowerHits3D", m_minShowerHits3D));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowerSlidingFitWindow", m_showerSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxCoincidenceTransverseSeparation", m_maxCoincidenceTransverseSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSpinePurity", m_minSpinePurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UnambiguousThreshold", m_unambiguousThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxConnectionDistance", m_maxConnectionDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinNConnectedHits", m_minNConnectedHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinElectronCompleteness", m_minElectronCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinElectronPurity", m_minElectronPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSeparationFromHit", m_maxSeparationFromHit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxProjectionSeparation", m_maxProjectionSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    LArMvaHelper::AlgorithmToolMap algorithmToolMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        LArMvaHelper::ProcessAlgorithmToolListToMap(*this, xmlHandle, "FeatureTools", m_algorithmToolNames, algorithmToolMap));

    for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMap)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMap));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
