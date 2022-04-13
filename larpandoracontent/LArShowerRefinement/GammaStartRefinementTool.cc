/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.cc
 *
 *  @brief  Implementation of the gamma start refinement tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "PandoraMonitoringApi.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h"

using namespace pandora;

namespace lar_content
{

GammaStartRefinementTool::GammaStartRefinementTool() : 
    m_counter(0), 
    m_maxXSeparation(3.f),
    m_maxSeparation(2.f),
    m_maxAngularDeviation(1.f),
    m_maxProjectedShowerStartSeparation(2.f)
{
}

GammaStartRefinementTool::~GammaStartRefinementTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition, 
    const CaloHitList *const pCaloHitListU, const CaloHitList *const pCaloHitListV, const CaloHitList *const pCaloHitListW)
{
    if (!LArPfoHelper::IsShower(pShowerPfo))
        return false;

    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

    if (caloHits3D.size() < 100)
        return false;

    if (!pAlgorithm->IsGamma(pShowerPfo, nuVertexPosition))
        return false;

    std::cout << "HAS PASSED TRUTH GAMMA HERE" << std::endl;

    // Keep track of 'taken' hits
    CaloHitList usedHitListU, usedHitListV, usedHitListW;
    ProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, protoShowerVectorU, usedHitListU);
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, protoShowerVectorV, usedHitListV);
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, protoShowerVectorW, usedHitListW);

    if (!this->IsShowerTruncatable(protoShowerVectorU, protoShowerVectorV, protoShowerVectorW))
        return false;

    // Now add in connection pathways that aren't present in some views? (terrible explanation)
    this->BuildHelperProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, protoShowerVectorU, usedHitListU);
    this->BuildHelperProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, protoShowerVectorV, usedHitListV);
    this->BuildHelperProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, protoShowerVectorW, usedHitListW);

    this->AssignShowerHits(pAlgorithm, pShowerPfo, TPC_VIEW_U, protoShowerVectorU);
    this->AssignShowerHits(pAlgorithm, pShowerPfo, TPC_VIEW_V, protoShowerVectorV);
    this->AssignShowerHits(pAlgorithm, pShowerPfo, TPC_VIEW_W, protoShowerVectorW);

    MatchedConnectionPathwayMap threeViewConnectionPathways, twoViewConnectionPathways, oneViewConnectionPathways;

    this->MatchConnectionPathways(pAlgorithm, nuVertexPosition, pShowerPfo, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, 
        threeViewConnectionPathways, twoViewConnectionPathways, oneViewConnectionPathways);

    if (threeViewConnectionPathways.empty() && oneViewConnectionPathways.empty())
        return false;

    // this is to say that the hybrid algorithm was active
    if (!threeViewConnectionPathways.empty())
    {
        object_creation::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["ActiveHybridGammaAlg"] = 1.f;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pShowerPfo, metadata));
    }

    // remove pathways
    this->RemoveThreeViewConnectionPathways(pAlgorithm, threeViewConnectionPathways, pShowerPfo, nuVertexPosition);
    this->RemoveOneViewConnectionPathways(pAlgorithm, oneViewConnectionPathways, pShowerPfo);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::RemoveThreeViewConnectionPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, 
    const MatchedConnectionPathwayMap &threeViewConnectionPathways, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition)
{
    IntVector toRefineIndices;

    for (auto &entry : threeViewConnectionPathways)
    {
        const int mapIndex(entry.first);
        const ProtoShowerVector &protoShowerVector(entry.second);

        /*
        unsigned int refinementCount(0);

        for (const ProtoShower &protoShower : protoShowerVector)
        {
            if (pAlgorithm->IsTrack(protoShower))
                ++refinementCount;
        }

        if (refinementCount >= 2)
        {
        */
            toRefineIndices.push_back(mapIndex);

            for (const ProtoShower &protoShower : protoShowerVector)
            {
                if (!protoShower.m_isHelper)
                {

                    //////////////////////////
                    /*
                    for (const CaloHit *const pCaloHit : protoShower.m_connectionPathway.m_pathwayHitList)
                    {
                        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "HIT TO REMOVE", GREEN, 2);
                    }
                    */
                    //////////////////////////

                    pAlgorithm->RemoveConnectionPathway(pShowerPfo, protoShower);
                }
            }
            //}
    }

    //////////////////////////
    //std::cout << "HITS TO REMOVE FROM THREE VIEW" << std::endl;
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    //////////////////////////

    /*
    CartesianVector showerVertex(0.f, 0.f, 0.f);
    if (this->FindShowerVertex(pAlgorithm, threeViewConnectionPathways, toRefineIndices, nuVertexPosition, pShowerPfo, showerVertex))
        pAlgorithm->SetGammaVertex(showerVertex, pShowerPfo);
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::FindShowerVertex(ShowerStartRefinementAlgorithm *const pAlgorithm, const MatchedConnectionPathwayMap &threeViewConnectionPathways, 
    const IntVector &toRefineIndices, const CartesianVector &nuVertexPosition, const ParticleFlowObject *const pShowerPfo, CartesianVector &showerVertex)
{
    bool found(false);
    float bestDistanceFromNu(std::numeric_limits<float>::max());

    for (auto &entry : threeViewConnectionPathways)
    {
        const int mapIndex(entry.first);

        if (std::find(toRefineIndices.begin(), toRefineIndices.end(), mapIndex) == toRefineIndices.end())
            continue;

        unsigned int viewsWithShowerHits(0);

        const ProtoShowerVector &protoShowerVector(entry.second);
        for (const ProtoShower &protoShower : protoShowerVector)
        {
            if (protoShower.m_showerCore.m_coreHitList.size() < 5)
                break;

            ++viewsWithShowerHits;
        }

        if (viewsWithShowerHits != 3)
            continue;

        CartesianVector showerStart3D(0.f, 0.f, 0.f);
        if (!LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, nuVertexPosition, protoShowerVector[0], protoShowerVector[1], protoShowerVector[2],
            m_maxXSeparation, showerStart3D))
        {
            continue;
        }

        const float distanceFromNu((showerStart3D - nuVertexPosition).GetMagnitude());

        if (distanceFromNu < bestDistanceFromNu)
        {
            found = true;
            bestDistanceFromNu = distanceFromNu;
            showerVertex = showerStart3D;
        }
    }

    //////////////////////////////////////
    /*
    if (found)
    {
        std::cout << "found a new gamma vertex!" << std::endl;
        PfoList visualize;
        visualize.push_back(pShowerPfo);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(pAlgorithm->GetPandora(), &visualize, "SHOWER PFO", VIOLET);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerVertex, "SHOWER VERTEX", VIOLET, 2);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }
    */
    //////////////////////////////////////

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::RemoveOneViewConnectionPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, 
    const MatchedConnectionPathwayMap &oneViewConnectionPathways, const ParticleFlowObject *const pShowerPfo)
{
    for (auto &entry : oneViewConnectionPathways)
    {
        const ProtoShowerVector &protoShowerVector(entry.second);

        for (const ProtoShower &protoShower : protoShowerVector)
        {
            if (!protoShower.m_isHelper) //&& pAlgorithm->IsTrack(protoShower))
            {
                //////////////////////////
                /*
                for (const CaloHit *const pCaloHit : protoShower.m_connectionPathway.m_pathwayHitList)
                {
                    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "HIT TO REMOVE", GREEN, 2);
                }
                */
                //////////////////////////

                pAlgorithm->RemoveConnectionPathway(pShowerPfo, protoShower);
            }
        }
    }

    //////////////////////////
    /*
    std::cout << "hits to remove from one view" << std::endl;
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::AssignShowerHits(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const HitType &hitType, ProtoShowerVector &protoShowerVector)
{
    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

    CaloHitList connectionPathwayHits;
    for (const ProtoShower &protoShower : protoShowerVector)
        connectionPathwayHits.insert(connectionPathwayHits.begin(), protoShower.m_connectionPathway.m_pathwayHitList.begin(), 
            protoShower.m_connectionPathway.m_pathwayHitList.end());

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        if (std::find(connectionPathwayHits.begin(), connectionPathwayHits.end(), pCaloHit) != connectionPathwayHits.end())
            continue;

        int bestProtoShower(-1);
        float bestT(std::numeric_limits<float>::max());

        for (unsigned int index = 0; index < protoShowerVector.size(); ++index)
        {
            const ProtoShower &protoShower(protoShowerVector[index]);
            const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);
            const CartesianVector &showerStartDirection(protoShower.m_showerCore.m_startDirection);

            const float t(showerStartDirection.GetCrossProduct(pCaloHit->GetPositionVector() - showerStartPosition).GetMagnitude());
            const float l(showerStartDirection.GetDotProduct(pCaloHit->GetPositionVector() - showerStartPosition));

            if ((l > 0.f) && (t < m_molliereRadius) && (t < bestT))
            {
                bestT = t;
                bestProtoShower = index;
            }
        }

        if (bestProtoShower >= 0)
            protoShowerVector[bestProtoShower].m_showerCore.m_coreHitList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType tpcView, ProtoShowerVector &protoShowerVector, CaloHitList &usedHitList)
{
    /////////////////////////////////

    std::cout << "Investigating " << (tpcView == TPC_VIEW_U ? "U" : tpcView == TPC_VIEW_V ? "V" : "W") << " view..." << std::endl;
    ClusterList clustersBegin;
    LArPfoHelper::GetClusters(pShowerPfo, tpcView, clustersBegin);
    PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &clustersBegin, "clustersBegin", BLUE);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    /////////////////////////////////

    if (!this->HasPathToNuVertex(pShowerPfo, nuVertexPosition, tpcView))
        return;

    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, tpcView, viewShowerHitList);

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, tpcView));

    /////////////////////////////////
    /*
    PfoList vPfo;
    vPfo.push_back(pShowerPfo);

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - projectedNuVertexPosition);

        if (displacementVector.GetMagnitudeSquared() > (m_pathwaySearchRegion * m_pathwaySearchRegion))
            continue;

        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "Region", RED, 2);
    }

    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, "projectedNuVertexPosition", GREEN, 2);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(pAlgorithm->GetPandora(), &vPfo, "Shower PFo", BLACK);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    AngularDecompositionMap angularDecompositionMap;
    this->FillAngularDecompositionMap(viewShowerHitList, projectedNuVertexPosition, angularDecompositionMap);

    if (angularDecompositionMap.empty())
        return;

    this->SmoothAngularDecompositionMap(angularDecompositionMap);

    // Obtain distribution peaks
    IntVector angularPeakVector;
    this->ObtainPeakVector(angularDecompositionMap, angularPeakVector);

    // Investigate peaks, searching for pathways (peaks indexed by theta0XZ bin lower edge)
    IntVector investigatedPeaks;

    for (unsigned int i = 0; i < angularPeakVector.size(); ++i)
    {
        int bestTheta0XZBin(0);

        if (!this->FindBestAngularPeak(angularDecompositionMap, angularPeakVector, investigatedPeaks, bestTheta0XZBin))
            break;

        investigatedPeaks.push_back(bestTheta0XZBin);

        const float theta0XZ(bestTheta0XZBin * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);

        //////////////////////////////////////
        /*
        std::cout << "bestTheta0XZBin: " << bestTheta0XZBin << std::endl;
        std::cout << "peakDirection: " << peakDirection << std::endl;
        const CartesianVector projection(projectedNuVertexPosition + (peakDirection * 14.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &projection, "direction", BLACK, 2, 1);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        //////////////////////////////////////

        CaloHitList showerSpineHitList;
        try
        {
            this->FindShowerSpine(pAlgorithm, viewShowerHitList, projectedNuVertexPosition, peakDirection, usedHitList, showerSpineHitList);
        }
        catch (...)
        {
            continue;
        }

        // Demand that spine is significant
        if (showerSpineHitList.size() < 20)
        {
            /////////////////////////////////
            /*
            std::cout << "Found shower spine is insignificant" << std::endl;
            */
            /////////////////////////////////
            continue;
        }

        /*
        PandoraMonitoringApi::VisualizeParticleFlowObjects(pAlgorithm->GetPandora(), &vPfo, "Shower PFo", BLACK);
        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "spine", RED, 2);
        }
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */

        // Demand that spine is not the entire pfo
        const float fractionCollected(static_cast<float>(showerSpineHitList.size()) / static_cast<float>(viewShowerHitList.size()));

        if (fractionCollected > 0.8f)
        {
            /////////////////////////////////
            //std::cout << "fraction of collected hits is most of the pfo!! ahh!!" << std::endl;
            /////////////////////////////////
            continue;
        }

        // Obtail longitudinal position of spine hits
        LongitudinalPositionMap longitudinalPositionMap;
        try
        {
            this->ObtainLongitudinalDecomposition(pAlgorithm, showerSpineHitList, longitudinalPositionMap);
        }
        catch(...)
        {
            /////////////////////////////////
            //std::cout << "the super fine fit failed probably a gap in the initial hits??" << std::endl;
            /////////////////////////////////
            continue;
        }

        // Obtain spine energy profile
        EnergySpectrumMap energySpectrumMap;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, longitudinalPositionMap, energySpectrumMap);

        usedHitList.insert(usedHitList.end(), showerSpineHitList.begin(), showerSpineHitList.end());

        // Find shower start position - and pathway!
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        this->FindShowerStart(pAlgorithm, projectedNuVertexPosition, peakDirection, longitudinalPositionMap, energySpectrumMap, showerSpineHitList, 
                              showerStartPosition, viewShowerHitList, isEndDownstream, protoShowerVector, false);

        /////////////////////////////////////
        /*
        for (const CaloHit *const pCaloHit : usedHitList)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "SHOWER SPINE", GREEN, 2);
        }
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        /////////////////////////////////////
    }

    this->FillOutPathways(pAlgorithm, viewShowerHitList, usedHitList, protoShowerVector);

    /////////////////////////////////

    for (const ProtoShower protoShower : protoShowerVector)
    {
        const CaloHitList &connectionPathway(protoShower.m_connectionPathway.m_pathwayHitList);
        const CaloHitList &showerCore(protoShower.m_showerCore.m_coreHitList);
        const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &connectionPathway, "ConnectionPathway", BLACK);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerCore, "ShowerCore", RED);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "ShowerStartPosition", BLACK, 2);
    }

    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, "Projected Nu Vertex", GREEN, 2);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    /////////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::BuildHelperProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType tpcView, ProtoShowerVector &protoShowerVector, CaloHitList &usedHitList)
{
    /////////////////////////////////
    /*
    std::cout << "Investigating " << (tpcView == TPC_VIEW_U ? "U" : tpcView == TPC_VIEW_V ? "V" : "W") << " view..." << std::endl;
    ClusterList clustersBegin;
    LArPfoHelper::GetClusters(pShowerPfo, tpcView, clustersBegin);
    PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &clustersBegin, "clustersBegin", RED);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    // can probably be more efficient here and use the shower start extrema
    const CaloHitList helperHitList(pAlgorithm->GetXIntervalHitsOfType(pShowerPfo, tpcView));

    if (helperHitList.empty())
        return;

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, tpcView));

    AngularDecompositionMap angularDecompositionMap;
    this->FillAngularDecompositionMap(helperHitList, projectedNuVertexPosition, angularDecompositionMap);

    if (angularDecompositionMap.empty())
        return;

    this->SmoothAngularDecompositionMap(angularDecompositionMap);

    // Obtain distribution peaks
    IntVector angularPeakVector;
    this->ObtainPeakVector(angularDecompositionMap, angularPeakVector);

    // Investigate peaks, searching for pathways (peaks indexed by theta0XZ bin lower edge)
    IntVector investigatedPeaks;

    for (unsigned int i = 0; i < angularPeakVector.size(); ++i)
    {
        int bestTheta0XZBin(0);

        if (!this->FindBestAngularPeak(angularDecompositionMap, angularPeakVector, investigatedPeaks, bestTheta0XZBin))
            break;

        investigatedPeaks.push_back(bestTheta0XZBin);

        const float theta0XZ(bestTheta0XZBin * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);

        CaloHitList showerSpineHitList;
        try
        {
            this->FindShowerSpine(pAlgorithm, helperHitList, projectedNuVertexPosition, peakDirection, usedHitList, showerSpineHitList);
        }
        catch (...)
        {
            continue;
        }

        // Demand that spine is significant
        if (showerSpineHitList.size() < 20)
            continue;

        // Obtail longitudinal position of spine hits
        LongitudinalPositionMap longitudinalPositionMap;
        try
        {
            this->ObtainLongitudinalDecomposition(pAlgorithm, showerSpineHitList, longitudinalPositionMap);
        }
        catch(...)
        {
            continue;
        }

        // Obtain spine energy profile
        EnergySpectrumMap energySpectrumMap;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, longitudinalPositionMap, energySpectrumMap);

        usedHitList.insert(usedHitList.end(), showerSpineHitList.begin(), showerSpineHitList.end());

        // Find shower start position - and pathway!
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        this->FindShowerStart(pAlgorithm, projectedNuVertexPosition, peakDirection, longitudinalPositionMap, energySpectrumMap, showerSpineHitList, 
                              showerStartPosition, helperHitList, isEndDownstream, protoShowerVector, true);
    }

    /////////////////////////////////
    /*
    for (const ProtoShower protoShower : protoShowerVector)
    {
        const CaloHitList &connectionPathway(protoShower.m_connectionPathway.m_pathwayHitList);
        const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &connectionPathway, "HELPER ConnectionPathway", GREEN);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, " HELPER ShowerStartPosition", GREEN, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::FillOutPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerPfoHits, 
    CaloHitList &unavailableHits, ProtoShowerVector &protoShowerVector)
{
    for (const CaloHit *const pCaloHit : showerPfoHits)
    {
        /////////////////////////////////
        //const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        /////////////////////////////////

        if (std::find(unavailableHits.begin(), unavailableHits.end(), pCaloHit) != unavailableHits.end())
            continue;

        int bestProtoShower(-1);
        float bestT(std::numeric_limits<float>::max());

        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            if ((std::abs(pMCParticle->GetParticleId()) == 11) || (pMCParticle->GetParticleId() == 22))
                continue;

            for (unsigned int i = 0; i < protoShowerVector.size(); ++i)
            {
                const CartesianVector pathwayDirection((protoShowerVector[i].m_showerCore.m_startPosition - protoShowerVector[i].m_connectionPathway.m_startPosition).GetUnitVector());
                const CartesianVector displacement(pCaloHit->GetPositionVector() - protoShowerVector[i].m_connectionPathway.m_startPosition);
                const float thisT(std::min(protoShowerVector[i].m_connectionPathway.m_startDirection.GetCrossProduct(displacement).GetMagnitude(),
                    LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), protoShowerVector[i].m_connectionPathway.m_pathwayHitList)));
                const float thisL(pathwayDirection.GetDotProduct(displacement));

                if (thisT > 1.f)
                    continue;

                if (thisL < -1.f)
                    continue;

                // only add hits in the initial region of interest
                if (thisL > 10.f)
                    continue;

                if (thisT < bestT)
                {
                    bestT = thisT;
                    bestProtoShower = i;
                }
            }

            if (bestProtoShower >= 0)
            {
                ////////////////////////////////////
                //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "BUFFED OUT HIT", VIOLET, 2);
                ////////////////////////////////////

                protoShowerVector[bestProtoShower].m_connectionPathway.m_pathwayHitList.push_back(pCaloHit);
                unavailableHits.push_back(pCaloHit);
            }
        }
        catch(...)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::IsShowerTruncatable(const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
    const ProtoShowerVector &protoShowerVectorW)
{
    if (protoShowerVectorU.empty() && protoShowerVectorV.empty() && protoShowerVectorW.empty())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::MatchConnectionPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const CartesianVector &nuVertexPosition, 
    const ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
    const ProtoShowerVector &protoShowerVectorW, MatchedConnectionPathwayMap &threeViewConnectionPathways, MatchedConnectionPathwayMap &twoViewConnectionPathways, 
    MatchedConnectionPathwayMap &oneViewConnectionPathways)
{
    // Only match connection pathways in 3D
    IntVector matchedProtoShowersU, matchedProtoShowersV, matchedProtoShowersW;

    // First try to find three view matches
    this->MatchShowerStart(pAlgorithm, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, matchedProtoShowersU, matchedProtoShowersV, 
        matchedProtoShowersW, threeViewConnectionPathways);

    ////////////////////////////
    /*
    std::cout << "THREE VIEW MATCHES - SHOWER START" << std::endl;
    for (auto &entry : threeViewConnectionPathways)
    {
        CartesianVector origin(0.f, 0.f, 0.f);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &origin, "ORIGIN", BLACK, 2);
    
        for (const ProtoShower &protoShower : entry.second)
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShower.m_showerCore.m_startPosition, "THREE MATCH", GREEN, 2);
    
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }
    threeViewConnectionPathways.clear();
    */
    ////////////////////////////

    this->MatchPeakDirection(pAlgorithm, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, matchedProtoShowersU, matchedProtoShowersV, 
      matchedProtoShowersW, threeViewConnectionPathways);

    ////////////////////////////

    std::cout << "THREE VIEW MATCHES" << std::endl;

    for (auto &entry : threeViewConnectionPathways)
    {
        CartesianVector origin(0.f, 0.f, 0.f);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &origin, "ORIGIN", BLACK, 2);

        for (const ProtoShower &protoShower : entry.second)
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShower.m_showerCore.m_startPosition, "THREE MATCH", GREEN, 2);

        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }
    //threeViewConnectionPathways.clear();

    ////////////////////////////

    /* // this can lead to mistakes.. if you get a three view match want to make sure that it is DEFINITELY correct
    this->MatchShowerStart2(pAlgorithm, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, matchedProtoShowersU, matchedProtoShowersV, 
        matchedProtoShowersW, threeViewConnectionPathways);

    ////////////////////////////
    std::cout << "THREE VIEW MATCHES - SHOWER START WITH PROJECTION" << std::endl;

    for (auto &entry : threeViewConnectionPathways)
    {
        CartesianVector origin(0.f, 0.f, 0.f);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &origin, "ORIGIN", BLACK, 2);
    
        for (const ProtoShower &protoShower : entry.second)
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShower.m_showerCore.m_startPosition, "THREE MATCH", GREEN, 2);
    
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }
    */
    // Now add in any unmatched connection pathways
    this->AddUnmatched(protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, matchedProtoShowersU, matchedProtoShowersV, 
        matchedProtoShowersW, oneViewConnectionPathways);

    ////////////////////////////

    std::cout << "ONE VIEW MATCHES" << std::endl;

    for (auto &entry : oneViewConnectionPathways)
    {
        CartesianVector origin(0.f, 0.f, 0.f);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &origin, "ORIGIN", BLACK, 2);

        for (const ProtoShower &protoShower : entry.second)
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShower.m_showerCore.m_startPosition, "ONE MATCH", RED, 2);

        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }

    ////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::MatchShowerStart(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShowerVector &protoShowerVectorU, 
    const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW, IntVector &matchedProtoShowersU, 
    IntVector &matchedProtoShowersV, IntVector &matchedProtoShowersW, MatchedConnectionPathwayMap &matchedConnectionPathwayMap)
{
    bool found(true);

    while(found)
    {
        found = false;

        int bestProtoShowerU(0), bestProtoShowerV(0), bestProtoShowerW(0);
        float bestMetric(std::numeric_limits<float>::max());

        for (unsigned int uIndex = 0; uIndex < protoShowerVectorU.size(); ++uIndex)
        {
            if (std::find(matchedProtoShowersU.begin(), matchedProtoShowersU.end(), uIndex) != matchedProtoShowersU.end())
                continue;

            const ProtoShower &protoShowerU(protoShowerVectorU[uIndex]);

            for (unsigned int vIndex = 0; vIndex < protoShowerVectorV.size(); ++vIndex)
            {
                if (std::find(matchedProtoShowersV.begin(), matchedProtoShowersV.end(), vIndex) != matchedProtoShowersV.end())
                    continue;

                const ProtoShower &protoShowerV(protoShowerVectorV[vIndex]);

                for (unsigned int wIndex = 0; wIndex < protoShowerVectorW.size(); ++wIndex)
                {
                    if (std::find(matchedProtoShowersW.begin(), matchedProtoShowersW.end(), wIndex) != matchedProtoShowersW.end())
                        continue;

                    const ProtoShower &protoShowerW(protoShowerVectorW[wIndex]);

                    if (protoShowerU.m_isHelper && protoShowerV.m_isHelper && protoShowerW.m_isHelper)
                        continue;

                    float metric(std::numeric_limits<float>::max());
                    if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, protoShowerU.m_showerCore.m_startPosition, protoShowerV.m_showerCore.m_startPosition, 
                        protoShowerW.m_showerCore.m_startPosition, m_maxXSeparation, m_maxSeparation, metric))
                        continue;

                    if (!LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, m_maxAngularDeviation))
                        continue;

                    if (metric < bestMetric)
                    {
                        found = true;
                        bestProtoShowerU = uIndex;
                        bestProtoShowerV = vIndex;
                        bestProtoShowerW = wIndex;
                        bestMetric = metric;
                    }
                }
            }
        }

        if (found)
        {
            matchedProtoShowersU.push_back(bestProtoShowerU);
            matchedProtoShowersV.push_back(bestProtoShowerV);
            matchedProtoShowersW.push_back(bestProtoShowerW);

            ProtoShowerVector matchedProtoShowers;
            matchedProtoShowers.push_back(protoShowerVectorU[bestProtoShowerU]);
            matchedProtoShowers.push_back(protoShowerVectorV[bestProtoShowerV]);
            matchedProtoShowers.push_back(protoShowerVectorW[bestProtoShowerW]);

            int mapInt(matchedConnectionPathwayMap.size());
            matchedConnectionPathwayMap[mapInt] = matchedProtoShowers;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::MatchPeakDirection(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShowerVector &protoShowerVectorU, 
    const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW, IntVector &matchedProtoShowersU, IntVector &matchedProtoShowersV, 
    IntVector &matchedProtoShowersW, MatchedConnectionPathwayMap &matchedConnectionPathwayMap)
{
    bool found(true);

    while(found)
    {
        found = false;

        int bestProtoShowerU(0), bestProtoShowerV(0), bestProtoShowerW(0);
        float bestMetric(std::numeric_limits<float>::max());

        for (unsigned int uIndex = 0; uIndex < protoShowerVectorU.size(); ++uIndex)
        {
            if (std::find(matchedProtoShowersU.begin(), matchedProtoShowersU.end(), uIndex) != matchedProtoShowersU.end())
                continue;

            const ProtoShower &protoShowerU(protoShowerVectorU[uIndex]);

            for (unsigned int vIndex = 0; vIndex < protoShowerVectorV.size(); ++vIndex)
            {
                if (std::find(matchedProtoShowersV.begin(), matchedProtoShowersV.end(), vIndex) != matchedProtoShowersV.end())
                    continue;

                const ProtoShower &protoShowerV(protoShowerVectorV[vIndex]);

                for (unsigned int wIndex = 0; wIndex < protoShowerVectorW.size(); ++wIndex)
                {
                    if (std::find(matchedProtoShowersW.begin(), matchedProtoShowersW.end(), wIndex) != matchedProtoShowersW.end())
                        continue;

                    const ProtoShower &protoShowerW(protoShowerVectorW[wIndex]);

                    if (protoShowerU.m_isHelper && protoShowerV.m_isHelper && protoShowerW.m_isHelper)
                        continue;

                    float metric(std::numeric_limits<float>::max());
                    if (!LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, protoShowerU.m_connectionPathway.m_startDirection, 
                        protoShowerV.m_connectionPathway.m_startDirection, protoShowerW.m_connectionPathway.m_startDirection, m_maxAngularDeviation, metric))
                    {
                        continue;
                    }

                    if (metric < bestMetric)
                    {
                        found = true;
                        bestProtoShowerU = uIndex;
                        bestProtoShowerV = vIndex;
                        bestProtoShowerW = wIndex;
                        bestMetric = metric;
                    }
                }
            }
        }

        if (found)
        {
            matchedProtoShowersU.push_back(bestProtoShowerU);
            matchedProtoShowersV.push_back(bestProtoShowerV);
            matchedProtoShowersW.push_back(bestProtoShowerW);

            ProtoShowerVector matchedProtoShowers;
            matchedProtoShowers.push_back(protoShowerVectorU[bestProtoShowerU]);
            matchedProtoShowers.push_back(protoShowerVectorV[bestProtoShowerV]);
            matchedProtoShowers.push_back(protoShowerVectorW[bestProtoShowerW]);

            int mapInt(matchedConnectionPathwayMap.size());
            matchedConnectionPathwayMap[mapInt] = matchedProtoShowers;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::AddUnmatched(const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, 
    const ProtoShowerVector &protoShowerVectorW, IntVector &matchedProtoShowersU, IntVector &matchedProtoShowersV, IntVector &matchedProtoShowersW, 
    MatchedConnectionPathwayMap &matchedConnectionPathwayMap)
{
    this->AddUnmatched(protoShowerVectorU, matchedProtoShowersU, matchedConnectionPathwayMap);
    this->AddUnmatched(protoShowerVectorV, matchedProtoShowersV, matchedConnectionPathwayMap);
    this->AddUnmatched(protoShowerVectorW, matchedProtoShowersW, matchedConnectionPathwayMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::AddUnmatched(const ProtoShowerVector &protoShowerVector, IntVector &matchedProtoShowers, 
    MatchedConnectionPathwayMap &matchedConnectionPathwayMap)
{
    for (unsigned int i = 0; i < protoShowerVector.size(); ++i)
    {
        if (std::find(matchedProtoShowers.begin(), matchedProtoShowers.end(), i) != matchedProtoShowers.end())
            continue;

        if (protoShowerVector[i].m_isHelper)
            continue;

        matchedProtoShowers.push_back(i);

        ProtoShowerVector unmatchedProtoShowers;
        unmatchedProtoShowers.push_back(protoShowerVector[i]);

        const int mapInt(matchedConnectionPathwayMap.size());
        matchedConnectionPathwayMap[mapInt] = unmatchedProtoShowers;
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxSeparation", m_maxSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxAngularDeviation", m_maxAngularDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxProjectedShowerStartSeparation", m_maxProjectedShowerStartSeparation));

    ShowerStartRefinementBaseTool::ReadSettings(xmlHandle);

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::FillTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition)
{
    CaloHitList caloHits2D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_U, caloHits2D);
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_V, caloHits2D);
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_W, caloHits2D);

    const CartesianVector nuVertexPositionU(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_U));
    const CartesianVector nuVertexPositionV(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_V));
    const CartesianVector nuVertexPositionW(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_W));

    ++m_counter;

    FloatVector deviationAngleU, deviationAngleV, deviationAngleW;
    FloatVector xCoordinateU, zCoordinateU;
    FloatVector xCoordinateV, zCoordinateV;
    FloatVector xCoordinateW, zCoordinateW;
    FloatVector hitEnergyU, hitEnergyV, hitEnergyW;

    CartesianPointVector initialHitPositionsU, initialHitPositionsV, initialHitPositionsW;

    for (const CaloHit *const pCaloHit : caloHits2D)
    {
        const HitType hitType(pCaloHit->GetHitType());
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector nuVertexPositionProjection(hitType == TPC_VIEW_U ? nuVertexPositionU : hitType == TPC_VIEW_V ? nuVertexPositionV : nuVertexPositionW);
        CartesianPointVector &initialHitPositions(hitType == TPC_VIEW_U ? initialHitPositionsU : hitType == TPC_VIEW_V ? initialHitPositionsV : initialHitPositionsW);

        if ((hitPosition - nuVertexPositionProjection).GetMagnitude() > 14.f)
            continue;

        initialHitPositions.push_back(pCaloHit->GetPositionVector());
    }

    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : caloHits2D)
    {
        const HitType hitType(pCaloHit->GetHitType());
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector nuVertexPositionProjection(hitType == TPC_VIEW_U ? nuVertexPositionU : hitType == TPC_VIEW_V ? nuVertexPositionV : nuVertexPositionW);

        if ((hitPosition - nuVertexPositionProjection).GetMagnitude() > 14.f)
            continue;

        FloatVector &hitEnergy(hitType == TPC_VIEW_U ? hitEnergyU : hitType == TPC_VIEW_V ? hitEnergyV : hitEnergyW);

        hitEnergy.push_back(pCaloHit->GetElectromagneticEnergy());

        FloatVector &xCoordinate(hitType == TPC_VIEW_U ? xCoordinateU : hitType == TPC_VIEW_V ? xCoordinateV : xCoordinateW);
        FloatVector &zCoordinate(hitType == TPC_VIEW_U ? zCoordinateU : hitType == TPC_VIEW_V ? zCoordinateV : zCoordinateW);

        xCoordinate.push_back(hitPosition.GetX());
        zCoordinate.push_back(hitPosition.GetZ());

        FloatVector &deviationAngle(hitType == TPC_VIEW_U ? deviationAngleU : hitType == TPC_VIEW_V ? deviationAngleV : deviationAngleW);

        float deviationAngleTemp = xAxis.GetOpeningAngle(hitPosition - nuVertexPositionProjection);

        if ((hitPosition.GetZ() - nuVertexPositionProjection.GetZ()) < 0.f)
            deviationAngleTemp *= -1.f;

        deviationAngle.push_back(deviationAngleTemp);
    }

    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ShowerCounter", m_counter));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleU", &deviationAngleU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleV", &deviationAngleV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleW", &deviationAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "XCoordinateU", &xCoordinateU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ZCoordinateU", &zCoordinateU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "XCoordinateV", &xCoordinateV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ZCoordinateV", &zCoordinateV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "XCoordinateW", &xCoordinateW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ZCoordinateW", &zCoordinateW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyU", &hitEnergyU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyV", &hitEnergyV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyW", &hitEnergyW));
    PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "ShowerDistribution"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::AngularDistributionTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const AngularDecompositionMap &angularDecompositionMapU, 
    const AngularDecompositionMap &angularDecompositionMapV, const AngularDecompositionMap &angularDecompositionMapW)
{
    IntVector deviationBinU, deviationBinV, deviationBinW;
    FloatVector deviationWeightU, deviationWeightV, deviationWeightW;

    for (const auto &entry : angularDecompositionMapU)
    {
        deviationBinU.push_back(entry.first);
        deviationWeightU.push_back(entry.second);
    }

    for (const auto &entry : angularDecompositionMapV)
    {
        deviationBinV.push_back(entry.first);
        deviationWeightV.push_back(entry.second);
    }

    for (const auto &entry : angularDecompositionMapW)
    {
        deviationBinW.push_back(entry.first);
        deviationWeightW.push_back(entry.second);
    }

    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationBinU", &deviationBinU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationBinV", &deviationBinV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationBinW", &deviationBinW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationWeightU", &deviationWeightU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationWeightV", &deviationWeightV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationWeightW", &deviationWeightW));
    PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "DeviationAngle"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content



    // Now try to find two view matches <-- i think when we reach this point the shower starts are just too crap...?
    /*
    this->MatchShowerStartTwoViews(pAlgorithm, pShowerPfo, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, matchedProtoShowersU, 
        matchedProtoShowersV, matchedProtoShowersW, twoViewConnectionPathways);
    
    ////////////////////////////
    std::cout << "TWO VIEW MATCHES - X OF SHOWER START" << std::endl;

    for (auto &entry : twoViewConnectionPathways)
    {
        CartesianVector origin(0.f, 0.f, 0.f);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &origin, "ORIGIN", BLACK, 2);

        for (const ProtoShower &protoShower : entry.second)
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShower.m_showerCore.m_startPosition, "TWO MATCH", VIOLET, 2);

        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }
    twoViewConnectionPathways.clear();
    ////////////////////////////
    
    this->MatchPeakDirectionTwoViews(pAlgorithm, nuVertexPosition, pShowerPfo, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, 
        matchedProtoShowersU, matchedProtoShowersV, matchedProtoShowersW, twoViewConnectionPathways);

    ////////////////////////////
    std::cout << "TWO VIEW MATCHES" << std::endl;

    for (auto &entry : twoViewConnectionPathways)
    {
        CartesianVector origin(0.f, 0.f, 0.f);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &origin, "ORIGIN", BLACK, 2);

        for (const ProtoShower &protoShower : entry.second)
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShower.m_showerCore.m_startPosition, "TWO MATCH", VIOLET, 2);

        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }
    */
    ////////////////////////////



//------------------------------------------------------------------------------------------------------------------------------------------
/*
void GammaStartRefinementTool::MatchShowerStartTwoViews(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW, IntVector &matchedProtoShowersU, 
    IntVector &matchedProtoShowersV, IntVector &matchedProtoShowersW, MatchedConnectionPathwayMap &matchedConnectionPathwayMap)
{
    bool found(true);

    while (found)
    {
        found = false;

        int bestProtoShower1(0), bestProtoShower2(0);
        HitType bestHitType1(TPC_VIEW_U), bestHitType2(TPC_VIEW_U);
        float bestSeparation(std::numeric_limits<float>::max());

        if (this->MatchShowerStartTwoViews(pAlgorithm, pShowerPfo, protoShowerVectorU, protoShowerVectorV, matchedProtoShowersU, matchedProtoShowersV, 
            bestProtoShower1, bestProtoShower2, bestSeparation))
        {
            found = true;
            bestHitType1 = TPC_VIEW_U; bestHitType2 = TPC_VIEW_V;
        }

        if (this->MatchShowerStartTwoViews(pAlgorithm, pShowerPfo, protoShowerVectorU, protoShowerVectorW, matchedProtoShowersU, matchedProtoShowersW, 
            bestProtoShower1, bestProtoShower2, bestSeparation))
        {
            found = true;
            bestHitType1 = TPC_VIEW_U; bestHitType2 = TPC_VIEW_W;
        }

        if (this->MatchShowerStartTwoViews(pAlgorithm, pShowerPfo, protoShowerVectorV, protoShowerVectorW, matchedProtoShowersV, matchedProtoShowersW, 
            bestProtoShower1, bestProtoShower2, bestSeparation))
        {
            found = true;
            bestHitType1 = TPC_VIEW_V; bestHitType2 = TPC_VIEW_W;
        }

        if (found)
        {
            ProtoShowerVector matchedProtoShowers;
            matchedProtoShowers.push_back(bestHitType1 == TPC_VIEW_U ? protoShowerVectorU[bestProtoShower1] : bestHitType1 == TPC_VIEW_V ? 
                protoShowerVectorV[bestProtoShower1] : protoShowerVectorW[bestProtoShower1]);
            matchedProtoShowers.push_back(bestHitType2 == TPC_VIEW_U ? protoShowerVectorU[bestProtoShower2] : bestHitType2 == TPC_VIEW_V ? 
                protoShowerVectorV[bestProtoShower2] : protoShowerVectorW[bestProtoShower2]);

            IntVector &matchedProtoShowers1(bestHitType1 == TPC_VIEW_U ? matchedProtoShowersU : bestHitType1 == TPC_VIEW_V ? matchedProtoShowersV : matchedProtoShowersW);
            IntVector &matchedProtoShowers2(bestHitType2 == TPC_VIEW_U ? matchedProtoShowersU : bestHitType2 == TPC_VIEW_V ? matchedProtoShowersV : matchedProtoShowersW);

            matchedProtoShowers1.push_back(bestProtoShower1);
            matchedProtoShowers2.push_back(bestProtoShower2);

            int mapInt(matchedConnectionPathwayMap.size());
            matchedConnectionPathwayMap[mapInt] = matchedProtoShowers;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// maybe should be done with the hits in the connection pathways?
bool GammaStartRefinementTool::MatchShowerStartTwoViews(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShowerVector &protoShowerVector1, const ProtoShowerVector &protoShowerVector2, IntVector &matchedProtoShowers1, IntVector &matchedProtoShowers2, 
    int &bestProtoShower1, int &bestProtoShower2, float &bestSeparation)
{
    bool found(false);

    for (unsigned int index1 = 0; index1 < protoShowerVector1.size(); ++index1)
    {
        if (std::find(matchedProtoShowers1.begin(), matchedProtoShowers1.end(), index1) != matchedProtoShowers1.end())
            continue;

        const ProtoShower &protoShower1(protoShowerVector1[index1]);
        const CartesianVector &showerStart1(protoShower1.m_showerCore.m_startPosition);
        //const HitType hitType1(protoShower1.m_connectionPathway.m_pathwayHitList.front()->GetHitType());

        for (unsigned int index2 = 0; index2 < protoShowerVector2.size(); ++index2)
        {
            if (std::find(matchedProtoShowers2.begin(), matchedProtoShowers2.end(), index2) != matchedProtoShowers2.end())
                continue;

            const ProtoShower &protoShower2(protoShowerVector2[index2]);
            const CartesianVector &showerStart2(protoShower2.m_showerCore.m_startPosition);
            //const HitType hitType2(protoShower2.m_connectionPathway.m_pathwayHitList.front()->GetHitType());
            
            if (protoShower1.m_isHelper && protoShower2.m_isHelper)
                continue;

            // Now assess compatability
            const float xSeparation(std::fabs(showerStart1.GetX() - showerStart2.GetX()));

            if (xSeparation > m_maxXSeparation)
                continue;

            // Check projection in third view
            HitType hitType3(TPC_VIEW_U);

            if (hitType1 == TPC_VIEW_U)
                hitType3 = (hitType2 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_V); 
            else if (hitType1 == TPC_VIEW_V)
                hitType3 = (hitType2 == TPC_VIEW_W ? TPC_VIEW_U : TPC_VIEW_W); 
            else
                hitType3 = (hitType2 == TPC_VIEW_U ? TPC_VIEW_V : TPC_VIEW_U); 

            ClusterList cluster3;
            LArPfoHelper::GetClusters(pShowerPfo, hitType3, cluster3);

            float chi2(0.f);
            CartesianVector projection3(0.f, 0.f, 0.f);

            LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType1, hitType2, showerStart1, showerStart2, projection3, chi2);
            
            const float closestDistance(LArClusterHelper::GetClosestDistance(projection3, cluster3));

            //////////////////////////////
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStart1, "showerStart1", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStart2, "showerStart2", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projection3, "PROJECTION", RED, 2);
            std::cout << "closestDistance: " << closestDistance << std::endl;
            PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            //////////////////////////////

            if (closestDistance > m_maxProjectedShowerStartSeparation)
                continue;

            if (xSeparation < bestSeparation)
            {
                found = true;
                bestProtoShower1 = index1;
                bestProtoShower2 = index2;
                bestSeparation = xSeparation;
            }
        }
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// maybe should be done with the hits in the connection pathways?
void GammaStartRefinementTool::MatchPeakDirectionTwoViews(ShowerStartRefinementAlgorithm *const pAlgorithm, const CartesianVector &nuVertexPosition, 
    const ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW, 
    IntVector &matchedProtoShowersU, IntVector &matchedProtoShowersV, IntVector &matchedProtoShowersW, MatchedConnectionPathwayMap &matchedConnectionPathwayMap)
{
    bool found(true);

    while (found)
    {
        found = false;

        // Use distance from nu vertex as a tie-breaker
        int bestProtoShower1(0), bestProtoShower2(0);
        HitType bestHitType1(TPC_VIEW_U), bestHitType2(TPC_VIEW_U);
        unsigned int bestHitCount(0); float bestClosestDistance(std::numeric_limits<float>::max()); 

        if (this->MatchPeakDirectionTwoViews(pAlgorithm, nuVertexPosition, pShowerPfo, protoShowerVectorU, protoShowerVectorV, matchedProtoShowersU, 
            matchedProtoShowersV, bestProtoShower1, bestProtoShower2, bestHitCount, bestClosestDistance))
        {
            found = true;
            bestHitType1 = TPC_VIEW_U; bestHitType2 = TPC_VIEW_V;
        }

        if (this->MatchPeakDirectionTwoViews(pAlgorithm, nuVertexPosition, pShowerPfo, protoShowerVectorU, protoShowerVectorW, matchedProtoShowersU, 
            matchedProtoShowersW, bestProtoShower1, bestProtoShower2, bestHitCount, bestClosestDistance))
        {
            found = true;
            bestHitType1 = TPC_VIEW_U; bestHitType2 = TPC_VIEW_W;
        }

        if (this->MatchPeakDirectionTwoViews(pAlgorithm, nuVertexPosition, pShowerPfo, protoShowerVectorV, protoShowerVectorW, matchedProtoShowersV, 
            matchedProtoShowersW, bestProtoShower1, bestProtoShower2, bestHitCount, bestClosestDistance))
        {
            found = true;
            bestHitType1 = TPC_VIEW_V; bestHitType2 = TPC_VIEW_W;
        }

        if (found)
        {
            ProtoShowerVector matchedProtoShowers;
            matchedProtoShowers.push_back(bestHitType1 == TPC_VIEW_U ? protoShowerVectorU[bestProtoShower1] : bestHitType1 == TPC_VIEW_V ? 
                protoShowerVectorV[bestProtoShower1] : protoShowerVectorW[bestProtoShower1]);
            matchedProtoShowers.push_back(bestHitType2 == TPC_VIEW_U ? protoShowerVectorU[bestProtoShower2] : bestHitType2 == TPC_VIEW_V ? 
                protoShowerVectorV[bestProtoShower2] : protoShowerVectorW[bestProtoShower2]);

            IntVector &matchedProtoShowers1(bestHitType1 == TPC_VIEW_U ? matchedProtoShowersU : bestHitType1 == TPC_VIEW_V ? matchedProtoShowersV : matchedProtoShowersW);
            IntVector &matchedProtoShowers2(bestHitType2 == TPC_VIEW_U ? matchedProtoShowersU : bestHitType2 == TPC_VIEW_V ? matchedProtoShowersV : matchedProtoShowersW);

            matchedProtoShowers1.push_back(bestProtoShower1);
            matchedProtoShowers2.push_back(bestProtoShower2);

            int mapInt(matchedConnectionPathwayMap.size());
            matchedConnectionPathwayMap[mapInt] = matchedProtoShowers;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::MatchPeakDirectionTwoViews(ShowerStartRefinementAlgorithm *const pAlgorithm, const CartesianVector &nuVertexPosition, 
    const ParticleFlowObject *const pShowerPfo, const ProtoShowerVector &protoShowerVector1, const ProtoShowerVector &protoShowerVector2, 
    IntVector &matchedProtoShowers1, IntVector &matchedProtoShowers2, int &bestProtoShower1, int &bestProtoShower2, unsigned int &bestProjectionHitCount, float &bestClosestDistance)
{
    bool found(false);

    for (unsigned int index1 = 0; index1 < protoShowerVector1.size(); ++index1)
    {
        if (std::find(matchedProtoShowers1.begin(), matchedProtoShowers1.end(), index1) != matchedProtoShowers1.end())
            continue;

        const ProtoShower &protoShower1(protoShowerVector1[index1]);
        const CartesianVector &showerPeakDirection1(protoShower1.m_connectionPathway.m_startDirection);
        const HitType hitType1(protoShower1.m_connectionPathway.m_pathwayHitList.front()->GetHitType());

        for (unsigned int index2 = 0; index2 < protoShowerVector2.size(); ++index2)
        {
            if (std::find(matchedProtoShowers2.begin(), matchedProtoShowers2.end(), index2) != matchedProtoShowers2.end())
                continue;

            const ProtoShower &protoShower2(protoShowerVector2[index2]);
            const CartesianVector &showerPeakDirection2(protoShower2.m_connectionPathway.m_startDirection);
            const HitType hitType2(protoShower2.m_connectionPathway.m_pathwayHitList.front()->GetHitType());
            
            if (protoShower1.m_isHelper && protoShower2.m_isHelper)
                continue;

            if (showerPeakDirection1.GetX() * showerPeakDirection2.GetX() < 0.f)
                continue;

            // Check projection in third view
            HitType hitType3(TPC_VIEW_U);

            if (hitType1 == TPC_VIEW_U)
                hitType3 = (hitType2 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_V); 
            else if (hitType1 == TPC_VIEW_V)
                hitType3 = (hitType2 == TPC_VIEW_W ? TPC_VIEW_U : TPC_VIEW_W); 
            else
                hitType3 = (hitType2 == TPC_VIEW_U ? TPC_VIEW_V : TPC_VIEW_U); 

            CaloHitList showerHitList3;
            LArPfoHelper::GetCaloHits(pShowerPfo, hitType3, showerHitList3);

            const CartesianVector projection3(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), hitType1, hitType2, showerPeakDirection1, showerPeakDirection2));
            const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, hitType3));

            // See if anything is there...
            CaloHitList projectionHits;
            for (const CaloHit *const pCaloHit : showerHitList3)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                const CartesianVector &displacementVector(hitPosition - projectedNuVertexPosition);
                const float l(projection3.GetDotProduct(displacementVector));
                const float t(projection3.GetCrossProduct(displacementVector).GetMagnitude());

                if ((l < m_growingFitInitialLength) && (l > 0.f) && (t < m_initialFitDistanceToLine))
                    projectionHits.push_back(pCaloHit);
            }


            ////////////////////////////////////////
            const CartesianVector &nuVertex1(protoShower1.m_connectionPathway.m_startPosition);
            const CartesianVector &nuVertex2(protoShower2.m_connectionPathway.m_startPosition);
            const CartesianVector &nuVertex3(projectedNuVertexPosition);
            const CartesianVector end1(nuVertex1 + (showerPeakDirection1 * 10.f));
            const CartesianVector end2(nuVertex2 + (showerPeakDirection2 * 10.f));
            const CartesianVector end3(nuVertex3 + (projection3 * 10.f));
            PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertex3, &end3, "PROJECTION", RED, 2, 2);
            PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertex1, &end1, "DIRECTION", BLACK, 2, 2);
            PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertex2, &end2, "DIRECTION", BLACK, 2, 2);

            for (const CaloHit *const pCaloHit : projectionHits)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "found hit", VIOLET, 2);
            }
            std::cout << "projectionHits.size(): " << projectionHits.size() << std::endl;

            PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ////////////////////////////////////////

            if (projectionHits.size() < m_minInitialHitsFound)
                continue;

            const float closestDistance(LArClusterHelper::GetClosestDistance(projectedNuVertexPosition, projectionHits));

            if ((projectionHits.size() > bestProjectionHitCount) || ((projectionHits.size() == bestProjectionHitCount) && (closestDistance < bestClosestDistance)))
            {
                found = true;
                bestProtoShower1 = index1;
                bestProtoShower2 = index2;
                bestProjectionHitCount = projectionHits.size();
                bestClosestDistance = closestDistance;
            }
        }
    }

    return found;
}
*/
