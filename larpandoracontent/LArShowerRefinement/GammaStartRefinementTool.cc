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
    m_showerCounter(0), 
    m_trackSearchWindow(5.f),
    m_minTrackBlipMean(3.f),
    m_electronFraction(0.3f)
{
}

GammaStartRefinementTool::~GammaStartRefinementTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition)
{
    // only apply gamma refinement algorithm to showers
    if (!LArPfoHelper::IsShower(pShowerPfo))
        return false;

    if (!this->HasPathToNuVertex(pShowerPfo, nuVertexPosition))
        return false;

    ////////////////////////////////
    // temporary  - actually works really well so keep this
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

    if (caloHits3D.size() < 100)
        return false;
    ////////////////////////////////

    this->FillTree(pAlgorithm, pShowerPfo, nuVertexPosition);
    //return true;

    ProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    /////////////////////////////////
    std::cout << "CONNECTION PATHWAYS U" << std::endl;

    ClusterList clustersBeginU;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_U, clustersBeginU);
    PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &clustersBeginU, "clustersBeginU", RED);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////

    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, protoShowerVectorU);

    /////////////////////////////////
    for (const ProtoShower protoShower : protoShowerVectorU)
    {
        const CaloHitList &connectionPathway(protoShower.m_connectionPathway.m_pathwayHitList);
        const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &connectionPathway, "ConnectionPathway", RED);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "ShowerStartPosition", RED, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////

    this->RemoveTrackPathways(pAlgorithm, pShowerPfo, protoShowerVectorU);

    /////////////////////////////////
    ClusterList clustersAfterU;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_U, clustersAfterU);
    PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &clustersAfterU, "clustersAfterU", BLUE);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////

    /////////////////////////////////
    std::cout << "CONNECTION PATHWAYS V" << std::endl;

    ClusterList clustersBeginV;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_V, clustersBeginV);
    PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &clustersBeginV, "clustersBeginV", RED);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////

    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, protoShowerVectorV);

    /////////////////////////////////

    for (const ProtoShower protoShower : protoShowerVectorV)
    {
        const CaloHitList &connectionPathway(protoShower.m_connectionPathway.m_pathwayHitList);
        const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &connectionPathway, "ConnectionPathway", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "ShowerStartPosition", BLACK, 2);
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    /////////////////////////////////

    this->RemoveTrackPathways(pAlgorithm, pShowerPfo, protoShowerVectorV);

    /////////////////////////////////

    PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &clustersBeginV, "clustersAfterV", BLUE);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    /////////////////////////////////

    /////////////////////////////////

    std::cout << "CONNECTION PATHWAYS W" << std::endl;

    ClusterList clustersBeginW;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_W, clustersBeginW);
    PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &clustersBeginW, "clustersBeginW", RED);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    /////////////////////////////////

    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, protoShowerVectorW);

    /////////////////////////////////

    for (const ProtoShower protoShower : protoShowerVectorW)
    {
        const CaloHitList &connectionPathway(protoShower.m_connectionPathway.m_pathwayHitList);
        const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &connectionPathway, "ConnectionPathway", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "ShowerStartPosition", BLACK, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    /////////////////////////////////

    this->RemoveTrackPathways(pAlgorithm, pShowerPfo, protoShowerVectorW);

    /////////////////////////////////

    PandoraMonitoringApi::VisualizeClusters(pAlgorithm->GetPandora(), &clustersBeginW, "clustersAfterW", BLUE);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    /////////////////////////////////

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType tpcView, ProtoShowerVector &protoShowerVector)
{
    /////////////////////////////////
    //std::cout << "Investigating " << (tpcView == TPC_VIEW_U ? "U" : tpcView == TPC_VIEW_V ? "V" : "W") << " view..." << std::endl;
    /////////////////////////////////

    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, tpcView, viewShowerHitList);

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, tpcView));

    // Fill angular decomposition map


    //////////////////////
    
    PfoList vPfo;
    vPfo.push_back(pShowerPfo);
    /*
    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - projectedNuVertexPosition);

        if (displacementVector.GetMagnitudeSquared() > (m_pathwaySearchRegion * m_pathwaySearchRegion))
            continue;

        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "Region", RED, 2);
    }

    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, "projectedNuVertexPosition", GREEN, 2);
    */
    PandoraMonitoringApi::VisualizeParticleFlowObjects(pAlgorithm->GetPandora(), &vPfo, "Shower PFo", BLACK);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////

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

    // Keep track of 'taken' hits
    CaloHitList unavailableHits;

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
        //std::cout << "bestTheta0XZBin: " << bestTheta0XZBin << std::endl;
        //std::cout << "peakDirection: " << peakDirection << std::endl;
        //const CartesianVector projection(projectedNuVertexPosition + (peakDirection * 14.f));
        //PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &projection, "direction", BLACK, 2, 1);
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        //////////////////////////////////////

        CaloHitList showerSpineHitList;
        try
        {
            this->FindShowerSpine(pAlgorithm, viewShowerHitList, projectedNuVertexPosition, peakDirection, unavailableHits, showerSpineHitList);
        }
        catch (...)
        {
            continue;
        }

        // Demand that spine is significant
        if (showerSpineHitList.size() < 20)
        {
            /////////////////////////////////
            //std::cout << "Found shower spine is insignificant" << std::endl;
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

        unavailableHits.insert(unavailableHits.end(), showerSpineHitList.begin(), showerSpineHitList.end());

        // Find shower start position - and pathway! - and possible post shower hits? (ISOBEL TODO)
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        this->FindShowerStart(pAlgorithm, projectedNuVertexPosition, peakDirection, longitudinalPositionMap, energySpectrumMap, showerSpineHitList, 
            showerStartPosition, viewShowerHitList, isEndDownstream, protoShowerVector);

        /////////////////////////////////////
        /*
        for (const CaloHit *const pCaloHit : unavailableHits)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "SHOWER SPINE", GREEN, 2);
        }
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        /////////////////////////////////////
    }

    this->FillOutPathways(pAlgorithm, viewShowerHitList, unavailableHits, protoShowerVector);

    /////////////////////////////////////
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------











//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::FillMCParticleToHitMap(const CaloHitList *const pCaloHitList, const HitType tpcView, LArMCParticleHelper::MCContributionMap &mcParticleToHitMap)
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (pCaloHit->GetHitType() != tpcView)
            continue;

        try 
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            mcParticleToHitMap[pMCParticle].push_back(pCaloHit);
        }
        catch(...)
        {
            continue;
        }
    }
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

void GammaStartRefinementTool::RemoveTrackPathways(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShowerVector &protoShowerVector)
{
    for (const ProtoShower &protoShower : protoShowerVector)
    {
        if (this->IsTrack(protoShower))
            this->RemoveConnectionPathway(pAlgorithm, pShowerPfo, protoShower);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::IsTrack(const ProtoShower &protoShower)
{
    int showerHits(0);
    int electronHits(0);
    int trackHits(0);

    for (const CaloHit *const pCaloHit : protoShower.m_connectionPathway.m_pathwayHitList)
    {
        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            const int pdg(std::abs(pMCParticle->GetParticleId()));

            if ((pdg == 11) || (pdg == 22))
            {
                ++showerHits;

                if (pdg == 11)
                    ++electronHits;
            }
            else
            {
                ++trackHits;
            }
        }
        catch(...)
        {
            continue;
        }
    }

    const float electronFraction(static_cast<float>(electronHits) / static_cast<float>(showerHits));
    const bool isElectron(electronFraction > m_electronFraction);

    if (isElectron)
        return false;

    ////////////////////////////////////
    //std::cout << "fraction: " << (static_cast<float>(showerHits) / static_cast<float>(showerHits + trackHits)) << std::endl;
    ////////////////////////////////////

    const float showerFraction(static_cast<float>(showerHits) / static_cast<float>(showerHits + trackHits));

    return (showerFraction < 0.8f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::RemoveConnectionPathway(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const ProtoShower &protoShower)
{
    const HitType tpcView(protoShower.m_connectionPathway.m_pathwayHitList.front()->GetHitType());

    ClusterList clusterList;
    LArPfoHelper::GetClusters(pShowerPfo, tpcView, clusterList);

    std::string clusterListName(tpcView == TPC_VIEW_U ? "ClustersU" : tpcView == TPC_VIEW_V ? "ClustersV" : "ClustersW");
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, clusterListName));

    for (const CaloHit *const pCaloHit : protoShower.m_connectionPathway.m_pathwayHitList)
    {
        ////////////////////////////////////
        //std::cout << "REMOVE HIT" << std::endl;
        ////////////////////////////////////

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*pAlgorithm, clusterList.front(), pCaloHit));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrackSearchWindow", m_trackSearchWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinTrackBlipMean", m_minTrackBlipMean));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ElectronFraction", m_electronFraction));

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
