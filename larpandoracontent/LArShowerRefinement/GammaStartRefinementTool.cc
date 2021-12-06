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
    m_microSlidingFitWindow(20),
    m_minSigmaDeviation(5.f),
    m_trackSearchWindow(5.f),
    m_nInitialEnergyBins(5),
    m_minTrackBlipMean(3.f),
    m_showerSlidingFitWindow(20),
    m_molliereRadius(9.f)
{
}

GammaStartRefinementTool::~GammaStartRefinementTool()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "DeviationAngle", "ShowerDistribution.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "BAD JAM" << std::endl;
    }
    
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ShowerDistribution", "ShowerDistribution.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "BAD JAM" << std::endl;
    }

    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "EnergyDistributionU", "ShowerDistribution.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "BAD JAM" << std::endl;
    }

    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "EnergyDistributionV", "ShowerDistribution.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "BAD JAM" << std::endl;
    }

    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "EnergyDistributionW", "ShowerDistribution.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "BAD JAM" << std::endl;
    }
    
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
    // temporary 
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

    if (caloHits3D.size() < 100)
        return false;
    ////////////////////////////////

    //this->FillTree(pAlgorithm, pShowerPfo, nuVertexPosition);

    // Investigate U View
    std::cout << "Investigating U view..." << std::endl;

    DeviationAngleMap deviationAngleMapU;
    this->FillAngularDecompositionMap(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, deviationAngleMapU);

    IntVector angularPeakVectorU;
    this->ObtainViewPeakVector(deviationAngleMapU, angularPeakVectorU);

    IntVector showerCounterVectorU;
    FloatVector projectionVectorU, energiesU;

    IntVector investigatedPeaksU;
    
    for (unsigned int i = 0; i < angularPeakVectorU.size(); ++i)
    {
        int bestAngularPeak(0);
        if (!this->FindBestAngularPeak(deviationAngleMapU, angularPeakVectorU, investigatedPeaksU, bestAngularPeak))
            break;

        investigatedPeaksU.push_back(bestAngularPeak);

        const CartesianVector peakDirection(std::cos(bestAngularPeak * pAlgorithm->m_binSize), 0.f, std::sin(bestAngularPeak * pAlgorithm->m_binSize));
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);
        std::cout << "bestAngularPeak U view: " << bestAngularPeak << std::endl;
        std::cout << "peakDirection: " << peakDirection << std::endl;

        const CartesianVector nuVertexPositionU(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_U));
        const CartesianVector projection(nuVertexPositionU + (peakDirection * 14.f));

        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexPositionU, &projection, "direction", BLACK, 2, 1);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        CaloHitList showerSpineHitList;
        this->FindShowerSpine(pAlgorithm, pShowerPfo, nuVertexPositionU, peakDirection, TPC_VIEW_U, showerSpineHitList);

        if (showerSpineHitList.size() < 20)
            continue;

        //look at energy distribution
        EnergySpectrumMap energySpectrumMapU;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, showerCounterVectorU, projectionVectorU, energiesU, energySpectrumMapU);

        CaloHitList showerPfoHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_U, showerPfoHitList);

        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        if (!this->FindShowerStart(pAlgorithm, energySpectrumMapU, showerSpineHitList, showerStartPosition, showerPfoHitList, isEndDownstream))
        {
            std::cout << "no shower start found!!!" << std::endl;
            continue;
        }

        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "shower start", BLACK, 2);
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexPositionU, &projection, "direction", BLACK, 2, 1);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerSpineHitList, "SPINE", BLACK);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }

    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionU", "ShowerCounter", &showerCounterVectorU));
    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionU", "ShowerCoreProjectionL", &projectionVectorU));
    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionU", "Energy", &energiesU));
    //PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "EnergyDistributionU"));
    
    // Investigate V View
    std::cout << "Investigating V view..." << std::endl;

    DeviationAngleMap deviationAngleMapV;
    this->FillAngularDecompositionMap(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, deviationAngleMapV);

    IntVector angularPeakVectorV;
    this->ObtainViewPeakVector(deviationAngleMapV, angularPeakVectorV);

    IntVector showerCounterVectorV;
    FloatVector projectionVectorV, energiesV;

    IntVector investigatedPeaksV;

    for (unsigned int i = 0; i < angularPeakVectorV.size(); ++i)
    {
        int bestAngularPeak(0);
        if (!this->FindBestAngularPeak(deviationAngleMapV, angularPeakVectorV, investigatedPeaksV, bestAngularPeak))
            break;

        investigatedPeaksV.push_back(bestAngularPeak);

        const CartesianVector peakDirection(std::cos(bestAngularPeak * pAlgorithm->m_binSize), 0.f, std::sin(bestAngularPeak * pAlgorithm->m_binSize));
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);
        std::cout << "bestAngularPeak V view: " << bestAngularPeak << std::endl;
        std::cout << "peakDirection: " << peakDirection << std::endl;

        const CartesianVector nuVertexPositionV(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_V));
        const CartesianVector projection(nuVertexPositionV + (peakDirection * 14.f));

        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexPositionV, &projection, "direction", BLACK, 2, 1);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        CaloHitList showerSpineHitList;
        this->FindShowerSpine(pAlgorithm, pShowerPfo, nuVertexPositionV, peakDirection, TPC_VIEW_V, showerSpineHitList);

        if (showerSpineHitList.size() < 20)
            continue;

        //look at energy distribution
        EnergySpectrumMap energySpectrumMapV;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, showerCounterVectorV, projectionVectorV, energiesV, energySpectrumMapV);

        CaloHitList showerPfoHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_V, showerPfoHitList);

        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        if (!this->FindShowerStart(pAlgorithm, energySpectrumMapV, showerSpineHitList, showerStartPosition, showerPfoHitList, isEndDownstream))
        {
            std::cout << "no shower start found!!!" << std::endl;
            continue;
        }

        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "shower start", BLACK, 2);
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexPositionV, &projection, "direction", BLACK, 2, 1);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerSpineHitList, "SPINE", BLACK);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }

    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionV", "ShowerCounter", &showerCounterVectorV));
    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionV", "ShowerCoreProjectionL", &projectionVectorV));
    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionV", "Energy", &energiesV));
    //PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "EnergyDistributionV"));

    // Investigate W View
    std::cout << "Investigating W view..." << std::endl;

    DeviationAngleMap deviationAngleMapW;
    this->FillAngularDecompositionMap(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, deviationAngleMapW);

    IntVector angularPeakVectorW;
    this->ObtainViewPeakVector(deviationAngleMapW, angularPeakVectorW);

    IntVector showerCounterVectorW;
    FloatVector projectionVectorW, energiesW;

    IntVector investigatedPeaksW;

    for (unsigned int i = 0; i < angularPeakVectorW.size(); ++i)
    {
        int bestAngularPeak(0);
        if (!this->FindBestAngularPeak(deviationAngleMapW, angularPeakVectorW, investigatedPeaksW, bestAngularPeak))
            break;

        investigatedPeaksW.push_back(bestAngularPeak);

        const CartesianVector peakDirection(std::cos(bestAngularPeak * pAlgorithm->m_binSize), 0.f, std::sin(bestAngularPeak * pAlgorithm->m_binSize));
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);
        std::cout << "bestAngularPeak W view: " << bestAngularPeak << std::endl;
        std::cout << "peakDirection: " << peakDirection << std::endl;

        const CartesianVector nuVertexPositionW(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_W));
        const CartesianVector projection(nuVertexPositionW + (peakDirection * 14.f));

        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexPositionW, &projection, "direction", BLACK, 2, 1);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        CaloHitList showerSpineHitList;
        this->FindShowerSpine(pAlgorithm, pShowerPfo, nuVertexPositionW, peakDirection, TPC_VIEW_W, showerSpineHitList);

        if (showerSpineHitList.size() < 20)
            continue;

        //look at energy distribution
        EnergySpectrumMap energySpectrumMapW;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, showerCounterVectorW, projectionVectorW, energiesW, energySpectrumMapW);

        CaloHitList showerPfoHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_W, showerPfoHitList);

        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        if (!this->FindShowerStart(pAlgorithm, energySpectrumMapW, showerSpineHitList, showerStartPosition, showerPfoHitList, isEndDownstream))
        {
            std::cout << "no shower start found!!!" << std::endl;
            continue;
        }

        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "shower start", BLACK, 2);
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexPositionW, &projection, "direction", BLACK, 2, 1);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerSpineHitList, "SPINE", BLACK);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    }

    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionW", "ShowerCounter", &showerCounterVectorW));
    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionW", "ShowerCoreProjectionL", &projectionVectorW));
    //PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "EnergyDistributionW", "Energy", &energiesW));
    //PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "EnergyDistributionW"));
    
    //this->AngularDistributionTree(pAlgorithm, deviationAngleMapU, deviationAngleMapV, deviationAngleMapW);
    /////////////
        /*
    ProtoShowerVector protoShowerVector;
    this->BuildProtoShowers(pShowerPfo, protoShowerVector);

    std::cout << "DDDDDDDDDDD" << std::endl;

    if (protoShowerVector.empty())
        return false;

    std::cout << "EEEEEEEEEEEE" << std::endl;

    // pfo splitting alg? (can be quite simple if we are able to build the protoShowers sufficiently...)
    // set ProtoShower parent pfo address to match!
    if (protoShowerVector.size() != 1)
        return false;

    for (const ProtoShower &protoShower : protoShowerVector)
    {
        if (this->IsElectronPathway(protoShower))
            continue;

        this->RemoveConnectionPathway(protoShower);

        // change metadata to say that we think that this gamma is a gamma (so that we don't extend it in the future tool)
    }
        */


    /*
    this->FillAngularDecompositionMap(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, deviationAngleMapV);
    this->FillAngularDecompositionMap(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, deviationAngleMapW);

    AngularPeakVector angularPeakVector;
    this->ObtainAngularPeakVector(pAlgorithm, deviationAngleMapU, deviationAngleMapV, deviationAngleMapW, angularPeakVector);
    */
    ///////
    /*
    // do this for U view
    int highestBin(0);
    float highestWeight(-std::numeric_limits<float>::max());

    for (const auto &entry : pAlgorithm->m_thetaMapU)
    {
        float weight = entry.second.size();

        if (weight > highestWeight)
        {
            highestWeight = weight;
            highestBin = entry.first;
        }
    }

    std::cout << "highest peak: " << highestBin * pAlgorithm->m_binSize << std::endl;

    const CartesianVector peakDirection(std::cos(highestBin * pAlgorithm->m_binSize), 0.f, std::sin(highestBin * pAlgorithm->m_binSize));

    std::cout << "peakDirection: " << peakDirection << std::endl;

    CaloHitList showerSpineHitList;
    this->FindShowerSpine(pAlgorithm, pShowerPfo, nuVertexPositionU, peakDirection, TPC_VIEW_U, showerSpineHitList);

    CartesianVector projection(nuVertexPositionU + (peakDirection * 14.f));

    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexPositionU, &projection, "direction", BLACK, 2, 1);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerSpineHitList, "SPINE", BLACK);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    // ISOBEL - SHOULD I ALSO HAVE A HIT CUT?
    if (!this->HasPathToNuVertex(pShowerPfo, nuVertexPosition))
        return false;

    std::cout << "CCCCCCCCCCCC" << std::endl;

    ProtoShowerVector protoShowerVector;
    this->BuildProtoShowers(pShowerPfo, protoShowerVector);

    std::cout << "DDDDDDDDDDD" << std::endl;

    if (protoShowerVector.empty())
        return false;

    std::cout << "EEEEEEEEEEEE" << std::endl;

    // pfo splitting alg? (can be quite simple if we are able to build the protoShowers sufficiently...)
    // set ProtoShower parent pfo address to match!
    if (protoShowerVector.size() != 1)
        return false;

    for (const ProtoShower &protoShower : protoShowerVector)
    {
        if (this->IsElectronPathway(protoShower))
            continue;

        this->RemoveConnectionPathway(protoShower);

        // change metadata to say that we think that this gamma is a gamma (so that we don't extend it in the future tool)
    }
    */

    return true;
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
    FloatVector xCoordinate, yCoordinate, zCoordinate;
    FloatVector lCoordinateU, tCoordinateU;
    FloatVector lCoordinateV, tCoordinateV;
    FloatVector lCoordinateW, tCoordinateW;
    FloatVector hitEnergyU, hitEnergyV, hitEnergyW;

    float lNuU(-999.f), tNuU(-999.f);
    float lNuV(-999.f), tNuV(-999.f);
    float lNuW(-999.f), tNuW(-999.f);

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

    const TwoDSlidingFitResult twoDSlidingFitU(&initialHitPositionsU, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult twoDSlidingFitV(&initialHitPositionsV, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult twoDSlidingFitW(&initialHitPositionsW, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    twoDSlidingFitU.GetLocalPosition(nuVertexPosition, lNuU, tNuU);
    twoDSlidingFitV.GetLocalPosition(nuVertexPosition, lNuV, tNuV);
    twoDSlidingFitW.GetLocalPosition(nuVertexPosition, lNuW, tNuW);

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
        xCoordinate.push_back(hitPosition.GetX());
        yCoordinate.push_back(hitPosition.GetY());
        zCoordinate.push_back(hitPosition.GetZ());

        const TwoDSlidingFitResult &twoDSlidingFit(hitType == TPC_VIEW_U ? twoDSlidingFitU : hitType == TPC_VIEW_V ? twoDSlidingFitV : twoDSlidingFitW);
        FloatVector &lCoordinate(hitType == TPC_VIEW_U ? lCoordinateU : hitType == TPC_VIEW_V ? lCoordinateV : lCoordinateW);
        FloatVector &tCoordinate(hitType == TPC_VIEW_U ? tCoordinateU : hitType == TPC_VIEW_V ? tCoordinateV : tCoordinateW);

        float l, t;
        twoDSlidingFit.GetLocalPosition(hitPosition, l, t);

        lCoordinate.push_back(l);
        tCoordinate.push_back(t);

        FloatVector &deviationAngle(hitType == TPC_VIEW_U ? deviationAngleU : hitType == TPC_VIEW_V ? deviationAngleV : deviationAngleW);

        float deviationAngleTemp = xAxis.GetOpeningAngle(hitPosition - nuVertexPositionProjection);

        if (hitPosition.GetZ() < 0.f)
            deviationAngleTemp *= -1.f;

        deviationAngle.push_back(deviationAngleTemp);
    }

    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ShowerCounter", m_counter));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleU", &deviationAngleU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleV", &deviationAngleV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleW", &deviationAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "XCoordinate", &xCoordinate));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "YCoordinate", &yCoordinate));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ZCoordinate", &zCoordinate));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateU", &lCoordinateU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateV", &lCoordinateV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateW", &lCoordinateW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "TCoordinateU", &tCoordinateU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "TCoordinateV", &tCoordinateV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "TCoordinateW", &tCoordinateW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyU", &hitEnergyU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyV", &hitEnergyV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyW", &hitEnergyW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateNuU", lNuU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateNuV", lNuV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateNuW", lNuW));
    PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "ShowerDistribution"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::AngularDistributionTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const DeviationAngleMap &deviationAngleMapU, 
    const DeviationAngleMap &deviationAngleMapV, const DeviationAngleMap &deviationAngleMapW)
{
    IntVector deviationBinU, deviationBinV, deviationBinW;
    FloatVector deviationWeightU, deviationWeightV, deviationWeightW;

    for (const auto &entry : deviationAngleMapU)
    {
        deviationBinU.push_back(entry.first);
        deviationWeightU.push_back(entry.second);
    }

    for (const auto &entry : deviationAngleMapV)
    {
        deviationBinV.push_back(entry.first);
        deviationWeightV.push_back(entry.second);
    }

    for (const auto &entry : deviationAngleMapW)
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

void GammaStartRefinementTool::FillAngularDecompositionMap(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const CartesianVector &nuVertexPosition, const HitType hitType, DeviationAngleMap &deviationAngleMap)
{
    CaloHitList caloHits2D;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, caloHits2D);

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const CartesianVector nuVertexProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, hitType));

    for (const CaloHit *const pCaloHit : caloHits2D)
    {
        const HitType caloHitType(pCaloHit->GetHitType());

        if (hitType != caloHitType)
            std::cout << "ISOBEL - WRONG HIT TYPE!!" << std::endl;

        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if ((hitPosition - nuVertexProjection).GetMagnitude() > 14.f)
            continue;

        float deviationAngle = xAxis.GetOpeningAngle(hitPosition - nuVertexProjection);

        if ((hitPosition - nuVertexProjection).GetZ() < 0.f)
            deviationAngle *= -1.f;

        const int thetaFactor(std::floor(deviationAngle / pAlgorithm->m_binSize));

        if (deviationAngleMap.find(thetaFactor) == deviationAngleMap.end())
            deviationAngleMap[thetaFactor] = 1;
        else
            deviationAngleMap[thetaFactor] += 1;
    }

    this->SmoothAngularDecompositionMap(deviationAngleMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::SmoothAngularDecompositionMap(DeviationAngleMap &deviationAngleMap)
{
    for (auto &entry : deviationAngleMap)
        std::cout << "binNo: " << entry.first << ", calo hit list size: " << entry.second << std::endl;

    const int lowestThetaFactor(deviationAngleMap.begin()->first);
    const int highestThetaFactor(deviationAngleMap.rbegin()->first);

    std::cout << "lowestThetaFactor: " << lowestThetaFactor << std::endl;
    std::cout << "highestThetaFactor: " << highestThetaFactor << std::endl;

    const int smoothingWindow = 3;
    const int loopMin = (-1) * (smoothingWindow - 1) / 2;
    const int loopMax = (smoothingWindow - 1) / 2;

    DeviationAngleMap deviationAngleMapTemp(deviationAngleMap);

    deviationAngleMap.clear();

    for (auto &entry : deviationAngleMapTemp)
    {
        float total(0.f);
        int i(entry.first), count(0);

        for (int j = loopMin; j <= loopMax; ++j)
        {
            int bin = i + j;

            //if (bin < lowestThetaFactor)
                //total += 0;
                //bin = i + (smoothingWindow - j);

            //if (bin > highestThetaFactor)
            //continue;
                //bin = i - (smoothingWindow - j);

            if (deviationAngleMapTemp.find(bin) == deviationAngleMapTemp.end())
                total += 0;
            else
                total += deviationAngleMapTemp.at(bin);

            ++count;
        }

        deviationAngleMap[i] = total / count;
    }

    std::cout << "AFTER SMOOTHING" << std::endl;
    for (auto &entry : deviationAngleMap)
        std::cout << "binNo: " << entry.first << ", calo hit list size: " << entry.second << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::ObtainAngularPeakVector(ShowerStartRefinementAlgorithm *const pAlgorithm, DeviationAngleMap &deviationAngleMapU, 
    DeviationAngleMap &deviationAngleMapV, DeviationAngleMap &deviationAngleMapW, AngularPeakVector &angularPeakVector)
{
    IntVector peakVectorU, peakVectorV, peakVectorW;

    this->ObtainViewPeakVector(deviationAngleMapU, peakVectorU);
    this->ObtainViewPeakVector(deviationAngleMapV, peakVectorV);
    this->ObtainViewPeakVector(deviationAngleMapW, peakVectorW);

    for (unsigned int uPeak = 0; uPeak < peakVectorU.size(); ++uPeak)
    {
        const CartesianVector uDirection(std::cos(uPeak * pAlgorithm->m_binSize), 0.f, std::sin(uPeak * pAlgorithm->m_binSize));

        for (unsigned int vPeak = 0; vPeak < peakVectorV.size(); ++vPeak)
        {
            const CartesianVector vDirection(std::cos(vPeak * pAlgorithm->m_binSize), 0.f, std::sin(vPeak * pAlgorithm->m_binSize));

            for (unsigned int wPeak = 0; wPeak < peakVectorW.size(); ++wPeak)
            {
                const CartesianVector wDirection(std::cos(wPeak * pAlgorithm->m_binSize), 0.f, std::sin(wPeak * pAlgorithm->m_binSize));

                const CartesianVector uPrediction(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vDirection, wDirection));
                const CartesianVector vPrediction(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, wDirection, uDirection));
                const CartesianVector wPrediction(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, uDirection, vDirection));

                const float chiSquared(uDirection.GetOpeningAngle(uPrediction) + vDirection.GetOpeningAngle(vPrediction) + wDirection.GetOpeningAngle(wPrediction));

                float m_maxChiSquared = 10;

                if (chiSquared > m_maxChiSquared)
                    continue;

                const float binWeightSum(deviationAngleMapU.at(uPeak) + deviationAngleMapV.at(vPeak) + deviationAngleMapW.at(wPeak));

                angularPeakVector.push_back(AngularPeak(uPeak, vPeak, wPeak, chiSquared, binWeightSum));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::ObtainViewPeakVector(DeviationAngleMap &deviationAngleMap, IntVector &viewPeakVector)
{
    for (auto &entry : deviationAngleMap)
    {
        const float bin(entry.first);

        //if ((deviationAngleMap.find(bin) == deviationAngleMap.begin()) || (deviationAngleMap.find(bin) == std::prev(deviationAngleMap.end())))
        //continue;

        const float binWeight(entry.second);

        // Find lower bin 
        int precedingBin(bin - 1);
        bool foundPreceeding(false);

        while(!foundPreceeding)
        {
            //if (deviationAngleMap.find(precedingBin) == deviationAngleMap.begin())
            //break;

            if ((deviationAngleMap.find(precedingBin) == deviationAngleMap.end()) ||
                (std::fabs(deviationAngleMap.at(precedingBin) - deviationAngleMap.at(bin)) > std::numeric_limits<float>::epsilon()))
            {
                foundPreceeding = true;
                break;
            }

            --precedingBin;
        }

        if (!foundPreceeding)
            continue;

        int followingBin(bin + 1);
        bool foundFollowing(false);

        while(!foundFollowing)
        {
            // not sure i want this here
            //if (deviationAngleMap.find(followingBin) ==  std::prev(deviationAngleMap.end()))
            //break;

            if ((deviationAngleMap.find(followingBin) == deviationAngleMap.end()) ||
                (std::fabs(deviationAngleMap.at(followingBin) - deviationAngleMap.at(bin)) > std::numeric_limits<float>::epsilon()))
            {
                foundFollowing = true;
                break;
            }

            ++followingBin;
        }

        if (!foundFollowing)
            continue;        


        const float precedingBinWeight(deviationAngleMap.find(precedingBin) == deviationAngleMap.end() ? 0.f : deviationAngleMap.at(precedingBin));
        const float followingBinWeight(deviationAngleMap.find(followingBin) == deviationAngleMap.end() ? 0.f : deviationAngleMap.at(followingBin));

        if ((binWeight < precedingBinWeight) || (binWeight < followingBinWeight))
            continue;

        viewPeakVector.push_back(bin);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::FindBestAngularPeak(DeviationAngleMap &deviationAngleMap, IntVector &viewPeakVector, IntVector &investigatedPeaks, int &bestAngularPeak)
{
    bool found(false);
    float bestWeight(0.f);

    for (unsigned int i = 0; i < viewPeakVector.size(); ++i)
    {
        if (std::find(investigatedPeaks.begin(), investigatedPeaks.end(), viewPeakVector[i]) != investigatedPeaks.end())
            continue;

        if (deviationAngleMap.at(viewPeakVector[i]) > bestWeight)
        {
            found = true;
            bestWeight = deviationAngleMap.at(viewPeakVector[i]);
            bestAngularPeak = viewPeakVector[i];
        }
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::GetEnergyDistribution(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerSpineHitList, 
    IntVector &showerCounterVector, FloatVector &projectionVector, FloatVector &energies, EnergySpectrumMap &energySpectrumMap)
{
    ++m_showerCounter;

    float e0(0.f);

    for (const CaloHit *pCaloHit : showerSpineHitList)
    {
        const float energy(1000.f * pCaloHit->GetElectromagneticEnergy());

        e0 += energy;
    }

    LongitudinalPositionMap longitudinalPositionMap;
    this->ObtainLongitudinalDecomposition(pAlgorithm, showerSpineHitList, longitudinalPositionMap);

    for (const CaloHit *pCaloHit : showerSpineHitList)
    {
        const float energy(1000.f * pCaloHit->GetElectromagneticEnergy());
        const float fractionalEnergy(energy / e0);
        const float projection(longitudinalPositionMap.at(pCaloHit));

        // For filling tree
        showerCounterVector.push_back(m_showerCounter);
        projectionVector.push_back(projection);
        energies.push_back(fractionalEnergy);

        int longitudinalIndex = std::floor(projection / m_energySpectrumBinSize);

        if (energySpectrumMap.find(longitudinalIndex) == energySpectrumMap.end())
            energySpectrumMap[longitudinalIndex] = fractionalEnergy;
        else
            energySpectrumMap[longitudinalIndex] += fractionalEnergy;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::ObtainLongitudinalDecomposition(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerSpineHitList, 
    LongitudinalPositionMap &longitudinalPositionMap)
{
    CartesianPointVector hitPositions;

    for (const CaloHit *const pCaloHit : showerSpineHitList)
        hitPositions.push_back(pCaloHit->GetPositionVector());

    const TwoDSlidingFitResult twoDSlidingFit(&hitPositions, m_microSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    /*
    const LayerFitResultMap &layerFitResultMap(twoDSlidingFit.GetLayerFitResultMap());

    CartesianPointVector visualize;

    for (const auto &entry : layerFitResultMap)
    {
        const LayerFitResult &layerFitResult(entry.second);

        CartesianVector layerPosition(0.f, 0.f, 0.f);
        float layerL(layerFitResult.GetL()), layerT(layerFitResult.GetFitT());

        twoDSlidingFit.GetGlobalPosition(layerL, layerT, layerPosition);
        visualize.push_back(layerPosition);
    }

    for (const CartesianVector &pos : visualize)
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &pos, "layer position", GREEN, 2);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */

    const LayerFitResultMap &layerFitResultMap(twoDSlidingFit.GetLayerFitResultMap());

    std::map<int, CaloHitList> layerToHitMap;

    for (const CaloHit *const pCaloHit : showerSpineHitList)
    {
        float hitL(0.f), hitT(0.f);
        const CartesianVector hitPosition(pCaloHit->GetPositionVector());

        twoDSlidingFit.GetLocalPosition(hitPosition, hitL, hitT);
        layerToHitMap[twoDSlidingFit.GetLayer(hitL)].push_back(pCaloHit);

        //if (layerFitResultMap.size() < 2)
        //BLAHHH
    }

    for (const auto &entry : layerToHitMap)
        std::cout << "layer: " << entry.first << ", nHits: " << entry.second.size() << std::endl;

    float runningDistance(0);

    for (auto iter = layerToHitMap.begin(); iter != layerToHitMap.end(); ++iter)
    {
        const int layer(iter->first);
        const float layerL(layerFitResultMap.at(layer).GetL());

        const int higherLayer(std::next(iter) == layerToHitMap.end() ? layer : std::next(iter)->first);
        const int middleLayer(layer);
        const int lowerLayer(iter == layerToHitMap.begin() ? layer : std::prev(iter)->first);

        CartesianVector lowerLayerPosition(0.f, 0.f, 0.f), middleLayerPosition(0.f, 0.f, 0.f), higherLayerPosition(0.f, 0.f, 0.f);

        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(lowerLayer).GetL(), layerFitResultMap.at(lowerLayer).GetFitT(), lowerLayerPosition);
        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(middleLayer).GetL(), layerFitResultMap.at(middleLayer).GetFitT(), middleLayerPosition);
        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(higherLayer).GetL(), layerFitResultMap.at(higherLayer).GetFitT(), higherLayerPosition);

        const float layerLength = std::next(iter) == layerToHitMap.end() ? 0.f : iter == layerToHitMap.begin() ? 0.f : (middleLayerPosition - lowerLayerPosition).GetMagnitude();

        for (const CaloHit *const pCaloHit : iter->second)
        {
            const CartesianVector hitPosition(pCaloHit->GetPositionVector());

            float hitL(0.f), hitT(0.f);
            twoDSlidingFit.GetLocalPosition(hitPosition, hitL, hitT);

            //int localLowerLayer(0), localHigherLayer(0);
            CartesianVector localLowerLayerPosition(0.f, 0.f, 0.f), localHigherLayerPosition(0.f, 0.f, 0.f);

            if (hitL < layerL)
            {
                //localLowerLayer = lowerLayer;
                localLowerLayerPosition = lowerLayerPosition;
                //localHigherLayer = iter == layerToHitMap.begin() ? higherLayer : middleLayer;
                localHigherLayerPosition = iter == layerToHitMap.begin() ? higherLayerPosition : middleLayerPosition;
            }
            else
            {
                //localLowerLayer = std::next(iter) == layerToHitMap.end() ? lowerLayer : middleLayer;
                localLowerLayerPosition = std::next(iter) == layerToHitMap.end() ? lowerLayerPosition : middleLayerPosition;
                //localHigherLayer = higherLayer;
                localHigherLayerPosition = higherLayerPosition;
            }

            CartesianVector displacement(higherLayerPosition - lowerLayerPosition);
            displacement = displacement.GetUnitVector();

            float longitudinalDisplament = (hitL > layerL ? layerLength : 0.f);
            longitudinalDisplament += (displacement.GetDotProduct(hitPosition - lowerLayerPosition) + runningDistance);

            longitudinalPositionMap[pCaloHit] = longitudinalDisplament;
        }

        runningDistance += layerLength;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

    bool GammaStartRefinementTool::FindShowerStart(ShowerStartRefinementAlgorithm *const pAlgorithm, const EnergySpectrumMap &energySpectrumMap, const CaloHitList &showerSpineHitList, 
                                                   CartesianVector &showerStartPosition, const CaloHitList &showerPfoHitList, const bool isEndDownstream)
{
    int nBins(0);
    float mean(0);

    auto iter(energySpectrumMap.begin());

    for (int i = 0 ; i < m_nInitialEnergyBins; ++i)
    {
        mean += iter->second;
        ++nBins;
        ++iter;
    }

    std::cout << "nBins: " << nBins << std::endl;

    mean /= static_cast<float>(nBins);

    float sigma(0);

    iter = energySpectrumMap.begin();

    for (int i = 0; i < m_nInitialEnergyBins; ++i)
    {
        sigma += std::pow(iter->second - mean, 2);
        ++iter;
    }

    sigma = std::sqrt(sigma / nBins);

    std::cout << "Initial Sigma: " << sigma << std::endl;

    bool found(false);
    int longitudinalStartBin(-999);

    for (unsigned int i = m_nInitialEnergyBins; i < energySpectrumMap.size(); ++i)
    {
        const float energyDeviation(iter->second - mean);

        std::cout << "energyDeviation / sigma: " << (energyDeviation / sigma) << std::endl;

        if ((energyDeviation / sigma > m_minSigmaDeviation) && !this->IsTrackBlip(energySpectrumMap, iter, mean, sigma) && 
            this->IsShowerTopology(pAlgorithm, (iter->first * m_energySpectrumBinSize), showerPfoHitList, showerSpineHitList, isEndDownstream))
        {
            found = true;
            break;
        }

        ++iter;
    }

    if (!found)
        return false;

    CartesianPointVector hitPositions;

    for (const CaloHit *const pCaloHit : showerSpineHitList)
        hitPositions.push_back(pCaloHit->GetPositionVector());

    longitudinalStartBin = iter->first;
    float longitudinalStartCoordinate(longitudinalStartBin * m_energySpectrumBinSize);

    const TwoDSlidingFitResult twoDSlidingFit(&hitPositions, m_microSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    twoDSlidingFit.GetGlobalFitPosition(longitudinalStartCoordinate, showerStartPosition);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::IsTrackBlip(const EnergySpectrumMap &energySpectrumMap, const EnergySpectrumMap::const_iterator &iter, const float initialMean, const float initialSigma)
{
    if (iter == std::prev(energySpectrumMap.end()))
        return true;

    //std::cout << "investigating whether bin: " << iter->first << " is a track blip" << std::endl;

    auto copyIter(std::next(iter));

    int nBins(0);
    float mean(0.f);

    for (int i = 0 ; i < m_trackSearchWindow; ++i)
    {
        mean += std::fabs(copyIter->second - initialMean) / initialSigma;
        ++nBins;
        ++copyIter;

        if (copyIter == energySpectrumMap.end())
            break;
    }

    mean /= static_cast<float>(nBins);


    return (mean < m_minTrackBlipMean);

    /*
    float sigma(0);

    copyIter = std::next(iter);

    for (int i = 0; i < m_trackSearchWindow; ++i)
    {
        sigma += std::pow(((copyIter->second - initialMean) / initialSigma) - mean, 2);
        ++copyIter;

        if (copyIter == energySpectrumMap.end())
            break;
    }

    sigma = std::sqrt(sigma / nBins);

    //std::cout << "sigma of post hits: " << sigma << std::endl;

    if (sigma > 2.f)
        return false;

    return true; 
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

    bool GammaStartRefinementTool::IsShowerTopology(ShowerStartRefinementAlgorithm *const pAlgorithm, const float longitudinalDistance, const CaloHitList &showerPfoHits, const CaloHitList &showerSpineHits, const bool isEndDownstream)
{
    CartesianPointVector hitPositions;

    for (const CaloHit *const pCaloHit : showerSpineHits)
        hitPositions.push_back(pCaloHit->GetPositionVector());

    const TwoDSlidingFitResult twoDSlidingFit(&hitPositions, m_microSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const LayerFitResultMap &layerFitResultMap(twoDSlidingFit.GetLayerFitResultMap());

    float runningDistance(0.f);

    CartesianVector previousLayerPosition(0.f, 0.f, 0.f);
    twoDSlidingFit.GetGlobalPosition(layerFitResultMap.begin()->second.GetL(), layerFitResultMap.begin()->second.GetFitT(), previousLayerPosition);

    int showerStartLayer(0);

    for (auto iter = std::next(layerFitResultMap.begin()); iter != layerFitResultMap.end(); ++iter)
    {
        const int layer(iter->first);

        CartesianVector layerPosition(0.f, 0.f, 0.f);
        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(layer).GetL(), layerFitResultMap.at(layer).GetFitT(), layerPosition);

        runningDistance += (layerPosition - previousLayerPosition).GetMagnitude();

        if (runningDistance > longitudinalDistance)
        {
            showerStartLayer = layer;
            break;
        }

        previousLayerPosition = layerPosition;
    }

    const float lCoordinate(layerFitResultMap.at(showerStartLayer).GetL()), tCoordinate(layerFitResultMap.at(showerStartLayer).GetFitT());
    const float localGradient(layerFitResultMap.at(showerStartLayer).GetGradient());

    CartesianVector showerStartPosition(0.f, 0.f, 0.f), showerStartDirection(0.f, 0.f, 0.f);

    twoDSlidingFit.GetGlobalPosition(lCoordinate, tCoordinate, showerStartPosition);
    twoDSlidingFit.GetGlobalDirection(localGradient, showerStartDirection);

    ///////////////////////////////////////////////////////
    CartesianVector positiveEdgeStart(0.f, 0.f, 0.f), positiveEdgeEnd(0.f, 0.f, 0.f), positiveEdgeDirection(0.f, 0.f, 0.f);
    CartesianVector negativeEdgeStart(0.f, 0.f, 0.f), negativeEdgeEnd(0.f, 0.f, 0.f), negativeEdgeDirection(0.f, 0.f, 0.f);

    bool isBetween(false);
    bool doesStraddle(false);
    // make this a loop? and if we cant find a place where the end is larger than abort. and demand shower length? in terms of fraction of shower length?
    if (this->CharacteriseShower(pAlgorithm, showerPfoHits, showerSpineHits, showerStartPosition, showerStartDirection, isEndDownstream, 1.0,
                                 positiveEdgeStart, positiveEdgeEnd, positiveEdgeDirection, negativeEdgeStart, negativeEdgeEnd, negativeEdgeDirection, isBetween, doesStraddle) != STATUS_CODE_SUCCESS)
    {
        std::cout << "failed to fit the shower with full length" << std::endl;
        return false;
    }

    if (!isBetween)
    {
        std::cout << "shower start is not inbetween start edge points :( " << std::endl;
        return false;
    }

    if (!doesStraddle)
    {
        std::cout << "shower does not straddle core " << std::endl;
        return false;
    }

    positiveEdgeStart = showerStartPosition;
    negativeEdgeStart = showerStartPosition;
    positiveEdgeDirection = positiveEdgeEnd - positiveEdgeStart;
    negativeEdgeDirection = negativeEdgeEnd - negativeEdgeStart;

    const float showerOpeningAngle(positiveEdgeDirection.GetOpeningAngle(negativeEdgeDirection) * 180.f / 3.14);
    const float positiveEdgeDeviation(positiveEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14);
    const float negativeEdgeDeviation(negativeEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14);

    std::cout << "edge opening angle: " << (positiveEdgeDirection.GetOpeningAngle(negativeEdgeDirection) * 180.f / 3.14) << std::endl;
    std::cout << "positive opening angle: " << (positiveEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14) << std::endl;
    std::cout << "negative opening angle: " << (negativeEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14) << std::endl;

    if (showerOpeningAngle < 5.f)
    {
        std::cout << "shower opening angle is too small" << std::endl;
        return false;
    }

    if (positiveEdgeDeviation > 45.f)
    {
        std::cout << "positive edge is not in the direction of the shower" << std::endl;
    }

    if (negativeEdgeDeviation > 45.f)
    {
        std::cout << "negative edge is not in the direction of the shower" << std::endl;
    }

    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &positiveEdgeStart, &positiveEdgeEnd, "positive direction", BLACK, 2, 1);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &negativeEdgeStart, &negativeEdgeEnd, "negative direction", BLACK, 2, 1);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    StatusCode GammaStartRefinementTool::CharacteriseShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerPfoHits, const CaloHitList &showerSpineHits, const CartesianVector &showerStartPosition, const CartesianVector &showerStartDirection, const bool isEndDownstream, float showerLengthFraction, CartesianVector &positiveEdgeStart, CartesianVector &positiveEdgeEnd, CartesianVector &/*positiveEdgeDirection*/, CartesianVector &negativeEdgeStart, CartesianVector &negativeEdgeEnd, CartesianVector &/*negativeEdgeDirection*/, bool &isBetween, bool &doesStraddle)
{
    // Find minL and maxL
    float minL(999), maxL(0.f);
    CaloHitList coreHitList;

    for (const CaloHit *const pCaloHit : showerPfoHits)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const float l(showerStartDirection.GetDotProduct(hitPosition - showerStartPosition));

        if ((isEndDownstream && (l < 0.f)) || (!isEndDownstream && (l > 0.f)))
            continue;

        if (std::find(showerSpineHits.begin(), showerSpineHits.end(), pCaloHit) == showerSpineHits.end())
            continue;

        if (l < minL)
            minL = l;

        if (l > maxL)
            maxL = l;

        coreHitList.push_back(pCaloHit);
    }

    // Get halo hits
    CartesianPointVector haloHits;
    CaloHitList haloHitList;

    for (const CaloHit *const pCaloHit : showerPfoHits)
    {
        if (std::find(coreHitList.begin(), coreHitList.end(), pCaloHit) != coreHitList.end())
            continue;

        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const float l(showerStartDirection.GetDotProduct(hitPosition - showerStartPosition));

        if ((isEndDownstream && (l < minL)) || (!isEndDownstream && (l > 0.f)))
            continue;

        if (l > (maxL * showerLengthFraction)) //can just do this if the lower end is less wide than the start end 
            continue;

        if (showerStartDirection.GetCrossProduct(hitPosition - showerStartPosition).GetMagnitude() > m_molliereRadius)
            continue;

        haloHits.push_back(pCaloHit->GetPositionVector());
        haloHitList.push_back(pCaloHit);
    }

    CartesianPointVector coordinateListP, coordinateListN;

    // Attempt to find end layer of shower...
    try 
    {
        // Perform shower sliding fit
        const TwoDSlidingShowerFitResult twoDShowerSlidingFit(&haloHits, m_showerSlidingFitWindow,  LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(twoDShowerSlidingFit.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(twoDShowerSlidingFit.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(twoDShowerSlidingFit.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        //////////////////////////////////////
        // Find edge positions for visualisation
        for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
        {
            LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
            LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

            CartesianVector positiveEdgePosition(0.f, 0.f, 0.f), negativeEdgePosition(0.f, 0.f, 0.f);
            if (layerFitResultMapP.end() != iterP)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterP->second.GetL(), iterP->second.GetFitT(), positiveEdgePosition);
                coordinateListP.push_back(positiveEdgePosition);
            }

            if (layerFitResultMapN.end() != iterN)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterN->second.GetL(), iterN->second.GetFitT(), negativeEdgePosition);
                coordinateListN.push_back(negativeEdgePosition);
            }
        }

        std::sort(coordinateListP.begin(), coordinateListP.end(), SortByDistanceToPoint(showerStartPosition));
        std::sort(coordinateListN.begin(), coordinateListN.end(), SortByDistanceToPoint(showerStartPosition));

        int countP(0);
        float pMaximumL(std::numeric_limits<float>::max());
        CartesianVector previousCoordinateP(0.f, 0.f, 0.f);

        for (const CartesianVector &coordinateP : coordinateListP)
        {
            if (countP != 0)
            {
                const float separation((coordinateP - previousCoordinateP).GetMagnitude());

                float thisT(0.f);
                twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(coordinateP, pMaximumL, thisT);

                if (separation > 5.f)
                {
                    std::cout << "gap in positive hits" << std::endl;
                    break;
                }
            }

            ++countP;
            previousCoordinateP = coordinateP;
        }

        int countN(0);
        float nMaximumL(std::numeric_limits<float>::max());
        CartesianVector previousCoordinateN(0.f, 0.f, 0.f);

        for (const CartesianVector &coordinateN : coordinateListN)
        {
            if (countN != 0)
            {
                const float separation((coordinateN - previousCoordinateN).GetMagnitude());

                float thisT(0.f);
                twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(coordinateN, nMaximumL, thisT);

                if (separation > 5.f)
                {
                    std::cout << "gap in negative hits" << std::endl;
                    break;
                }
            }

            ++countN;
            previousCoordinateN = coordinateN;
        }

        // now refil the core/halo hits...
        CaloHitList coreHitListTemp(coreHitList);
        CaloHitList haloHitListTemp(haloHitList);

        coreHitList.clear(); haloHits.clear(); coordinateListP.clear(); coordinateListN.clear(); 

        for (const CaloHit *const pCaloHit : coreHitListTemp)
        {
            float thisL(0.f), thisT(0.f);
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(hitPosition, thisL, thisT);

            if (thisL > std::max(pMaximumL, nMaximumL))
                continue;

            coreHitList.push_back(pCaloHit);
        }

        for (const CaloHit *const pCaloHit : haloHitListTemp)
        {
            float thisL(0.f), thisT(0.f);
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(hitPosition, thisL, thisT);

            if (thisL > std::max(pMaximumL, nMaximumL))
                continue;

            haloHits.push_back(pCaloHit->GetPositionVector());
        }
    }
    catch (const StatusCodeException &)
    {
        return STATUS_CODE_FAILURE;
    }

    //////////////////////////////////////
    CartesianVector end(showerStartPosition + (showerStartDirection * 10.f));
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, &end, "direction", VIOLET, 2, 1);

    for (auto &entry : coreHitList)
    {
        CartesianVector hitPosition(entry->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "core", BLUE, 2);
    }

    for (auto &entry : haloHits)
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &entry, "halo", RED, 2);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    //////////////////////////////////////

    // now for reals...
    try
    {
        // Perform shower sliding fit
        const TwoDSlidingShowerFitResult twoDShowerSlidingFit(&haloHits, m_showerSlidingFitWindow,  LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(twoDShowerSlidingFit.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(twoDShowerSlidingFit.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(twoDShowerSlidingFit.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        //////////////////////////////////////
        // Find edge positions for visualisation
        CartesianPointVector coordinateListP, coordinateListN;

        int layerCount(0);
        bool isFirstBetween(false), isLastBetween(false);

        for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
        {
            LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
            LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

            CartesianVector positiveEdgePosition(0.f, 0.f, 0.f), negativeEdgePosition(0.f, 0.f, 0.f);
            if (layerFitResultMapP.end() != iterP)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterP->second.GetL(), iterP->second.GetFitT(), positiveEdgePosition);
                coordinateListP.push_back(positiveEdgePosition);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &positiveEdgePosition, "positive edge", GREEN, 2);
            }

            if (layerFitResultMapN.end() != iterN)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterN->second.GetL(), iterN->second.GetFitT(), negativeEdgePosition);
                coordinateListN.push_back(negativeEdgePosition);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &negativeEdgePosition, "negative edge", RED, 2);
            }

            if ((layerFitResultMapP.end() != iterP) && (layerFitResultMapN.end() != iterN))
            {
                ++layerCount;

                const CartesianVector positiveDisplacement((positiveEdgePosition - showerStartPosition).GetUnitVector());
                const float positiveOpeningAngleFromCore(showerStartDirection.GetOpeningAngle(positiveDisplacement));

                if ((positiveOpeningAngleFromCore * 180 / 3.14) > 45.f)
                    continue;

                const CartesianVector positiveClockwiseRotation(positiveDisplacement.GetZ() * std::sin(positiveOpeningAngleFromCore) + 
                    positiveDisplacement.GetX() * std::cos(positiveOpeningAngleFromCore), 0.f, positiveDisplacement.GetZ() * std::cos(positiveOpeningAngleFromCore) - 
                    positiveDisplacement.GetX() * std::sin(positiveOpeningAngleFromCore));
                const CartesianVector positiveAnticlockwiseRotation(positiveDisplacement.GetZ() * std::sin(-1.f * positiveOpeningAngleFromCore) + 
                    positiveDisplacement.GetX() * std::cos(-1.f * positiveOpeningAngleFromCore), 0.f, positiveDisplacement.GetZ() * std::cos(-1.f * positiveOpeningAngleFromCore) - 
                    positiveDisplacement.GetX() * std::sin(-1.f * positiveOpeningAngleFromCore));

                const float positiveClockwiseT((positiveClockwiseRotation - showerStartDirection).GetMagnitude());
                const float positiveAnticlockwiseT((positiveAnticlockwiseRotation - showerStartDirection).GetMagnitude());
                const bool isPositiveClockwise(positiveClockwiseT < positiveAnticlockwiseT);

                const CartesianVector negativeDisplacement((negativeEdgePosition - showerStartPosition).GetUnitVector());
                const float negativeOpeningAngleFromCore(showerStartDirection.GetOpeningAngle(negativeDisplacement));

                if ((negativeOpeningAngleFromCore * 180 / 3.14) > 45.f)
                    continue;

                const CartesianVector negativeClockwiseRotation(negativeDisplacement.GetZ() * std::sin(negativeOpeningAngleFromCore) + 
                    negativeDisplacement.GetX() * std::cos(negativeOpeningAngleFromCore), 0.f, negativeDisplacement.GetZ() * std::cos(negativeOpeningAngleFromCore) - 
                    negativeDisplacement.GetX() * std::sin(negativeOpeningAngleFromCore));
                const CartesianVector negativeAnticlockwiseRotation(negativeDisplacement.GetZ() * std::sin(-1.f * negativeOpeningAngleFromCore) + 
                    negativeDisplacement.GetX() * std::cos(-1.f * negativeOpeningAngleFromCore), 0.f, negativeDisplacement.GetZ() * std::cos(-1.f * negativeOpeningAngleFromCore) - 
                    negativeDisplacement.GetX() * std::sin(-1.f * negativeOpeningAngleFromCore));

                const float negativeClockwiseT((negativeClockwiseRotation - showerStartDirection).GetMagnitude());
                const float negativeAnticlockwiseT((negativeAnticlockwiseRotation - showerStartDirection).GetMagnitude());
                const bool isNegativeClockwise(negativeClockwiseT < negativeAnticlockwiseT);

                if (!doesStraddle)
                    doesStraddle = (isPositiveClockwise != isNegativeClockwise);

                if (layerCount == 1)
                    isFirstBetween = doesStraddle;

                isLastBetween = (isPositiveClockwise != isNegativeClockwise);
                /*
                if (doesStraddle)
                {
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerStartDirection, "showerStartDirection", VIOLET, 2);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positiveClockwiseRotation, "positiveClockwiseRotation", BLUE, 2);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positiveAnticlockwiseRotation, "positiveAnticlockwiseRotation", BLUE, 2);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &negativeClockwiseRotation, "negativeClockwiseRotation", BLUE, 2);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &negativeAnticlockwiseRotation, "negativeAnticlockwiseRotation", BLUE, 2);
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                }
                */
            }
        }

        isBetween = (isFirstBetween || isLastBetween);
        //////////////////////////////////////



        // Find extremal coordinates

        float showerStartL(0.f), showerStartT(0.f);
        twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(showerStartPosition, showerStartL, showerStartT);

        const int startLayer(twoDShowerSlidingFit.GetShowerFitResult().GetLayer(showerStartL));
        const int endLayer(twoDShowerSlidingFit.GetShowerFitResult().GetMaxLayer());

        // positive
        int positiveStartLayer(0);
        int positiveEndLayer(0);

        for (auto &entry : layerFitResultMapN)
        {
            const int bestStartSeparation(std::abs(startLayer - positiveStartLayer));
            const int thisStartSeparation(std::abs(startLayer - entry.first));

            if ((std::abs(bestStartSeparation - thisStartSeparation) > 0) && (thisStartSeparation < bestStartSeparation))
                positiveStartLayer = entry.first;

            const int bestEndSeparation(std::abs(endLayer - positiveEndLayer));
            const int thisEndSeparation(std::abs(endLayer - entry.first));

            if ((std::abs(bestEndSeparation - thisEndSeparation) > 0) && (thisEndSeparation < bestEndSeparation))
                positiveEndLayer = entry.first;
        }

        if (std::abs(startLayer - positiveStartLayer) > 2)
        {
            std::cout << "positiveStartLayer too far away" << std::endl;
            return STATUS_CODE_FAILURE;
        }

        // negative
        int negativeStartLayer(0);
        int negativeEndLayer(0);

        for (auto &entry : layerFitResultMapN)
        {
            const int bestStartSeparation(std::abs(startLayer - negativeStartLayer));
            const int thisStartSeparation(std::abs(startLayer - entry.first));

            if ((std::abs(bestStartSeparation - thisStartSeparation) > 0) && (thisStartSeparation < bestStartSeparation))
                negativeStartLayer = entry.first;

            const int bestEndSeparation(std::abs(endLayer - negativeEndLayer));
            const int thisEndSeparation(std::abs(endLayer - entry.first));

            if ((std::abs(bestEndSeparation - thisEndSeparation) > 0) && (thisEndSeparation < bestEndSeparation))
                negativeEndLayer = entry.first;
        }

        if (std::abs(startLayer - negativeStartLayer) > 2)
        {
            std::cout << "negativeStartLayer too far away" << std::endl;
            return STATUS_CODE_FAILURE;
        }

        // should this be with respect to the shower direction??

        const float showerStartPositiveLocalL(layerFitResultMapP.at(positiveStartLayer).GetL()), showerStartPositiveLocalT(layerFitResultMapP.at(positiveStartLayer).GetFitT());
        const float showerEndPositiveLocalL(layerFitResultMapP.at(positiveEndLayer).GetL()), showerEndPositiveLocalT(layerFitResultMapP.at(positiveEndLayer).GetFitT());
        const float showerStartNegativeLocalL(layerFitResultMapN.at(negativeStartLayer).GetL()), showerStartNegativeLocalT(layerFitResultMapN.at(negativeStartLayer).GetFitT());
        const float showerEndNegativeLocalL(layerFitResultMapN.at(negativeEndLayer).GetL()), showerEndNegativeLocalT(layerFitResultMapN.at(negativeEndLayer).GetFitT());

        twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(showerStartPositiveLocalL, showerStartPositiveLocalT, positiveEdgeStart);
        twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(showerEndPositiveLocalL, showerEndPositiveLocalT, positiveEdgeEnd);
        twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(showerStartNegativeLocalL, showerStartNegativeLocalT, negativeEdgeStart);
        twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(showerEndNegativeLocalL, showerEndNegativeLocalT, negativeEdgeEnd);
        /*
        std::cout << "startLayer: " << startLayer << std::endl;
        std::cout << "positiveStartLayer: " << positiveStartLayer << std::endl;
        std::cout << "Start (L, T): (" << showerStartPositiveLocalL << ", " << showerStartPositiveLocalT << ")" << std::endl;
        std::cout << "End (L, T): (" << showerEndPositiveLocalL << ", " << showerEndPositiveLocalT << ")" << std::endl;
        std::cout << "negativeStartLayer: " << negativeStartLayer << std::endl;
        std::cout << "Start (L, T): (" << showerStartNegativeLocalL << ", " << showerStartNegativeLocalT << ")" << std::endl;
        std::cout << "End (L, T): (" << showerEndNegativeLocalL << ", " << showerEndNegativeLocalT << ")" << std::endl;

        const float showerStartPositiveT(showerStartDirection.GetCrossProduct(positiveEdgeStart - showerStartPosition).GetMagnitude());
        const float showerEndPositiveT(showerStartDirection.GetCrossProduct(positiveEdgeEnd - showerStartPosition).GetMagnitude());
        const float showerStartNegativeT(showerStartDirection.GetCrossProduct(negativeEdgeStart - showerStartPosition).GetMagnitude());
        const float showerEndNegativeT(showerStartDirection.GetCrossProduct(negativeEdgeEnd - showerStartPosition).GetMagnitude());


        // work out which side we are on?? 
        const float positiveStartDeviationAngle(showerStartDirection.GetOpeningAngle(positiveEdgeStart));
        const float negativeStartDeviationAngle(showerStartDirection.GetOpeningAngle(negativeEdgeStart));
        const float showerStartDeviationAngle(showerStartDirection.GetOpeningAngle(showerStartPosition));
        const float showerStartToPositive(showerStartPosition.GetOpeningAngle(positiveEdgeStart));
        const float showerStartToNegative(showerStartPosition.GetOpeningAngle(negativeEdgeStart));

        const float positiveStartSum(std::fabs(positiveStartDeviationAngle - (showerStartDeviationAngle + showerStartToPositive)));
        const float negativeStartSum(std::fabs(negativeStartDeviationAngle - (showerStartDeviationAngle + showerStartToNegative)));
        const bool isStartPositive(positiveStartSum < negativeStartSum);
        const float startT(isStartPositive ? showerStartPositiveT : showerStartNegativeT);

        std::cout << "positiveStartSum: " << positiveStartSum << std::endl;
        std::cout << "negativeStartSum: " << negativeStartSum << std::endl;

        if ((positiveStartSum > ()) && (negativeStartSum > std::numeric_limits<float>::epsilon()))
        {
            isBetween = false;
        }
        else if (showerStartT > (startT + 1.f))
        {
            const float positiveEndDeviationAngle(showerStartDirection.GetOpeningAngle(positiveEdgeEnd));
            const float negativeEndDeviationAngle(showerStartDirection.GetOpeningAngle(negativeEdgeEnd));
            const float showerEndDeviationAngle(showerStartDirection.GetOpeningAngle(showerStartPosition));
            const float showerEndToPositive(showerStartPosition.GetOpeningAngle(positiveEdgeEnd));
            const float showerEndToNegative(showerStartPosition.GetOpeningAngle(negativeEdgeEnd));

            const float positiveEndSum(std::fabs(positiveEndDeviationAngle - (showerEndDeviationAngle + showerEndToPositive)));
            const float negativeEndSum(std::fabs(negativeEndDeviationAngle - (showerEndDeviationAngle + showerEndToNegative)));
            const bool isEndPositive(positiveEndSum < negativeEndSum);
            const float endT(isEndPositive ? showerEndPositiveT : showerEndNegativeT);

            std::cout << "positiveEndSum: " << positiveEndSum << std::endl;
            std::cout << "negativeEndSum: " << negativeEndSum << std::endl;

            if ((negativeEndSum > std::numeric_limits<float>::epsilon()) && (negativeEndSum > std::numeric_limits<float>::epsilon()))
                isBetween = false;

            if (showerStartT > endT)
                isBetween = false;
        }
        else
        {
            isBetween = true;
        }
        */

        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &positiveEdgeStart, "Positive Edge", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &negativeEdgeStart, "Negative Edge", BLACK, 2);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        // Obtain direction
        //const float localPositiveGradient(layerFitResultMapP.at(positiveStartLayer).GetGradient());
        //const float localNegativeGradient(layerFitResultMapN.at(negativeStartLayer).GetGradient());
        //twoDShowerSlidingFit.GetNegativeEdgeFitResult().GetGlobalDirection(localNegativeGradient, negativeEdgeDirection);
        //twoDShowerSlidingFit.GetPositiveEdgeFitResult().GetGlobalDirection(localPositiveGradient, positiveEdgeDirection);

        // Swap start and end if needs be
        if ((showerStartPosition - positiveEdgeEnd).GetMagnitudeSquared() < (showerStartPosition - positiveEdgeStart).GetMagnitude())
        {
            const CartesianVector startTemp(positiveEdgeStart), endTemp(positiveEdgeEnd);

            positiveEdgeStart = endTemp;
            positiveEdgeEnd = startTemp;
        }

        if ((showerStartPosition - negativeEdgeEnd).GetMagnitudeSquared() < (showerStartPosition - negativeEdgeStart).GetMagnitude())
        {
            const CartesianVector startTemp(negativeEdgeStart), endTemp(negativeEdgeEnd);

            negativeEdgeStart = endTemp;
            negativeEdgeEnd = startTemp;
        }
    }
    catch (const StatusCodeException &)
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::RemoveConnectionPathway(const ProtoShower &protoShower)
{
    std::cout << protoShower.m_showerCore.m_startPosition.GetX() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MicroSlidingFitWindow", m_microSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinSigmaDeviation", m_minSigmaDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrackSearchWindow", m_trackSearchWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "NInitialEnergyBins", m_nInitialEnergyBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinTrackBlipMean", m_minTrackBlipMean));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShowerSlidingFitWindow", m_showerSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MolliereRadius", m_molliereRadius));

    ShowerStartRefinementBaseTool::ReadSettings(xmlHandle);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
