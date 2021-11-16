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

#include "larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h"

using namespace pandora;

namespace lar_content
{

GammaStartRefinementTool::GammaStartRefinementTool() : m_counter(0), m_showerCounter(0)
{
}

GammaStartRefinementTool::~GammaStartRefinementTool()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ShowerDistribution", "ShowerDistribution.root", "UPDATE"));
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

    DeviationAngleMap deviationAngleMapU;
    this->FillAngularDecompositionMap(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, deviationAngleMapU);

    IntVector angularPeakVectorU;
    this->ObtainViewPeakVector(deviationAngleMapU, angularPeakVectorU);

    bool changes(true);

    IntVector investigatedPeaks;

    while(changes)
    {
        changes = false;

        const int bestAngularPeak(this->FindBestAngularPeak(deviationAngleMapU, angularPeakVectorU));

        const CartesianVector peakDirection(std::cos(bestAngularPeak * pAlgorithm->m_binSize), 0.f, std::sin(bestAngularPeak * pAlgorithm->m_binSize));

        std::cout << "peakDirection: " << peakDirection << std::endl;

        const CartesianVector nuVertexPositionU(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_U));

        CaloHitList showerSpineHitList;
        LongitudinalPositionMap longitudinalPositionMap;
        this->FindShowerSpine(pAlgorithm, pShowerPfo, nuVertexPositionU, peakDirection, TPC_VIEW_U, showerSpineHitList, longitudinalPositionMap);

        CartesianVector projection(nuVertexPositionU + (peakDirection * 14.f));

        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexPositionU, &projection, "direction", BLACK, 2, 1);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerSpineHitList, "SPINE", BLACK);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        //look at energy distribution
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, longitudinalPositionMap);

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
    }




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

        if (hitPosition.GetZ() < 0.f)
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

    const int smoothingWindow = 5;
    const int loopMin = (1) * (smoothingWindow - 1) / 2;
    const int loopMax = 1 + (smoothingWindow - 1) / 2;

    DeviationAngleMap deviationAngleMapTemp(deviationAngleMap);

    deviationAngleMap.clear();

    for (auto &entry : deviationAngleMapTemp)
    {
        float total(0.f);
        int i(entry.first);

        for (int j = loopMin; j <= loopMax; ++j)
        {
            int bin = i + j;

            if (bin < lowestThetaFactor)
                bin = i + (smoothingWindow - j);

            if (bin > highestThetaFactor)
                bin = i - (smoothingWindow - j);

            if (deviationAngleMapTemp.find(bin) == deviationAngleMapTemp.end())
                total += 0;
            else
                total += deviationAngleMapTemp.at(bin);
        }

        deviationAngleMap[i] = total / smoothingWindow;
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

        if ((deviationAngleMap.find(bin) == deviationAngleMap.begin()) || (deviationAngleMap.find(bin) == std::prev(deviationAngleMap.end())))
            continue;

        const float binWeight(entry.second);

        /*
        if (binWeight < 2.f)
            continue;
        */
        const int precedingBin(bin - 1), followingBin(bin + 1);
        const float precedingBinWeight(deviationAngleMap.find(precedingBin) == deviationAngleMap.end() ? 0.f : deviationAngleMap.at(precedingBin));
        const float followingBinWeight(deviationAngleMap.find(followingBin) == deviationAngleMap.end() ? 0.f : deviationAngleMap.at(followingBin));

        if ((binWeight < precedingBinWeight) || (binWeight < followingBinWeight))
            continue;

        viewPeakVector.push_back(bin);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int GammaStartRefinementTool::FindBestAngularPeak(DeviationAngleMap &deviationAngleMap, IntVector &viewPeakVector)
{
    //ISOBEL - HAVE A CHECK OF 'IF NOT FOUND?'

    float bestWeight(0.f);
    int bestAngularPeak(0);

    for (unsigned int i = 0; i < viewPeakVector.size(); ++i)
    {
        // ISOBEL - HAVE A CHECK ON IF NOT FOUND

        if (deviationAngleMap.at(viewPeakVector[i]) > bestWeight)
        {
            bestWeight = deviationAngleMap.at(viewPeakVector[i]);
            bestAngularPeak = viewPeakVector[i];
        }
    }

    return bestAngularPeak;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::GetEnergyDistribution(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerSpineHitList, const LongitudinalPositionMap &longitudinalPositionMap)
{
    ++m_showerCounter;

    FloatVector projectionVector, energies;

    for (const CaloHit *pCaloHit : showerSpineHitList)
    {
        const float projection(longitudinalPositionMap.at(pCaloHit));

        projectionVector.emplace_back(projection);
        energies.emplace_back(1000.f * pCaloHit->GetElectromagneticEnergy());
    }

    float e0(0.f);

    for (const float energy : energies)
        e0 += energy;

    for (unsigned int i = 0; i < energies.size(); ++i)
        energies[i] = energies[i] / e0;

    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ShowerCounter", m_showerCounter));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ShowerCoreProjectionL", &projectionVector));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "Energy", &energies));
    PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "ShowerDistribution"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::RemoveConnectionPathway(const ProtoShower &protoShower)
{
    std::cout << protoShower.m_showerCore.m_startPosition.GetX() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
