/**
 *  @file   larpandoracontent/LArControlFlow/TrackDirectionTool.cc
 *
 *  @brief  Implementation of the track direction finding Tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArControlFlow/TrackDirectionTool.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

//----------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

TrackDirectionTool::TrackDirectionTool() :
    m_slidingFitWindow(5),
    m_minClusterCaloHits(5),
    m_minClusterLength(1.f),
    m_numberTrackEndHits(100000),
    m_enableFragmentRemoval(true),
    m_enableSplitting(true),
    m_tableInitialEnergy(2000.f),
    m_tableStepSize(0.5f),
    m_writeTable(false),
    m_lookupTableFileName("lookuptable.root"),
    m_treeName("lookuptable")
{
}

//--------------------------------------------------------------------
StatusCode TrackDirectionTool::Initialize()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
void TrackDirectionTool::FindDirections(const pandora::ParticleFlowObject *const pPfo, float &downProbability, const MasterAlgorithm *const)
{
    try
    {
        TrackDirectionTool::DirectionFitObject fitResult = this->GetPfoDirection(pPfo);
        downProbability = -1;
        downProbability = fitResult.GetProbability();
    }
    catch (...)
    {
    }
}

//----------------------------------------------------------------------------------------------------------

TrackDirectionTool::DirectionFitObject TrackDirectionTool::GetClusterDirection(const Cluster *const pTargetClusterW)
{
    try
    {

        if (LArClusterHelper::GetClusterHitType(pTargetClusterW) != TPC_VIEW_W)
        {
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        DirectionFitObject finalDirectionFitObject;

        this->AddToSlidingFitCache(pTargetClusterW);
        this->GetCalorimetricDirection(pTargetClusterW, finalDirectionFitObject);
        this->ComputeProbability(finalDirectionFitObject);
        this->SetEndpoints(finalDirectionFitObject, pTargetClusterW);
        ;

        return finalDirectionFitObject;
    }
    catch (StatusCodeException &statusCodeException)
    {
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackDirectionTool::DirectionFitObject TrackDirectionTool::GetPfoDirection(const pandora::ParticleFlowObject *const pPfo)
{

    try
    {
        const pandora::Vertex *const pVertex = LArPfoHelper::GetVertex(pPfo);
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        LArTrackStateVector trackStateVector;
        LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, m_slidingFitWindow, slidingFitPitch, trackStateVector);

        const Cluster *const pClusterW = GetTargetClusterFromPFO(pPfo, trackStateVector);

        CaloHitList totalcalohits;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, totalcalohits);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, totalcalohits);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, totalcalohits);

        if (pClusterW->GetNCaloHits() <= m_minClusterCaloHits)
        {
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }

        DirectionFitObject finalDirectionFitObject = GetClusterDirection(pClusterW);
        this->ComputeProbability(finalDirectionFitObject);

        //If the PFO is 3D, then 3D endpoints should be set
        if (LArPfoHelper::IsThreeD(pPfo))
            SetEndpoints(finalDirectionFitObject, trackStateVector);

        return finalDirectionFitObject;
    }

    catch (StatusCodeException &statusCodeException)
    {
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::WriteLookupTableToTree(LookupTable &lookupTable)
{
    std::vector<int> mapVector1, reverseMapVector2;
    std::vector<double> mapVector2, reverseMapVector1;

    for (auto &element : lookupTable.GetMap())
    {
        mapVector1.push_back(element.first);
        mapVector2.push_back(element.second);
    }

    for (auto &element : lookupTable.GetReverseMap())
    {
        reverseMapVector1.push_back(element.first);
        reverseMapVector2.push_back(element.second);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mapVector1", &mapVector1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mapVector2", &mapVector2));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "reverseMapVector1", &reverseMapVector1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "reverseMapVector2", &reverseMapVector2));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "maxRange", lookupTable.GetMaxRange()));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
}

//-----------------------------------------------------------------------------------------------------------------------

const Cluster *TrackDirectionTool::GetTargetClusterFromPFO(const ParticleFlowObject *pPfo, const LArTrackStateVector &trackStateVector)
{
    ClusterList clusterListW;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterListW);

    if (trackStateVector.size() != 0)
    {
        TrackState firstTrackState(*(trackStateVector.begin())), lastTrackState(trackStateVector.back());
        const pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
        const pandora::CartesianVector endPosition(lastTrackState.GetPosition());

        const pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
        const pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

        float currentEndpointDistance(std::numeric_limits<float>::max());
        ClusterList bestClusterList;

        for (const Cluster *const pCluster : clusterListW)
        {
            CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
            LArClusterHelper::GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);

            const pandora::CartesianVector lowZClusterVector(innerCoordinate.GetZ() < outerCoordinate.GetZ() ? innerCoordinate : outerCoordinate);
            const pandora::CartesianVector highZClusterVector(innerCoordinate.GetZ() > outerCoordinate.GetZ() ? innerCoordinate : outerCoordinate);

            if (innerCoordinate.GetY() != 0 || outerCoordinate.GetY() != 0)
                continue;

            float endpointDistance(std::abs(lowZVector.GetZ() - lowZClusterVector.GetZ()));

            if (endpointDistance < currentEndpointDistance)
            {
                currentEndpointDistance = endpointDistance;
                bestClusterList.clear();
                bestClusterList.push_back(pCluster);
            }
        }

        if (bestClusterList.size() == 0)
        {
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }

        const Cluster *const pCluster(*(bestClusterList.begin()));
        return pCluster;
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_FAILURE);
        return 0;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetEndpoints(DirectionFitObject &fitResult, const Cluster *const pCluster)
{
    CartesianVector lowZVector(0.f, 0.f, 0.f), highZVector(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pCluster, lowZVector, highZVector);

    fitResult.SetBeginpoint(lowZVector);
    fitResult.SetEndpoint(highZVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetEndpoints(DirectionFitObject &fitResult, const LArTrackStateVector &trackStateVector)
{
    TrackState firstTrackState(*(trackStateVector.begin())), lastTrackState(trackStateVector.back());
    const pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
    const pandora::CartesianVector endPosition(lastTrackState.GetPosition());

    const pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
    const pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

    fitResult.SetBeginpoint(lowZVector);
    fitResult.SetEndpoint(highZVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FillHitChargeVector(const Cluster *const pCluster, HitChargeVector &hitChargeVector)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    const TwoDSlidingFitResult &slidingFit(this->GetCachedSlidingFit(pCluster));

    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CaloHit *const pCaloHit(*hitIter);
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        float hitWidth(pCaloHit->GetCellSize1());

        float caloHitEnergy(pCaloHit->GetInputEnergy());
        caloHitEnergy *= 273.5;          //ADC to electron
        caloHitEnergy *= 23.6 / 1000000; //ionisation energy per electron in MeV
        caloHitEnergy /= 0.62;

        float rL(0.f), rT(0.f);
        slidingFit.GetLocalPosition(caloHitPosition, rL, rT);
        if (rL == 0.)
            continue;

        float calibratedUncertainty(
            std::sqrt((0.00419133 * (caloHitEnergy / hitWidth) * (caloHitEnergy / hitWidth)) + (0.00967141 * (caloHitEnergy / hitWidth)))); //70%
        HitCharge hitCharge(pCaloHit, rL, hitWidth, caloHitEnergy, calibratedUncertainty);
        hitChargeVector.push_back(hitCharge);
    }

    std::sort(hitChargeVector.begin(), hitChargeVector.end(), SortHitChargeVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TrackInnerFilter(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector)
{
    //Fill endpoint protected area into filtered vector and put all other hits in a separate vector
    float endpointProtectionRange(0.05);
    filteredHitChargeVector.insert(filteredHitChargeVector.begin(), hitChargeVector.begin(),
        hitChargeVector.begin() + endpointProtectionRange * hitChargeVector.size());
    filteredHitChargeVector.insert(filteredHitChargeVector.begin(),
        hitChargeVector.begin() + (1.0 - endpointProtectionRange) * hitChargeVector.size(), hitChargeVector.end());

    HitChargeVector innerHitChargeVector(hitChargeVector.begin() + endpointProtectionRange * hitChargeVector.size(),
        hitChargeVector.begin() + (1.0 - endpointProtectionRange) * hitChargeVector.size());

    int nNeighboursToConsider(5);
    this->SetNearestNeighbourValues(innerHitChargeVector, nNeighboursToConsider);

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortByDistanceToNN);
    filteredHitChargeVector.insert(filteredHitChargeVector.begin(), innerHitChargeVector.begin(),
        innerHitChargeVector.begin() + 0.72 * innerHitChargeVector.size()); //lots of testing has been done to optimise percentage
    std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetNearestNeighbourValues(HitChargeVector &innerHitChargeVector, int &nNeighboursToConsider)
{
    float trackLength(0.f);
    this->GetTrackLength(innerHitChargeVector, trackLength);

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortHitChargeVectorByChargeOverWidth);
    float ChargeOverWidthRange(
        (*(std::prev(innerHitChargeVector.end(), 1))).GetChargeOverWidth() - (*innerHitChargeVector.begin()).GetChargeOverWidth());

    for (HitCharge &hitCharge1 : innerHitChargeVector)
    {
        std::vector<float> distancesToNN;

        for (HitCharge &hitCharge2 : innerHitChargeVector)
        {
            if (&hitCharge1 == &hitCharge2)
                continue;

            float ChargeOverWidthDistance(
                (trackLength / ChargeOverWidthRange) * (std::abs(hitCharge1.GetChargeOverWidth() - hitCharge2.GetChargeOverWidth())));
            float Ldistance(std::abs(hitCharge1.GetLongitudinalPosition() - hitCharge2.GetLongitudinalPosition()));
            float distanceToNN(std::sqrt(ChargeOverWidthDistance * ChargeOverWidthDistance + Ldistance * Ldistance));

            distancesToNN.push_back(distanceToNN);
        }

        std::sort(distancesToNN.begin(), distancesToNN.end());
        float nearestNeighboursDistanceSum(std::accumulate(distancesToNN.begin(), distancesToNN.begin() + nNeighboursToConsider, 0.f));
        hitCharge1.SetDistanceToNN(nearestNeighboursDistanceSum);
    }

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortHitChargeVectorByRL);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SimpleTrackEndFilter(HitChargeVector &hitChargeVector)
{

    float lowerBound(0.9), upperBound(2.2);

    while (((*(hitChargeVector.begin())).GetChargeOverWidth() / (*(std::next(hitChargeVector.begin(), 1))).GetChargeOverWidth() <= lowerBound ||
               (*(hitChargeVector.begin())).GetChargeOverWidth() / (*(std::next(hitChargeVector.begin(), 1))).GetChargeOverWidth() >= upperBound) &&
           (hitChargeVector.size() > 1))
        hitChargeVector.erase(hitChargeVector.begin());

    while (((*(std::prev(hitChargeVector.end(), 1))).GetChargeOverWidth() / (*(std::prev(hitChargeVector.end(), 2))).GetChargeOverWidth() <= lowerBound ||
               (*(std::prev(hitChargeVector.end(), 1))).GetChargeOverWidth() / (*(std::prev(hitChargeVector.end(), 2))).GetChargeOverWidth() >= upperBound) &&
           (hitChargeVector.size() > 1))
        hitChargeVector.pop_back();

    for (HitChargeVector::const_iterator iter = hitChargeVector.begin(); iter != hitChargeVector.end();)
    {
        if ((*iter).m_intails == true && hitChargeVector.size() > 1)
        {
            iter = hitChargeVector.erase(iter);
        }
        else
        {
            ++iter;
        }
    }

    //Get track length and Q over W span for last step
    float trackLength(0.f), minQoverW(1e6), maxQoverW(0.f);
    this->GetTrackLength(hitChargeVector, trackLength);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetChargeOverWidth() < minQoverW)
            minQoverW = hitCharge.GetChargeOverWidth();

        if (hitCharge.GetChargeOverWidth() > maxQoverW)
            maxQoverW = hitCharge.GetChargeOverWidth();
    }

    //If there is only one hit in a wide charge range, remove it
    for (HitChargeVector::const_iterator iter = hitChargeVector.begin(); iter != hitChargeVector.end();)
    {
        bool nearbyCharge(false);

        for (HitCharge &hitCharge : hitChargeVector)
        {
            if (std::abs(hitCharge.GetLongitudinalPosition() - (*iter).GetLongitudinalPosition()) <= 0.025 * trackLength)
                continue;

            if (std::abs(hitCharge.GetChargeOverWidth() - (*iter).GetChargeOverWidth()) <= 0.1 * (maxQoverW - minQoverW))
            {
                nearbyCharge = true;
                break;
            }
        }

        if (!nearbyCharge && hitChargeVector.size() > 1)
        {
            iter = hitChargeVector.erase(iter);
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SplitHitCollectionBySize(
    HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &smallHitChargeVector, HitChargeVector &largeHitChargeVector)
{
    HitChargeVector leftHitCollection, rightHitCollection;

    for (const HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() <= splitPosition)
            leftHitCollection.push_back(hitCharge);
        else
            rightHitCollection.push_back(hitCharge);
    }

    if (leftHitCollection.size() <= rightHitCollection.size())
    {
        smallHitChargeVector = leftHitCollection;
        largeHitChargeVector = rightHitCollection;
    }
    else
    {
        smallHitChargeVector = rightHitCollection;
        largeHitChargeVector = leftHitCollection;
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SplitHitCollectionByLeftRight(
    HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &leftHitChargeVector, HitChargeVector &rightHitChargeVector)
{
    for (const HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() <= splitPosition)
            leftHitChargeVector.push_back(hitCharge);
        else
            rightHitChargeVector.push_back(hitCharge);
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetTrackLength(HitChargeVector &hitChargeVector, float &trackLength)
{
    trackLength = 0.f;

    for (HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() > trackLength)
            trackLength = hitCharge.GetLongitudinalPosition();
    }
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetAverageQoverWTrackBody(HitChargeVector &hitChargeVector, float &averageChargeTrackBody)
{
    HitChargeVector tempHitChargeVector;

    for (HitCharge &hitCharge : hitChargeVector)
        tempHitChargeVector.push_back(hitCharge);

    std::sort(tempHitChargeVector.begin(), tempHitChargeVector.end(), SortHitChargeVectorByChargeOverWidth);

    int nEntries(0);

    for (int q = 0; q < tempHitChargeVector.size(); q++)
    {
        if (q <= 0.1 * tempHitChargeVector.size() || q >= 0.6 * tempHitChargeVector.size())
            continue;

        averageChargeTrackBody += tempHitChargeVector.at(q).GetChargeOverWidth();
        nEntries++;
    }

    averageChargeTrackBody /= nEntries;
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetQoverWRange(HitChargeVector &hitChargeVector, float &QoverWRange)
{
    HitChargeVector tempHitChargeVector;

    for (HitCharge &hitCharge : hitChargeVector)
        tempHitChargeVector.push_back(hitCharge);

    std::sort(tempHitChargeVector.begin(), tempHitChargeVector.end(), SortHitChargeVectorByChargeOverWidth);

    float minQoverW(1e6), maxQoverW(0.f);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetChargeOverWidth() > maxQoverW)
            maxQoverW = hitCharge.GetChargeOverWidth();

        if (hitCharge.GetChargeOverWidth() < minQoverW)
            minQoverW = hitCharge.GetChargeOverWidth();
    }

    QoverWRange = (maxQoverW - minQoverW);
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FitHitChargeVector(HitChargeVector &hitChargeVector, TrackDirectionTool::DirectionFitObject &fitResult, int numberHitsToConsider)
{
    float particleForwardsChiSquared(0.f), particleBackwardsChiSquared(0.f);
    int numberHits(std::min(2 * numberHitsToConsider, (int)hitChargeVector.size())), particleForwardsFitStatus(-1), particleBackwardsFitStatus(-1);
    HitChargeVector forwardsFitPoints, backwardsFitPoints;

    if (hitChargeVector.size() != 0)
    {
        this->PerformFits(hitChargeVector, forwardsFitPoints, backwardsFitPoints, numberHitsToConsider, particleForwardsChiSquared,
            particleBackwardsChiSquared, particleForwardsFitStatus, particleBackwardsFitStatus);
    }

    float mean_dEdx(0.f);
    HitChargeVector thisHitChargeVector = hitChargeVector;
    for (HitCharge &hitCharge : thisHitChargeVector)
        mean_dEdx += hitCharge.GetChargeOverWidth();
    mean_dEdx /= thisHitChargeVector.size();

    std::sort(thisHitChargeVector.begin(), thisHitChargeVector.end(), SortHitChargeVectorByRL);
    std::sort(forwardsFitPoints.begin(), forwardsFitPoints.end(), SortHitChargeVectorByRL);
    std::sort(backwardsFitPoints.begin(), backwardsFitPoints.end(), SortHitChargeVectorByRL);

    DirectionFitObject finalDirectionFitObject(thisHitChargeVector, forwardsFitPoints, backwardsFitPoints, numberHits, mean_dEdx,
        particleForwardsChiSquared, particleBackwardsChiSquared);

    SplitObject tefObject(fitResult.GetTEFObject());
    finalDirectionFitObject.SetTEFObject(tefObject);

    fitResult = finalDirectionFitObject;
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FitHitChargeVector(HitChargeVector &hitChargeVector1, HitChargeVector &hitChargeVector2,
    TrackDirectionTool::DirectionFitObject &fitResult1, TrackDirectionTool::DirectionFitObject &fitResult2, int numberHitsToConsider)
{
    this->FitHitChargeVector(hitChargeVector1, fitResult1, numberHitsToConsider);
    this->FitHitChargeVector(hitChargeVector2, fitResult2, numberHitsToConsider);
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ComputeProbability(DirectionFitObject &fitResult)
{
    float deltaChiSquaredPerHit(fitResult.GetDeltaChiSquaredPerHit());

    if (deltaChiSquaredPerHit < 0.0 && deltaChiSquaredPerHit > -40.0)
    {
        float beta = 0.01;
        float alpha = 0.30465 + (0.00082051 * fitResult.GetNHits()) - (0.12857 * fitResult.GetMinChiSquaredPerHit());
        if (alpha < 0)
        {
            alpha = 0.000000000000000000000000000000001;
        }
        float pmax = 0.90706 + (0.00011538 * fitResult.GetNHits()) - (0.032143 * fitResult.GetMinChiSquaredPerHit());
        float x = std::abs(deltaChiSquaredPerHit);
        float paa = (1.0 / (2.0 * alpha));
        float pab = ((2.0 * alpha * pmax) + (2.0 * beta * pmax) - alpha - beta);
        float pac = (pow(((alpha + beta) / beta), beta / alpha));
        float p0 = 0.5 + paa * ((pab) * (pac));
        float pc = 0.5 + (p0 - 0.5) * (1.0 - exp(-alpha * x)) * (exp(-beta * x));

        fitResult.SetProbability(pc);
    }
    else
    {
        float probability(0.5);
        fitResult.SetProbability(probability);
        return;
    }
}

//---------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::PerformFits(HitChargeVector &hitChargeVector, HitChargeVector &forwardsFitPoints, HitChargeVector &backwardsFitPoints,
    int numberHitsToConsider, float &forwardsChiSquared, float &backwardsChiSquared, int & /*fitStatus1*/, int & /*fitStatus2*/)
{

    lar_content::TrackDirectionTool::LookupTable globalMuonLookupTable;
    if (globalMuonLookupTable.GetMap().empty())
    {
        globalMuonLookupTable.SetInitialEnergy(m_tableInitialEnergy);
        globalMuonLookupTable.SetBinWidth(m_tableStepSize);

        FillLookupTable(globalMuonLookupTable, 105.7);
    }

    double TotalCharge = (0.f);
    double TotalHitWidth = (0.f);
    float trackLength(0.f);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        TotalHitWidth += hitCharge.GetHitWidth();
        TotalCharge += hitCharge.GetCharge();
    }

    this->GetTrackLength(hitChargeVector, trackLength);

    lar_content::TrackDirectionTool::HitChargeVector binnedHitChargeVector;
    BinHitChargeVector(hitChargeVector, binnedHitChargeVector);

    //forwards fit
    double particleMass(105.7);
    double maxScale = 0;
    if (trackLength != 0)
    {
        maxScale = globalMuonLookupTable.GetMaxRange() / trackLength;
        if (maxScale < 1.1)
        {
            maxScale = 1.1;
        }
    }
    LookupTable lookupTable = globalMuonLookupTable;

    const int nParameters = 3;
    const std::string parName[nParameters] = {"ENDENERGY", "SCALE", "EXTRA"};
    const double vstart[nParameters] = {2.1, 1.0, 1.0};
    const double step[nParameters] = {0.5, 5, 0.5};
    const double highphysbound[nParameters] = {25, maxScale, 1.0e1};
    std::list<double> chiSquaredList;
    std::vector<double> p0list;
    std::vector<double> p1list;
    std::vector<double> p2list;

    double M = 105.7;
    for (double p0 = vstart[0]; p0 < highphysbound[0]; p0 = p0 + step[0])
    {
        double Ee(p0);
        double Le(GetLengthfromEnergy(lookupTable, Ee));
        for (double p1 = vstart[1]; p1 < highphysbound[1]; p1 = p1 + step[1])
        {
            double L(p1 * trackLength);
            double Ls(Le - L);
            double Es(GetEnergyfromLength(lookupTable, Ls));
            double alpha((Es - Ee) / TotalCharge), beta(L / TotalHitWidth);
            for (double p2 = vstart[2]; p2 < highphysbound[2]; p2 = p2 + step[2])
            {

                double chiSquared(0.0);

                for (lar_content::TrackDirectionTool::HitCharge hitCharge : binnedHitChargeVector)
                {
                    double L_i(Ls + (p1 * hitCharge.GetLongitudinalPosition())); //length
                    double E_i(GetEnergyfromLength(lookupTable, L_i));           //energy

                    double dEdx_2D(p2 * (beta / alpha) * BetheBloch(E_i, M)); //energy per length, calculated
                    double ChargeOverWidth(hitCharge.GetChargeOverWidth());   //energy per length, experimental

                    chiSquared +=
                        ((ChargeOverWidth - dEdx_2D) * (ChargeOverWidth - dEdx_2D)) / (hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
                }

                if (chiSquaredList.size() > 0)
                {
                    if (chiSquared < chiSquaredList.back())
                    {
                        chiSquaredList.push_back(chiSquared);
                        p0list.push_back(p0);
                        p1list.push_back(p1);
                        p2list.push_back(p2);
                        chiSquaredList.pop_front();
                        p0list.erase(p0list.begin());
                        p1list.erase(p1list.begin());
                        p2list.erase(p2list.begin());
                    }
                }
                else
                {
                    chiSquaredList.push_back(chiSquared);
                    p0list.push_back(p0);
                    p1list.push_back(p1);
                    p2list.push_back(p2);
                }
            }
        }
    }

    int indexvalue = 0;
    double outpar[nParameters];

    outpar[0] = p0list[indexvalue];
    outpar[1] = p1list[indexvalue];
    outpar[2] = p2list[indexvalue];

    //backwards fit
    //--------------------------------------------------------------------
    const int nParameters2 = 3;
    const std::string parName2[nParameters2] = {"ENDENERGY", "SCALE", "EXTRA"};
    const double vstart2[nParameters2] = {2.1, 1.0, 1.0};
    const double step2[nParameters] = {0.5, 5, 0.5};
    const double highphysbound2[nParameters2] = {25, maxScale, 1.0e1};
    std::list<double> chiSquaredList2;
    std::vector<double> p0list2;
    std::vector<double> p1list2;
    std::vector<double> p2list2;

    for (double p02 = vstart2[0]; p02 < highphysbound2[0]; p02 = p02 + step2[0])
    {
        double Ee(p02);
        double Le(GetLengthfromEnergy(lookupTable, Ee));
        for (double p12 = vstart2[1]; p12 < highphysbound2[1]; p12 = p12 + step2[1])
        {
            double L(p12 * trackLength);
            double Ls(Le - L);
            double Es(GetEnergyfromLength(lookupTable, Ls));
            double alpha((Es - Ee) / TotalCharge), beta(L / TotalHitWidth);
            for (double p22 = vstart2[2]; p22 < highphysbound2[2]; p22 = p22 + step2[2])
            {
                double chiSquared2(0.0);
                for (lar_content::TrackDirectionTool::HitCharge hitCharge : binnedHitChargeVector)
                {
                    double L_i(Ls + (p12 * hitCharge.GetLongitudinalPosition())); //length
                    double E_i(GetEnergyfromLength(lookupTable, L_i));            //energy

                    double dEdx_2D(p22 * (beta / alpha) * BetheBloch(E_i, M)); //energy per length, calculated
                    double ChargeOverWidth(hitCharge.GetChargeOverWidth());    //energy per length, experimental

                    chiSquared2 +=
                        ((ChargeOverWidth - dEdx_2D) * (ChargeOverWidth - dEdx_2D)) / (hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
                }

                if (chiSquaredList2.size() > 0)
                {
                    if (chiSquared2 < chiSquaredList2.back())
                    {
                        chiSquaredList2.push_back(chiSquared2);
                        p0list2.push_back(p02);
                        p1list2.push_back(p12);
                        p2list2.push_back(p22);
                        chiSquaredList2.pop_front();
                        p0list2.erase(p0list2.begin());
                        p1list2.erase(p1list2.begin());
                        p2list2.erase(p2list2.begin());
                    }
                }
                else
                {
                    chiSquaredList2.push_back(chiSquared2);
                    p0list2.push_back(p02);
                    p1list2.push_back(p12);
                    p2list2.push_back(p22);
                }
            }
        }
    }

    int indexvalue2 = 0;
    double outpar2[nParameters];

    outpar2[0] = p0list2[indexvalue2];
    outpar2[1] = p1list2[indexvalue2];
    outpar2[2] = p2list2[indexvalue2];

    //--------------------------------------------------------------------------

    double f_Ee(outpar[0]), f_L(outpar[1] * trackLength);
    double f_Le(GetLengthfromEnergy(lookupTable, f_Ee));
    double f_Ls = f_Le - f_L;

    double f_Es = GetEnergyfromLength(lookupTable, f_Ls);
    double f_deltaE = f_Es - f_Ee;

    double f_alpha = f_deltaE / TotalCharge;
    double f_beta = f_L / TotalHitWidth;

    double b_Ee(outpar2[0]), b_L(outpar2[1] * trackLength);
    double b_Le(GetLengthfromEnergy(lookupTable, b_Ee));
    double b_Ls = b_Le - b_L;

    double b_Es = GetEnergyfromLength(lookupTable, b_Ls);
    double b_deltaE = b_Es - b_Ee;

    double b_alpha = b_deltaE / TotalCharge;
    double b_beta = b_L / TotalHitWidth;

    //--------------------------------------------------------------------------

    int nHitsConsidered(0);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        double f_L_i = f_Ls + (outpar[1] * hitCharge.GetLongitudinalPosition());
        double f_E_i = GetEnergyfromLength(lookupTable, f_L_i);
        double f_dEdx_2D = outpar[2] * (f_beta / f_alpha) * BetheBloch(f_E_i, particleMass);

        double b_L_i = b_Ls + (outpar2[1] * (trackLength - hitCharge.GetLongitudinalPosition()));
        double b_E_i = GetEnergyfromLength(lookupTable, b_L_i);
        double b_dEdx_2D = outpar2[2] * (b_beta / b_alpha) * BetheBloch(b_E_i, particleMass);

        double Q_fit_f(f_dEdx_2D * hitCharge.GetHitWidth());
        double Q_fit_b(b_dEdx_2D * hitCharge.GetHitWidth());

        float forwardsDelta(hitCharge.GetChargeOverWidth() - f_dEdx_2D), backwardsDelta(hitCharge.GetChargeOverWidth() - b_dEdx_2D);

        float f_sigma(std::sqrt((0.00164585 * f_dEdx_2D * f_dEdx_2D) + (0.0201838 * f_dEdx_2D))); //80%
        float b_sigma(std::sqrt((0.00164585 * b_dEdx_2D * b_dEdx_2D) + (0.0201838 * b_dEdx_2D))); //80%

        float lp(hitCharge.GetLongitudinalPosition()), hw(hitCharge.GetHitWidth());
        float f_Q_fit_f(Q_fit_f), f_Q_fit_b(Q_fit_b);
        HitCharge forwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, f_Q_fit_f, f_sigma);
        forwardsFitPoints.push_back(forwardsRecoHitCharge);
        HitCharge backwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, f_Q_fit_b, b_sigma);
        backwardsFitPoints.push_back(backwardsRecoHitCharge);

        float forwardsHitChisquared((forwardsDelta * forwardsDelta) / (f_sigma * f_sigma));
        float backwardsHitChisquared((backwardsDelta * backwardsDelta) / (b_sigma * b_sigma));

        float Q_fit_forwards(Q_fit_f), Q_fit_backwards(Q_fit_b);

        hitCharge.SetForwardsFitCharge(Q_fit_forwards);
        hitCharge.SetForwardsSigma(f_sigma);
        hitCharge.SetForwardsDelta(forwardsDelta);
        hitCharge.SetForwardsChiSquared(forwardsHitChisquared);

        hitCharge.SetBackwardsFitCharge(Q_fit_backwards);
        hitCharge.SetBackwardsSigma(b_sigma);
        hitCharge.SetBackwardsDelta(backwardsDelta);
        hitCharge.SetBackwardsChiSquared(backwardsHitChisquared);

        if (!((hitChargeVector.size() >= 2 * numberHitsToConsider) && nHitsConsidered > numberHitsToConsider &&
                nHitsConsidered < hitChargeVector.size() - numberHitsToConsider))
        {
            forwardsChiSquared += forwardsHitChisquared;
            backwardsChiSquared += backwardsHitChisquared;
        }

        nHitsConsidered++;
    }

    std::sort(forwardsFitPoints.begin(), forwardsFitPoints.end(), SortHitChargeVectorByRL);
    std::sort(backwardsFitPoints.begin(), backwardsFitPoints.end(), SortHitChargeVectorByRL);
}

//---------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetCalorimetricDirection(const Cluster *pTargetClusterW, DirectionFitObject &directionFitObject)
{
    HitChargeVector hitChargeVector;
    this->FillHitChargeVector(pTargetClusterW, hitChargeVector);

    HitChargeVector filteredHitChargeVector;
    this->TrackInnerFilter(hitChargeVector, filteredHitChargeVector);
    this->SimpleTrackEndFilter(filteredHitChargeVector);

    if (pTargetClusterW->GetNCaloHits() < 1.5 * m_minClusterCaloHits || LArClusterHelper::GetLength(pTargetClusterW) < m_minClusterLength)
    {
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    this->FitHitChargeVector(filteredHitChargeVector, directionFitObject);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::AddToSlidingFitCache(const Cluster *const pCluster)
{

    TwoDSlidingFitResultMap slidingFitResultMap;
    if (slidingFitResultMap.size() != 0)
    {
        std::unordered_map<const pandora::Cluster *, lar_content::TwoDSlidingFitResult>::const_iterator got = slidingFitResultMap.find(pCluster);
        if (got != slidingFitResultMap.end())
        {
            return;
        }
    }

    if (slidingFitResultMap.find(pCluster) != slidingFitResultMap.end())
    {
        return;
    }

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFit(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFit)).second)
    {
        std::cout << "Sliding fit failure" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &TrackDirectionTool::GetCachedSlidingFit(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortHitChargeVectorByRL(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetLongitudinalPosition() < hitCharge2.GetLongitudinalPosition();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortHitChargeVectorByChargeOverWidth(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetChargeOverWidth() < hitCharge2.GetChargeOverWidth();
}

//----------------------------------------------------------------------------------------------------------------------------------

bool TrackDirectionTool::SortByDistanceToNN(HitCharge &hitCharge1, HitCharge &hitCharge2)
{
    return hitCharge1.GetDistanceToNN() < hitCharge2.GetDistanceToNN();
}

//----------------------------------------------------------------------------------------------------------------------------------
void TrackDirectionTool::BinHitChargeVector(HitChargeVector &hitChargeVector, HitChargeVector &binnedHitChargeVector)
{
    float binSize = (hitChargeVector.size() > 50 ? (0.5 + (hitChargeVector.size() - 50) * 2.5 / 300) : 0.0);

    if (binSize <= 0.0)
    {
        binnedHitChargeVector = hitChargeVector;
        return;
    }
    else if (binSize > 4.0)
        binSize = 4.0;

    float trackLength(0.f);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() > trackLength)
            trackLength = hitCharge.GetLongitudinalPosition();
    }

    for (float i = binSize; i <= trackLength; i += binSize)
    {
        int nHitsBin(0);
        float meanBinPosition(0.f), meanBinWidth(0.f);
        float meanBinCharge(0.f), sumSquaredSigmas(0.f);

        for (HitCharge &hitCharge : hitChargeVector)
        {
            if (hitCharge.GetLongitudinalPosition() > i)
                break;

            if (hitCharge.GetLongitudinalPosition() < i && hitCharge.GetLongitudinalPosition() >= (i - binSize))
            {
                nHitsBin++;
                meanBinPosition += hitCharge.GetLongitudinalPosition();
                meanBinWidth += hitCharge.GetHitWidth();

                meanBinCharge += hitCharge.GetCharge();
                sumSquaredSigmas += (hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
            }
        }

        if (nHitsBin == 0)
            continue;

        meanBinPosition /= nHitsBin;
        meanBinWidth /= nHitsBin;
        meanBinCharge /= nHitsBin;
        sumSquaredSigmas /= (nHitsBin * nHitsBin);

        float binUncertainty(std::sqrt(sumSquaredSigmas));
        HitCharge binnedHitCharge(NULL, meanBinPosition, meanBinWidth, meanBinCharge, binUncertainty);
        if (binnedHitCharge.GetCharge() > 0.1)
            binnedHitChargeVector.push_back(binnedHitCharge);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

double TrackDirectionTool::DensityCorrection(double &T, double &M)
{

    double p = std::sqrt((T * T) + 2 * T * M);
    double gamma = std::sqrt(1 + ((p / M) * (p / M)));
    double beta = std::sqrt(1 - 1 / (gamma * gamma));
    double X = std::log10(beta * gamma);
    const double C = 5.2146;
    const double a = 0.19559;
    const double m = 3.0;
    const double X1 = 3.0;
    const double X0 = 0.2000;
    const double delta0 = 0.0;

    if (X < X0)
        return delta0;
    else if ((X > X0) && (X < X1))
        return 2 * X * std::log(10) - C + (a * (std::pow((X1 - X), m)));
    else
        return 2 * X * std::log(10) + C;
}

//----------------------------------------------------------------------------------------------------------------------------------

double TrackDirectionTool::BetheBloch(double &T, double &M)
{
    const double K = 0.307075; // constant K in MeV cm mol^-1
    const double z = 1;        // charge in e
    const double Z = 18;       // Atomic number Z
    const double A = 39.948;   // Atomic mass in g mol-1
    const double m_e = 0.511;  // Mass of electron in MeV
    const double rho = 1.396;  // Density of material in g cm^-3 (here: argon density)
    const double I = 0.000188; // Ionisation energy in MeV

    double p(std::sqrt((T * T) + 2 * T * M));
    double gamma(std::sqrt(1 + ((p / M) * (p / M))));
    double beta(std::sqrt(1 - 1 / (gamma * gamma)));

    double T_max(2 * m_e * (beta * gamma * beta * gamma) / (1 + 2 * gamma * m_e / M + ((m_e / M) * (m_e / M))));
    return rho * ((K * z * z * Z) / A) *
           (0.5 * std::log(2 * m_e * T_max * (beta * gamma * beta * gamma) / (I * I)) - (beta * beta) - (0.5 * DensityCorrection(p, M))) /
           (beta * beta); //in MeV/cm
}
//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FillLookupTable(LookupTable &lookupTable, double M)
{
    std::map<int, double> lookupMap;
    std::map<double, int> reverseLookupMap;

    double currentEnergy(lookupTable.GetInitialEnergy()), binWidth(lookupTable.GetBinWidth());
    int maxBin(0);

    for (double n = 0; n < 100000; ++n)
    {
        double currentdEdx = BetheBloch(currentEnergy, M);

        if ((currentdEdx * binWidth) >= currentEnergy)
        {
            double maxRange = (n * binWidth) + (currentEnergy / currentdEdx);
            lookupTable.SetMaxRange(maxRange);
            maxBin = n;

            lookupMap.insert(std::pair<int, double>(n, 0.0));
            reverseLookupMap.insert(std::pair<double, int>(0.0, n));
            break;
        }
        else
        {
            lookupMap.insert(std::pair<int, double>(n, currentEnergy));
            reverseLookupMap.insert(std::pair<double, int>(currentEnergy, n));
        }

        currentEnergy -= (currentdEdx * binWidth);
    }

    //remove redundant entries to make lookup much faster
    for (std::map<int, double>::iterator it = lookupMap.begin(); it != lookupMap.end(); it++)
    {
        double n(it->first);
        double val(it->second);
        double distanceFromMaxRange((maxBin - n) * binWidth);

        if (n == 0 || n == maxBin)
            continue;
        else
        {
            double distanceMagnitude(floor(distanceFromMaxRange / 2.0));
            double samplingDistance((1 + distanceMagnitude) * binWidth);

            if (!(remainder(distanceFromMaxRange, samplingDistance) == 0.0))
            {
                lookupMap.erase(n);
                reverseLookupMap.erase(val);
            }
        }
    }

    lookupTable.SetMap(lookupMap);
    lookupTable.SetReverseMap(reverseLookupMap);

    if (lookupTable.GetMaxRange() == 0.f)
        std::cout << "Warning: the lookup table max range has not been correctly set." << std::endl;
}
//----------------------------------------------------------------------------------------------------------------------------------

double TrackDirectionTool::GetEnergyfromLength(LookupTable &lookupTable, double &trackLength)
{
    std::map<int, double> lookupMap(lookupTable.GetMap());
    double binWidth(lookupTable.GetBinWidth());

    if (trackLength >= lookupTable.GetMaxRange())
        return 0.5; //0 energy means infinite dE/dx
    else if (trackLength <= 1.0)
        return lookupTable.GetInitialEnergy();

    int n(std::floor(trackLength / binWidth));
    std::map<int, double>::iterator nextEntryIterator(lookupMap.upper_bound(n));
    std::map<int, double>::iterator previousEntryIterator(std::prev(nextEntryIterator, 1));

    double leftLength(previousEntryIterator->first * binWidth), rightLength(nextEntryIterator->first * binWidth);
    double leftEnergy(previousEntryIterator->second), rightEnergy(nextEntryIterator->second);
    double lengthDifference(rightLength - leftLength);
    double energyDifference(leftEnergy - rightEnergy);

    double finalEnergy(leftEnergy - (((trackLength - leftLength) / (lengthDifference)) * (energyDifference)));

    //very small energy values lead to huge dE/dx values: truncate
    if (finalEnergy <= 2.0)
        return 2.0;
    else
        return finalEnergy;
}

//----------------------------------------------------------------------------------------------------------------------------------

double TrackDirectionTool::GetLengthfromEnergy(LookupTable &lookupTable, double &currentEnergy)
{
    std::map<double, int> reverseLookupMap(lookupTable.GetReverseMap());
    double binWidth(lookupTable.GetBinWidth());

    if (currentEnergy <= 0.0)
        return lookupTable.GetMaxRange();
    else if (currentEnergy >= lookupTable.GetInitialEnergy())
        return 0.0;

    std::map<double, int>::iterator nextEntryIterator(reverseLookupMap.upper_bound(currentEnergy));
    std::map<double, int>::iterator previousEntryIterator(std::prev(nextEntryIterator, 1));

    double upperEnergy(nextEntryIterator->first), lowerEnergy(previousEntryIterator->first);
    double leftLength(nextEntryIterator->second * binWidth), rightLength(previousEntryIterator->second * binWidth);
    double lengthDifference(rightLength - leftLength);
    double energyDifference(upperEnergy - lowerEnergy);

    return leftLength + (((upperEnergy - currentEnergy) / energyDifference) * lengthDifference);
}

//----------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackDirectionTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NumberTrackEndHits", m_numberTrackEndHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnableFragmentRemoval", m_enableFragmentRemoval));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnableSplitting", m_enableSplitting));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTable", m_writeTable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_lookupTableFileName));

    return STATUS_CODE_SUCCESS;
}

//----------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
