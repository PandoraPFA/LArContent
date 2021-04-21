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
    m_slidingFitWindow(20),
    m_minClusterCaloHits(20),
    m_minClusterLength(10.f),
    m_tableInitialEnergy(2000.f),
    m_tableStepSize(0.5f),
    m_lookupTableFileName("lookuptable.root"),
    m_treeName("lookuptable"),
    m_endpointProtectionRange(0.05),
    m_nNeighboursToConsider(5),
    m_lowerBound(0.9),
    m_upperBound(2.2),
    m_endFilterMultiplierTrack(0.025),
    m_endFilterMultiplierCharge(0.1), 
    m_sigmaFitMultiplier(0.00164585),
    m_ADCToElectron(273.5),
    m_ionEPerElectron(2.36e-5),
    m_energyScale(0.62), 
    m_uncertaintyCalibration1(0.00419133),
    m_uncertaintyCalibration2(0.00967141),
    m_innerHitChargeMultiplier(0.72)
{
}

//--------------------------------------------------------------------
StatusCode TrackDirectionTool::Initialize()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
float TrackDirectionTool::FindDirections(const pandora::Algorithm *const, const pandora::ParticleFlowObject *const pPfo)
{
    float downProbability(-1);
    try
    {
        TrackDirectionTool::DirectionFitObject fitResult = this->GetPfoDirection(pPfo);
        downProbability = fitResult.GetProbability();
    }
    catch (...)
    {
      throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    return downProbability;
}

//----------------------------------------------------------------------------------------------------------

void TrackDirectionTool::GetClusterDirection(const Cluster *const pTargetClusterW, TrackDirectionTool::DirectionFitObject &finalDirectionFitObject)
{
    try
    {

        if (LArClusterHelper::GetClusterHitType(pTargetClusterW) != TPC_VIEW_W)
        {
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        this->AddToSlidingFitCache(pTargetClusterW);
        this->GetCalorimetricDirection(pTargetClusterW, finalDirectionFitObject);
        this->ComputeProbability(finalDirectionFitObject);
        this->SetEndPoints(finalDirectionFitObject, pTargetClusterW);
        
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
        const pandora::Vertex *const pVertex(LArPfoHelper::GetVertex(pPfo));
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        LArTrackStateVector trackStateVector;
        LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, m_slidingFitWindow, slidingFitPitch, trackStateVector);

        const Cluster *const pClusterW(GetTargetClusterFromPFO(pPfo, trackStateVector));

        if (pClusterW->GetNCaloHits() <= m_minClusterCaloHits)
        {
            throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
        }

        DirectionFitObject finalDirectionFitObject;
	this->GetClusterDirection(pClusterW, finalDirectionFitObject);
        this->ComputeProbability(finalDirectionFitObject);

        //If the PFO is 3D, then 3D endpoints should be set
        if (LArPfoHelper::IsThreeD(pPfo))
            SetEndPoints(finalDirectionFitObject, trackStateVector);

        return finalDirectionFitObject;
    }

    catch (StatusCodeException &statusCodeException)
    {
        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::WriteLookupTableToTree(LookupTable &lookupTable) const
{
    IntVector mapVector1, reverseMapVector2;
    FloatVector mapVector2, reverseMapVector1;

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

    if (!trackStateVector.empty())
    {
        TrackState firstTrackState(trackStateVector.front()), lastTrackState(trackStateVector.back());
        const pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
        const pandora::CartesianVector endPosition(lastTrackState.GetPosition());

        const pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
        const pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

        float currentEndpointDistance(std::numeric_limits<float>::max());
        const Cluster *pBestCluster(nullptr);

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
		pBestCluster = pCluster;
            }
        }

	if (pBestCluster == nullptr)
        {
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }

	const Cluster *const pCluster(pBestCluster);
        return pCluster;
    }
    else{
      throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetEndPoints(DirectionFitObject &fitResult, const Cluster *const pCluster)
{
    CartesianVector lowZVector(0.f, 0.f, 0.f), highZVector(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pCluster, lowZVector, highZVector);

    fitResult.SetBeginPoint(lowZVector);
    fitResult.SetEndPoint(highZVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetEndPoints(DirectionFitObject &fitResult, const LArTrackStateVector &trackStateVector)
{
    TrackState firstTrackState(trackStateVector.front()), lastTrackState(trackStateVector.back());
    const pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
    const pandora::CartesianVector endPosition(lastTrackState.GetPosition());

    const pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
    const pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

    fitResult.SetBeginPoint(lowZVector);
    fitResult.SetEndPoint(highZVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FillHitChargeVector(const Cluster *const pCluster, HitChargeVector &hitChargeVector)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    const TwoDSlidingFitResult &slidingFit(this->GetCachedSlidingFit(pCluster));

    for (const CaloHit *pCaloHit : caloHitList)
    {
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        float hitWidth(pCaloHit->GetCellSize1());

        float caloHitEnergy(pCaloHit->GetInputEnergy());
        caloHitEnergy *= m_ADCToElectron;
        caloHitEnergy *= m_ionEPerElectron;
        caloHitEnergy /= m_energyScale;

        float rL(0.f), rT(0.f);
        slidingFit.GetLocalPosition(caloHitPosition, rL, rT);
        if (rL < std::numeric_limits<float>::epsilon())
            continue;

        float calibratedUncertainty(
            std::sqrt((m_uncertaintyCalibration1 * (caloHitEnergy / hitWidth) * (caloHitEnergy / hitWidth)) + (m_uncertaintyCalibration2 * (caloHitEnergy / hitWidth)))); //70%
        HitCharge hitCharge(pCaloHit, rL, hitWidth, caloHitEnergy, calibratedUncertainty);
        hitChargeVector.push_back(hitCharge);
    }

    std::sort(hitChargeVector.begin(), hitChargeVector.end(), SortHitChargeVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::TrackInnerFilter(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector)
{
    //Fill endpoint protected area into filtered vector and put all other hits in a separate vector
    const int offset{static_cast<int>(m_endpointProtectionRange * hitChargeVector.size())};
    const int offsetBack{static_cast<int>((1.f - m_endpointProtectionRange) * hitChargeVector.size())};

    filteredHitChargeVector.insert(filteredHitChargeVector.begin(), hitChargeVector.begin(),
        hitChargeVector.begin() + offset);
    filteredHitChargeVector.insert(filteredHitChargeVector.begin(),
        hitChargeVector.begin() + offsetBack, hitChargeVector.end());

    HitChargeVector innerHitChargeVector(hitChargeVector.begin() + offset,
        hitChargeVector.begin() + offsetBack);

    this->SetNearestNeighbourValues(innerHitChargeVector, m_nNeighboursToConsider);

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortByDistanceToNN);
    filteredHitChargeVector.insert(filteredHitChargeVector.begin(), innerHitChargeVector.begin(),
        innerHitChargeVector.begin() + m_innerHitChargeMultiplier * innerHitChargeVector.size()); //lots of testing has been done to optimise percentage
    std::sort(filteredHitChargeVector.begin(), filteredHitChargeVector.end(), SortHitChargeVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::SetNearestNeighbourValues(HitChargeVector &innerHitChargeVector, const int nNeighboursToConsider)
{
    float trackLength(this->GetTrackLength(innerHitChargeVector));

    std::sort(innerHitChargeVector.begin(), innerHitChargeVector.end(), SortHitChargeVectorByChargeOverWidth);
    float ChargeOverWidthRange(innerHitChargeVector.back().GetChargeOverWidth() - innerHitChargeVector.front().GetChargeOverWidth());

    for (HitCharge &hitCharge1 : innerHitChargeVector)
    {
        FloatVector distancesToNN;

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
    for (auto iter = hitChargeVector.begin(); hitChargeVector.size() > 1;)
    {
	auto next{std::next(iter)};
	const float ratio{iter->GetChargeOverWidth() / next->GetChargeOverWidth()};    

	if (ratio <= m_lowerBound || ratio >= m_upperBound)
	  iter = hitChargeVector.erase(iter);
	else
	  break;
    }
    
    for (auto iter = hitChargeVector.end(); hitChargeVector.size() > 1; iter = hitChargeVector.end())
      {
	auto prev{std::prev(iter, 1)};
	auto prev2{std::prev(iter, 2)};
	const float ratio{prev->GetChargeOverWidth() / prev2->GetChargeOverWidth()};   

	if (ratio <= m_lowerBound || ratio >= m_upperBound)
	  hitChargeVector.erase(iter);
	else
	  break;
      }

    for (HitChargeVector::const_iterator iter = hitChargeVector.begin(); iter != hitChargeVector.end();)
    {
        if ((*iter).GetInTails() && hitChargeVector.size() > 1)
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
    this->GetTrackLength(hitChargeVector);

    for (HitCharge &hitCharge : hitChargeVector)
    {
      float chargeOverWidth = hitCharge.GetChargeOverWidth();

        if (chargeOverWidth < minQoverW)
            minQoverW = chargeOverWidth;

        if (chargeOverWidth > maxQoverW)
            maxQoverW = chargeOverWidth;
    }

    //If there is only one hit in a wide charge range, remove it
    for (HitChargeVector::const_iterator iter = hitChargeVector.begin(); iter != hitChargeVector.end();)
    {
        bool nearbyCharge(false);

        for (HitCharge &hitCharge : hitChargeVector)
        {
            if (std::abs(hitCharge.GetLongitudinalPosition() - (*iter).GetLongitudinalPosition()) <= m_endFilterMultiplierTrack * trackLength)
                continue;

            if (std::abs(hitCharge.GetChargeOverWidth() - (*iter).GetChargeOverWidth()) <= m_endFilterMultiplierCharge * (maxQoverW - minQoverW))
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
    const HitChargeVector &hitChargeVector, float splitPosition, HitChargeVector &smallHitChargeVector, HitChargeVector &largeHitChargeVector)
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
    const HitChargeVector &hitChargeVector, float splitPosition, HitChargeVector &leftHitChargeVector, HitChargeVector &rightHitChargeVector)
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

float TrackDirectionTool::GetTrackLength(const HitChargeVector &hitChargeVector) const
{
    float trackLength = 0.f;

    for (const HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() > trackLength)
            trackLength = hitCharge.GetLongitudinalPosition();
    }

    return trackLength;
}

//--------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FitHitChargeVector(HitChargeVector &hitChargeVector, TrackDirectionTool::DirectionFitObject &fitResult, int numberHitsToConsider)
{
    float particleForwardsChiSquared(0.f), particleBackwardsChiSquared(0.f);
    int numberHits(std::min(2 * numberHitsToConsider, (int)hitChargeVector.size()));
    HitChargeVector forwardsFitPoints, backwardsFitPoints;

    if (!hitChargeVector.empty())
        this->PerformFits(hitChargeVector, forwardsFitPoints, backwardsFitPoints, numberHitsToConsider, particleForwardsChiSquared, particleBackwardsChiSquared);

    float dEdxMean(0.f);
    HitChargeVector thisHitChargeVector = hitChargeVector;
    if (!thisHitChargeVector.empty()) {
      for (HitCharge &hitCharge : thisHitChargeVector)
	dEdxMean += hitCharge.GetChargeOverWidth();
      dEdxMean /= thisHitChargeVector.size();
    }

    std::sort(thisHitChargeVector.begin(), thisHitChargeVector.end(), SortHitChargeVectorByRL);
    std::sort(forwardsFitPoints.begin(), forwardsFitPoints.end(), SortHitChargeVectorByRL);
    std::sort(backwardsFitPoints.begin(), backwardsFitPoints.end(), SortHitChargeVectorByRL);

    DirectionFitObject finalDirectionFitObject(thisHitChargeVector, forwardsFitPoints, backwardsFitPoints, numberHits, dEdxMean,
        particleForwardsChiSquared, particleBackwardsChiSquared);

    fitResult = finalDirectionFitObject;
}

//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::ComputeProbability(DirectionFitObject &fitResult)
{
    float deltaChiSquaredPerHit(fitResult.GetDeltaChiSquaredPerHit());

    float maxDeltaChiSquaredPerHit(0.0);
    float minDeltaChiSquaredPerHit(-40.0);

    if (deltaChiSquaredPerHit < maxDeltaChiSquaredPerHit && deltaChiSquaredPerHit > minDeltaChiSquaredPerHit)
    {
        const float alphaFitParamConst(0.30465f);
	const float alphaFitParamN(0.00082051f);
	const float alphaFitParamChi(0.12857f);
	const float pmaxFitParamConst(0.90706f);
	const float pmaxFitParamN(0.00011538f);
	const float pmaxFitParamChi(0.032143f);

        float beta(0.01);
        float alpha = alphaFitParamConst + (alphaFitParamN * fitResult.GetNHits()) - (alphaFitParamChi * fitResult.GetMinChiSquaredPerHit());
        if (alpha < std::numeric_limits<float>::epsilon())
        {
	  alpha = 1e-33;
        }
        float pmax = pmaxFitParamConst + (pmaxFitParamN * fitResult.GetNHits()) - (pmaxFitParamChi * fitResult.GetMinChiSquaredPerHit());
        float absDeltaChiSquaredPerHit = std::abs(deltaChiSquaredPerHit);

        float p0PartA = (1.0 / (2.0 * alpha));
        float p0PartB = ((2.0 * alpha * pmax) + (2.0 * beta * pmax) - alpha - beta);
        float p0PartC = (pow(((alpha + beta) / beta), beta / alpha));
        float p0 = 0.5 + p0PartA * ((p0PartB) * (p0PartC));
        float pc = 0.5 + (p0 - 0.5) * (1.0 - exp(-alpha * absDeltaChiSquaredPerHit)) * (exp(-beta * absDeltaChiSquaredPerHit));

        fitResult.SetProbability(pc);
    }
    else
    {
      fitResult.SetProbability(0.5f);
    }
}

//---------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::PerformFits(HitChargeVector &hitChargeVector, HitChargeVector &forwardsFitPoints, HitChargeVector &backwardsFitPoints,
    int numberHitsToConsider, float &forwardsChiSquared, float &backwardsChiSquared)
{
    const float particleMass(105.7);
    LookupTable lookupTable;
    if (lookupTable.GetMap().empty())
    {
        lookupTable.SetInitialEnergy(m_tableInitialEnergy);
        lookupTable.SetBinWidth(m_tableStepSize);

        FillLookupTable(lookupTable, particleMass);
    }

    float totalCharge = (0.f);
    float totalHitWidth = (0.f);
    float trackLength(0.f);
    bool isForwardsFit(true);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        totalHitWidth += hitCharge.GetHitWidth();
        totalCharge += hitCharge.GetCharge();
    }

    trackLength = this->GetTrackLength(hitChargeVector);

    HitChargeVector binnedHitChargeVector;
    BinHitChargeVector(hitChargeVector, binnedHitChargeVector);

    //forwards fit
    float maxScale(0);
    if (trackLength > std::numeric_limits<float>::epsilon())
    {
        maxScale = lookupTable.GetMaxRange() / trackLength;
        if (maxScale < 1.1)
        {
            maxScale = 1.1;
        }
    }

    const int nParameters(3);
    float holdp0Value(0.0f);
    float holdp1Value(0.0f);
    float holdp2Value(0.0f);

    this->GetFitParameters(isForwardsFit, lookupTable, binnedHitChargeVector, trackLength, totalCharge, totalHitWidth, maxScale, particleMass, holdp0Value, holdp1Value, holdp2Value);

    float outpar[nParameters];

    outpar[0] = holdp0Value;
    outpar[1] = holdp1Value;
    outpar[2] = holdp2Value;
    
    //backwards fit
    //--------------------------------------------------------------------
    isForwardsFit = false;
    float holdp0Value2(0.0f);
    float holdp1Value2(0.0f);
    float holdp2Value2(0.0f);

    this->GetFitParameters(isForwardsFit, lookupTable, binnedHitChargeVector, trackLength, totalCharge, totalHitWidth, maxScale, particleMass, holdp0Value2, holdp1Value2, holdp2Value2);

    float outpar2[nParameters];

    outpar2[0] = holdp0Value2;
    outpar2[1] = holdp1Value2;
    outpar2[2] = holdp2Value2;

    //--------------------------------------------------------------------------

    float f_Ee(outpar[0]), f_L(outpar[1] * trackLength);
    float f_Le(GetLengthfromEnergy(lookupTable, f_Ee));
    float f_Ls = f_Le - f_L;

    float f_Es = GetEnergyfromLength(lookupTable, f_Ls);
    float f_deltaE = f_Es - f_Ee;

    float f_alpha = f_deltaE / totalCharge;
    float f_beta = f_L / totalHitWidth;

    float b_Ee(outpar2[0]), b_L(outpar2[1] * trackLength);
    float b_Le(GetLengthfromEnergy(lookupTable, b_Ee));
    float b_Ls = b_Le - b_L;

    float b_Es = GetEnergyfromLength(lookupTable, b_Ls);
    float b_deltaE = b_Es - b_Ee;

    float b_alpha = b_deltaE / totalCharge;
    float b_beta = b_L / totalHitWidth;

    //--------------------------------------------------------------------------

    int nHitsConsidered(0);

    for (HitCharge &hitCharge : hitChargeVector)
    {
        float f_L_i = f_Ls + (outpar[1] * hitCharge.GetLongitudinalPosition());
        float f_E_i = GetEnergyfromLength(lookupTable, f_L_i);
        float f_dEdx_2D = outpar[2] * (f_beta / f_alpha) * BetheBloch(f_E_i, particleMass);

        float b_L_i = b_Ls + (outpar2[1] * (trackLength - hitCharge.GetLongitudinalPosition()));
        float b_E_i = GetEnergyfromLength(lookupTable, b_L_i);
        float b_dEdx_2D = outpar2[2] * (b_beta / b_alpha) * BetheBloch(b_E_i, particleMass);

        float qFitForwards(f_dEdx_2D * hitCharge.GetHitWidth());
        float qFitBackwards(b_dEdx_2D * hitCharge.GetHitWidth());

        float forwardsDelta(hitCharge.GetChargeOverWidth() - f_dEdx_2D), backwardsDelta(hitCharge.GetChargeOverWidth() - b_dEdx_2D);

        float sigmaForwards(std::sqrt((m_sigmaFitMultiplier * f_dEdx_2D * f_dEdx_2D) + (0.0201838 * f_dEdx_2D))); //80%
        float sigmaBackwards(std::sqrt((m_sigmaFitMultiplier * b_dEdx_2D * b_dEdx_2D) + (0.0201838 * b_dEdx_2D))); //80%

        float lp(hitCharge.GetLongitudinalPosition()), hw(hitCharge.GetHitWidth());
        float f_Q_fit_f(qFitForwards), f_Q_fit_b(qFitBackwards);
        HitCharge forwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, f_Q_fit_f, sigmaForwards);
        forwardsFitPoints.push_back(forwardsRecoHitCharge);
        HitCharge backwardsRecoHitCharge(hitCharge.GetCaloHit(), lp, hw, f_Q_fit_b, sigmaBackwards);
        backwardsFitPoints.push_back(backwardsRecoHitCharge);

        float forwardsHitChiSquared((forwardsDelta * forwardsDelta) / (sigmaForwards * sigmaForwards));
        float backwardsHitChiSquared((backwardsDelta * backwardsDelta) / (sigmaBackwards * sigmaBackwards));

        float Q_fit_forwards(qFitForwards), Q_fit_backwards(qFitBackwards);

        hitCharge.SetForwardsFitCharge(Q_fit_forwards);
        hitCharge.SetForwardsSigma(sigmaForwards);
        hitCharge.SetForwardsDelta(forwardsDelta);
        hitCharge.SetForwardsChiSquared(forwardsHitChiSquared);

        hitCharge.SetBackwardsFitCharge(Q_fit_backwards);
        hitCharge.SetBackwardsSigma(sigmaBackwards);
        hitCharge.SetBackwardsDelta(backwardsDelta);
        hitCharge.SetBackwardsChiSquared(backwardsHitChiSquared);

        if (!((hitChargeVector.size() >= 2 * numberHitsToConsider) && nHitsConsidered > numberHitsToConsider &&
                nHitsConsidered < hitChargeVector.size() - numberHitsToConsider))
        {
            forwardsChiSquared += forwardsHitChiSquared;
            backwardsChiSquared += backwardsHitChiSquared;
        }

        nHitsConsidered++;
    }

    std::sort(forwardsFitPoints.begin(), forwardsFitPoints.end(), SortHitChargeVectorByRL);
    std::sort(backwardsFitPoints.begin(), backwardsFitPoints.end(), SortHitChargeVectorByRL);
}

//---------------------------------------------------------------------------------------------------------------------------------

  void TrackDirectionTool::GetFitParameters(bool isForwardsFit, LookupTable lookupTable, HitChargeVector binnedHitChargeVector, const float trackLength, const float totalCharge, const float totalHitWidth, const float maxScale, const float particleMass, float &holdp0Value, float &holdp1Value, float &holdp2Value)
{
    const int nParameters(3);
    const std::string parName[nParameters] = {"ENDENERGY", "SCALE", "EXTRA"};
    const float vstart[nParameters] = {2.1, 1.0, 1.0};
    const float step[nParameters] = {0.5, 5.0, 0.5};
    const float highphysbound[nParameters] = {25.0, maxScale, 1.0e1};
    float holdChiSquaredValue(std::numeric_limits<float>::max());

    for (float p0 = vstart[0]; p0 < highphysbound[0]; p0 = p0 + step[0])
    {
        float Ee(p0);
        float Le(GetLengthfromEnergy(lookupTable, Ee));
        for (float p1 = vstart[1]; p1 < highphysbound[1]; p1 = p1 + step[1])
        {
            float L(p1 * trackLength);
            float Ls(Le - L);
            float Es(GetEnergyfromLength(lookupTable, Ls));
            float alpha((Es - Ee) / totalCharge), beta(L / totalHitWidth);
            for (float p2 = vstart[2]; p2 < highphysbound[2]; p2 = p2 + step[2])
            {
                float chiSquared(0.0);

                for (lar_content::TrackDirectionTool::HitCharge hitCharge : binnedHitChargeVector)
                {
		    float L_i(0.0f);
		    if(isForwardsFit) {
                      L_i = (Ls + (p1 * hitCharge.GetLongitudinalPosition())); //length
		    }
		    else if(!isForwardsFit) {
		      L_i = (Ls + (p1 * (trackLength - hitCharge.GetLongitudinalPosition()))); //length
         	    }
                    float E_i(GetEnergyfromLength(lookupTable, L_i));           //energy

                    float dEdx2D(p2 * (beta / alpha) * BetheBloch(E_i, particleMass)); //energy per length, calculated
                    float ChargeOverWidth(hitCharge.GetChargeOverWidth());   //energy per length, experimental

                    chiSquared +=
                        ((ChargeOverWidth - dEdx2D) * (ChargeOverWidth - dEdx2D)) / (hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
                }
		if (chiSquared < holdChiSquaredValue)
                {
		  holdChiSquaredValue = chiSquared;
		  holdp0Value = p0;
		  holdp1Value = p1;
		  holdp2Value = p2;
		}
            }
        }
    }
    
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

float TrackDirectionTool::DensityCorrection(const float T, const float M)
{

    const float p = std::sqrt((T * T) + 2 * T * M);
    const float gamma = std::sqrt(1 + ((p / M) * (p / M)));
    const float beta = std::sqrt(1 - 1 / (gamma * gamma));
    const float X = std::log10(beta * gamma);
    const float C = 5.2146;
    const float a = 0.19559;
    const float m = 3.0;
    const float X1 = 3.0;
    const float X0 = 0.2000;
    const float delta0 = 0.0;

    if (X < X0)
        return delta0;
    else if ((X > X0) && (X < X1))
        return 2 * X * std::log(10) - C + (a * (std::pow((X1 - X), m)));
    else
        return 2 * X * std::log(10) + C;
}

//----------------------------------------------------------------------------------------------------------------------------------

float TrackDirectionTool::BetheBloch(float &T, const float M)
{
    const float K = 0.307075; // constant K in MeV cm mol^-1
    const float z = 1;        // charge in e
    const float Z = 18;       // Atomic number Z
    const float A = 39.948;   // Atomic mass in g mol-1
    const float m_e = 0.511;  // Mass of electron in MeV
    const float rho = 1.396;  // Density of material in g cm^-3 (here: argon density)
    const float I = 0.000188; // Ionisation energy in MeV

    const float p(std::sqrt((T * T) + 2 * T * M));
    const float gamma(std::sqrt(1 + ((p / M) * (p / M))));
    const float beta(std::sqrt(1 - 1 / (gamma * gamma)));

    const float T_max(2 * m_e * (beta * gamma * beta * gamma) / (1 + 2 * gamma * m_e / M + ((m_e / M) * (m_e / M))));
    return rho * ((K * z * z * Z) / A) *
           (0.5 * std::log(2 * m_e * T_max * (beta * gamma * beta * gamma) / (I * I)) - (beta * beta) - (0.5 * DensityCorrection(p, M))) /
           (beta * beta); //in MeV/cm
}
//----------------------------------------------------------------------------------------------------------------------------------

void TrackDirectionTool::FillLookupTable(LookupTable &lookupTable, const float M)
{
    std::map<int, float> lookupMap;
    std::map<float, int> reverseLookupMap;

    float currentEnergy(lookupTable.GetInitialEnergy());
    float binWidth(lookupTable.GetBinWidth());
    int maxBin(0);

    for (float n = 0; n < 100000; ++n)
    {
        float currentdEdx = BetheBloch(currentEnergy, M);

        if ((currentdEdx * binWidth) >= currentEnergy)
        {
            float maxRange = (n * binWidth) + (currentEnergy / currentdEdx);
            lookupTable.SetMaxRange(maxRange);
            maxBin = n;

            lookupMap.insert(std::pair<int, float>(n, 0.0));
            reverseLookupMap.insert(std::pair<float, int>(0.0, n));
            break;
        }
        else
        {
            lookupMap.insert(std::pair<int, float>(n, currentEnergy));
            reverseLookupMap.insert(std::pair<float, int>(currentEnergy, n));
        }

        currentEnergy -= (currentdEdx * binWidth);
    }

    //remove redundant entries to make lookup much faster
    for (std::map<int, float>::iterator it = lookupMap.begin(); it != lookupMap.end(); it++)
    {
        const float n(it->first);

        if (n == 0 || n == maxBin)
            continue;
	 else
	 {
	    const float val(it->second);
	    const float distanceFromMaxRange((maxBin - n) * binWidth);
	   
            const float distanceMagnitude(floor(distanceFromMaxRange / 2.0));
            const float samplingDistance((1 + distanceMagnitude) * binWidth);

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

float TrackDirectionTool::GetEnergyfromLength(LookupTable &lookupTable, const float &trackLength)
{
    std::map<int, float> lookupMap(lookupTable.GetMap());
    const float binWidth(lookupTable.GetBinWidth());

    if (trackLength >= lookupTable.GetMaxRange())
        return 0.5; //0 energy means infinite dE/dx
    else if (trackLength <= 1.0)
        return lookupTable.GetInitialEnergy();

    const int n(std::floor(trackLength / binWidth));
    std::map<int, float>::iterator nextEntryIterator(lookupMap.upper_bound(n));
    std::map<int, float>::iterator previousEntryIterator(std::prev(nextEntryIterator, 1));

    const float leftLength(previousEntryIterator->first * binWidth), rightLength(nextEntryIterator->first * binWidth);
    const float leftEnergy(previousEntryIterator->second), rightEnergy(nextEntryIterator->second);
    const float lengthDifference(rightLength - leftLength);
    const float energyDifference(leftEnergy - rightEnergy);

    const float finalEnergy(leftEnergy - (((trackLength - leftLength) / (lengthDifference)) * (energyDifference)));

    //very small energy values lead to huge dE/dx values: truncate
    if (finalEnergy <= 2.0)
        return 2.0;
    else
        return finalEnergy;
}

//----------------------------------------------------------------------------------------------------------------------------------

float TrackDirectionTool::GetLengthfromEnergy(LookupTable &lookupTable, const float &currentEnergy)
{
    std::map<float, int> reverseLookupMap(lookupTable.GetReverseMap());
    const float binWidth(lookupTable.GetBinWidth());

    if (currentEnergy <= 0.0)
        return lookupTable.GetMaxRange();
    else if (currentEnergy >= lookupTable.GetInitialEnergy())
        return 0.0;

    std::map<float, int>::iterator nextEntryIterator(reverseLookupMap.upper_bound(currentEnergy));
    std::map<float, int>::iterator previousEntryIterator(std::prev(nextEntryIterator, 1));

    const float upperEnergy(nextEntryIterator->first), lowerEnergy(previousEntryIterator->first);
    const float leftLength(nextEntryIterator->second * binWidth), rightLength(previousEntryIterator->second * binWidth);
    const float lengthDifference(rightLength - leftLength);
    const float energyDifference(upperEnergy - lowerEnergy);

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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_lookupTableFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EndpointProtectionRange", m_endpointProtectionRange));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "nNeighboursToConsider", m_nNeighboursToConsider));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LowerBound", m_lowerBound));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UpperBound", m_upperBound));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EndFilterMultiplierTrack", m_endFilterMultiplierTrack));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EndFilterMultiplierCharge", m_endFilterMultiplierCharge));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SigmaFitMultiplier", m_sigmaFitMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ADCToElectron", m_ADCToElectron));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IonEPerElectron", m_ionEPerElectron));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScale", m_energyScale));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UncertaintyCalibration1", m_uncertaintyCalibration1));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UncertaintyCalibration2", m_uncertaintyCalibration2));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InnerHitChargeMultiplier", m_innerHitChargeMultiplier));

    return STATUS_CODE_SUCCESS;
}

//----------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
