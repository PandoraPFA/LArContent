/**
 *  @file   larpandoracontent/LArObjects/LArTwoDSlidingFitResult.cc
 *
 *  @brief  Implementation of the lar two dimensional sliding fit result class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

TwoDSlidingFitResult::TwoDSlidingFitResult(const Cluster *const pCluster, const unsigned int layerFitHalfWindow, const float layerPitch) :
    m_pCluster(pCluster),
    m_layerFitHalfWindow(layerFitHalfWindow),
    m_layerPitch(layerPitch),
    m_axisIntercept(0.f, 0.f, 0.f),
    m_axisDirection(0.f, 0.f, 0.f),
    m_orthoDirection(0.f, 0.f, 0.f)
{
    // Get a list of hits coordinates from the cluster
    CartesianPointVector coordinateVector;
    LArClusterHelper::GetCoordinateVector(pCluster, coordinateVector);

    // Calculate the sliding fit result
    this->CalculateAxes(coordinateVector);
    this->FillLayerFitContributionMap(coordinateVector);
    this->PerformSlidingLinearFit();
    this->FindSlidingFitSegments();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingFitResult::TwoDSlidingFitResult(const ClusterList &clusterList, const unsigned int layerFitHalfWindow, const float layerPitch) :
    m_layerFitHalfWindow(layerFitHalfWindow),
    m_layerPitch(layerPitch),
    m_axisIntercept(0.f, 0.f, 0.f),
    m_axisDirection(0.f, 0.f, 0.f),
    m_orthoDirection(0.f, 0.f, 0.f)
{
    // Get a list of hits coordinates from the cluster
    CartesianPointVector coordinateVector;
    for (const Cluster * const pCluster : clusterList)
    {
        CartesianPointVector clusterCoordinateVector;
        LArClusterHelper::GetCoordinateVector(pCluster, clusterCoordinateVector);
        
        coordinateVector.insert(coordinateVector.end(), std::make_move_iterator(clusterCoordinateVector.begin()), 
                                std::make_move_iterator(clusterCoordinateVector.end()));
    }

    // Calculate the sliding fit result
    this->CalculateAxes(coordinateVector);
    this->FillLayerFitContributionMap(coordinateVector);
    this->PerformSlidingLinearFit();
    this->FindSlidingFitSegments();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingFitResult::TwoDSlidingFitResult(const Cluster *const pCluster, const unsigned int layerFitHalfWindow, const float layerPitch,
        const CartesianVector &axisIntercept, const CartesianVector &axisDirection, const CartesianVector &orthoDirection) :
    m_pCluster(pCluster),
    m_layerFitHalfWindow(layerFitHalfWindow),
    m_layerPitch(layerPitch),
    m_axisIntercept(axisIntercept),
    m_axisDirection(axisDirection),
    m_orthoDirection(orthoDirection)
{
    // Get a list of hits coordinates from the cluster
    CartesianPointVector coordinateVector;
    LArClusterHelper::GetCoordinateVector(pCluster, coordinateVector);

    // Calculate the sliding fit result
    this->FillLayerFitContributionMap(coordinateVector);
    this->PerformSlidingLinearFit();
    this->FindSlidingFitSegments();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingFitResult::TwoDSlidingFitResult(const Cluster *const pCluster, const unsigned int layerFitHalfWindow, const float layerPitch,
        const CartesianVector &axisIntercept, const CartesianVector &axisDirection, const CartesianVector &orthoDirection, 
        const LayerFitContributionMap &layerFitContributionMap) :
    m_pCluster(pCluster),
    m_layerFitHalfWindow(layerFitHalfWindow),
    m_layerPitch(layerPitch),
    m_axisIntercept(axisIntercept),
    m_axisDirection(axisDirection),
    m_orthoDirection(orthoDirection),
    m_layerFitContributionMap(layerFitContributionMap)
{
    this->PerformSlidingLinearFit();
    this->FindSlidingFitSegments();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetLayerFitHalfWindowLength() const
{
    return (static_cast<float>(m_layerFitHalfWindow)) * m_layerPitch;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TwoDSlidingFitResult::GetMinLayer() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_layerFitResultMap.begin()->first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TwoDSlidingFitResult::GetMaxLayer() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_layerFitResultMap.rbegin()->first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TwoDSlidingFitResult::GetLayer(const float rL) const
{
    if (m_layerPitch < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return static_cast<int>(std::floor(rL / m_layerPitch));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetL(const int layer) const
{
    return (static_cast<float>(layer) + 0.5f) * m_layerPitch;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetLocalPosition(const CartesianVector &position, float &rL, float &rT) const
{
    const CartesianVector displacement(position - m_axisIntercept);

    rL = displacement.GetDotProduct(m_axisDirection);
    rT = displacement.GetDotProduct(m_orthoDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetLocalDirection(const CartesianVector &direction, float&dTdL) const
{
    float pL(0.f), pT(0.f);
    this->GetLocalPosition((direction + m_axisIntercept), pL, pT);

    if (std::fabs(pL) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    dTdL = pT / pL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalPosition(const float rL, const float rT, CartesianVector &position) const
{  
    position = m_axisIntercept + (m_axisDirection * rL) + (m_orthoDirection * rT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalDirection(const float dTdL, CartesianVector &direction) const
{
    const float pL(1.f / std::sqrt(1.f + dTdL * dTdL));
    const float pT(dTdL / std::sqrt(1.f + dTdL * dTdL));

    CartesianVector globalCoordinates(0.f, 0.f, 0.f);
    this->GetGlobalPosition(pL, pT, globalCoordinates);
    direction = (globalCoordinates - m_axisIntercept).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalMinLayerPosition() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin();
    CartesianVector position(0.f, 0.f, 0.f);
    this->GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);
    return position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalMaxLayerPosition() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    LayerFitResultMap::const_reverse_iterator iter = m_layerFitResultMap.rbegin();
    CartesianVector position(0.f, 0.f, 0.f);
    this->GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);
    return position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalMinLayerDirection() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin();
    CartesianVector direction(0.f, 0.f, 0.f);
    this->GetGlobalDirection(iter->second.GetGradient(), direction);
    return direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalMaxLayerDirection() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    LayerFitResultMap::const_reverse_iterator iter = m_layerFitResultMap.rbegin();
    CartesianVector direction(0.f, 0.f, 0.f);
    this->GetGlobalDirection(iter->second.GetGradient(), direction);
    return direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetMinLayerRms() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin();
    return iter->second.GetRms();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetMaxLayerRms() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    LayerFitResultMap::const_reverse_iterator iter = m_layerFitResultMap.rbegin();
    return iter->second.GetRms();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetFitRms(const float rL) const
{
    LayerInterpolation layerInterpolation;
    const StatusCode statusCode(this->LongitudinalInterpolation(rL, layerInterpolation));

    if (STATUS_CODE_SUCCESS != statusCode)
        throw StatusCodeException(statusCode);

    return this->GetFitRms(layerInterpolation);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetCosScatteringAngle(const float rL) const
{
    const float halfWindowLength(this->GetLayerFitHalfWindowLength());

    CartesianVector firstDirection(0.f,0.f,0.f);
    CartesianVector secondDirection(0.f,0.f,0.f);

    const StatusCode firstStatusCode(this->GetGlobalFitDirection(rL - halfWindowLength, firstDirection));
    const StatusCode secondStatusCode(this->GetGlobalFitDirection(rL + halfWindowLength, secondDirection));

    if (STATUS_CODE_SUCCESS != firstStatusCode)
        throw StatusCodeException(firstStatusCode);

    if (STATUS_CODE_SUCCESS != secondStatusCode)
        throw StatusCodeException(secondStatusCode);

    return firstDirection.GetDotProduct(secondDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetGlobalFitPosition(const float rL, CartesianVector &position) const
{
    LayerInterpolation layerInterpolation;
    const StatusCode statusCode(this->LongitudinalInterpolation(rL, layerInterpolation));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    position = this->GetGlobalFitPosition(layerInterpolation);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetGlobalFitDirection(const float rL, CartesianVector &direction) const
{
    LayerInterpolation layerInterpolation;
    const StatusCode statusCode(this->LongitudinalInterpolation(rL, layerInterpolation));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    direction = this->GetGlobalFitDirection(layerInterpolation);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetGlobalFitPositionAtX(const float x, CartesianVector &position) const
{
    LayerInterpolationList layerInterpolationList;
    const StatusCode statusCode(this->TransverseInterpolation(x, layerInterpolationList));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    if (layerInterpolationList.size() != 1)
        return STATUS_CODE_NOT_FOUND;

    position = this->GetGlobalFitPosition(layerInterpolationList.at(0));
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetGlobalFitDirectionAtX(const float x, CartesianVector &direction) const
{
    LayerInterpolationList layerInterpolationList;
    const StatusCode statusCode(this->TransverseInterpolation(x, layerInterpolationList));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    if (layerInterpolationList.size() != 1)
        return STATUS_CODE_NOT_FOUND;

    direction = this->GetGlobalFitDirection(layerInterpolationList.at(0));
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetGlobalFitProjection(const CartesianVector &inputPosition, CartesianVector &projectedPosition) const
{
    float rL(0.f), rT(0.f);
    this->GetLocalPosition(inputPosition, rL, rT);
    return this->GetGlobalFitPosition(rL, projectedPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetGlobalFitPositionListAtX(const float x, CartesianPointVector &positionList) const
{
    LayerInterpolationList layerInterpolationList;
    const StatusCode statusCode(this->TransverseInterpolation(x, layerInterpolationList));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    for (LayerInterpolationList::const_iterator iter = layerInterpolationList.begin(), iterEnd = layerInterpolationList.end();
        iter != iterEnd; ++iter)
    {
        positionList.push_back(this->GetGlobalFitPosition(*iter));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetTransverseProjection(const float x, const FitSegment &fitSegment, CartesianVector &position) const
{
    LayerInterpolation layerInterpolation;
    const StatusCode statusCode(this->TransverseInterpolation(x, fitSegment, layerInterpolation));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    position = this->GetGlobalFitPosition(layerInterpolation);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetTransverseProjection(const float x, const FitSegment &fitSegment, CartesianVector &position,
    CartesianVector &direction) const
{
    LayerInterpolation layerInterpolation;
    const StatusCode statusCode(this->TransverseInterpolation(x, fitSegment, layerInterpolation));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    position = this->GetGlobalFitPosition(layerInterpolation);
    direction = this->GetGlobalFitDirection(layerInterpolation);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetExtrapolatedPosition(const float rL, CartesianVector &position) const
{
    const StatusCode statusCode(this->GetGlobalFitPosition(rL, position));

    if (STATUS_CODE_NOT_FOUND != statusCode)
        return statusCode;

    const int thisLayer(this->GetLayer(rL));
    const int minLayer(this->GetMinLayer());
    const int maxLayer(this->GetMaxLayer());

    if (thisLayer <= minLayer)
    {
        position = (this->GetGlobalMinLayerPosition() + this->GetGlobalMinLayerDirection() * (rL - this->GetL(minLayer)));
    }
    else if (thisLayer >= maxLayer)
    {
        position = (this->GetGlobalMaxLayerPosition() + this->GetGlobalMaxLayerDirection() * (rL - this->GetL(maxLayer)));
    }
    else
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}   

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetExtrapolatedDirection(const float rL, CartesianVector &direction) const
{
    const StatusCode statusCode(this->GetGlobalFitDirection(rL, direction));

    if (STATUS_CODE_NOT_FOUND != statusCode)
        return statusCode;  

    const int thisLayer(this->GetLayer(rL));
    const int minLayer(this->GetMinLayer());
    const int maxLayer(this->GetMaxLayer());

    if (thisLayer <= minLayer)
    {
        direction = this->GetGlobalMinLayerDirection();
    }
    else if (thisLayer >= maxLayer)
    {
        direction = this->GetGlobalMaxLayerDirection();
    }
    else
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetExtrapolatedPositionAtX(const float x, CartesianVector &position) const
{
    const StatusCode statusCode(this->GetGlobalFitPositionAtX(x, position));

    if (STATUS_CODE_NOT_FOUND != statusCode)
        return statusCode;

    const float minLayerX(this->GetGlobalMinLayerPosition().GetX());
    const float maxLayerX(this->GetGlobalMaxLayerPosition().GetX());

    const int minLayer(this->GetMinLayer());
    const int maxLayer(this->GetMaxLayer());

    const int innerLayer((minLayerX < maxLayerX) ? minLayer : maxLayer);
    const int outerLayer((minLayerX < maxLayerX) ? maxLayer : minLayer);

    const CartesianVector innerVertex((innerLayer == minLayer) ? this->GetGlobalMinLayerPosition() : this->GetGlobalMaxLayerPosition());
    const CartesianVector innerDirection((innerLayer == minLayer) ? this->GetGlobalMinLayerDirection() * -1.f : this->GetGlobalMaxLayerDirection());

    const CartesianVector outerVertex((outerLayer == minLayer) ? this->GetGlobalMinLayerPosition() : this->GetGlobalMaxLayerPosition());
    const CartesianVector outerDirection((outerLayer == minLayer) ? this->GetGlobalMinLayerDirection() * -1.f : this->GetGlobalMaxLayerDirection());

    if (innerDirection.GetX() > -std::numeric_limits<float>::epsilon() || outerDirection.GetX() < +std::numeric_limits<float>::epsilon() ||
        outerVertex.GetX() - innerVertex.GetX() < +std::numeric_limits<float>::epsilon())
    {
        return STATUS_CODE_NOT_FOUND;
    }
    else if (x >= outerVertex.GetX())
    {
        position = outerVertex + outerDirection * ((x - outerVertex.GetX()) / outerDirection.GetX());
    }
    else if (x <= innerVertex.GetX())
    {
        position = innerVertex + innerDirection * ((x - innerVertex.GetX()) / innerDirection.GetX());
    }
    else
    {
        return STATUS_CODE_NOT_FOUND;
    }

    // TODO How to assign an uncertainty on the extrapolated position?
    return STATUS_CODE_SUCCESS;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

const FitSegment &TwoDSlidingFitResult::GetFitSegment(const float rL) const
{
    int layer(this->GetLayer(rL));

    for (FitSegmentList::const_iterator iter = m_fitSegmentList.begin(), iterEnd = m_fitSegmentList.end(); iter != iterEnd; ++iter)
    {
        const FitSegment &fitSegment = *iter;

        if (layer >= fitSegment.GetStartLayer() && layer <= fitSegment.GetEndLayer())
            return fitSegment;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

// Private member functions start here
//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::CalculateAxes(const CartesianPointVector &coordinateVector)
{
    // Use extremal coordinates to define axis intercept and direction
    CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(coordinateVector, innerCoordinate, outerCoordinate);
    m_axisIntercept = innerCoordinate;
    m_axisDirection = (outerCoordinate - innerCoordinate).GetUnitVector();

    // Use y-axis to generate an orthogonal axis (assuming that cluster occupies X-Z plane)
    CartesianVector yAxis(0.f, 1.f, 0.f);
    m_orthoDirection = yAxis.GetCrossProduct(m_axisDirection).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::FillLayerFitContributionMap(const CartesianPointVector &coordinateVector)
{
    if (m_layerPitch < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((m_axisDirection.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon()) ||
        (m_orthoDirection.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (!m_layerFitContributionMap.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (CartesianPointVector::const_iterator iter = coordinateVector.begin(), iterEnd = coordinateVector.end(); iter != iterEnd; ++iter)
    {
        float rL(0.f), rT(0.f);
        this->GetLocalPosition(*iter, rL, rT);
        m_layerFitContributionMap[this->GetLayer(rL)].AddPoint(rL, rT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::PerformSlidingLinearFit()
{
    if (!m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if ((m_layerPitch < std::numeric_limits<float>::epsilon()) || (m_layerFitContributionMap.empty()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    unsigned int slidingNPoints(0);
    double slidingSumT(0.), slidingSumL(0.), slidingSumTT(0.), slidingSumLT(0.), slidingSumLL(0.);

    const LayerFitContributionMap &layerFitContributionMap(this->GetLayerFitContributionMap());
    const int innerLayer(layerFitContributionMap.begin()->first);
    const int layerFitHalfWindow(static_cast<int>(this->GetLayerFitHalfWindow()));

    for (int iLayer = innerLayer; iLayer < innerLayer + layerFitHalfWindow; ++iLayer)
    {
        LayerFitContributionMap::const_iterator lyrIter = layerFitContributionMap.find(iLayer);

        if (layerFitContributionMap.end() != lyrIter)
        {
            slidingSumT += lyrIter->second.GetSumT();
            slidingSumL += lyrIter->second.GetSumL();
            slidingSumTT += lyrIter->second.GetSumTT();
            slidingSumLT += lyrIter->second.GetSumLT();
            slidingSumLL += lyrIter->second.GetSumLL();
            slidingNPoints += lyrIter->second.GetNPoints();
        }
    }

    const int outerLayer(layerFitContributionMap.rbegin()->first);

    for (int iLayer = innerLayer; iLayer <= outerLayer; ++iLayer)
    {
        const int fwdLayer(iLayer + layerFitHalfWindow);
        LayerFitContributionMap::const_iterator fwdIter = layerFitContributionMap.find(fwdLayer);

        if (layerFitContributionMap.end() != fwdIter)
        {
            slidingSumT += fwdIter->second.GetSumT();
            slidingSumL += fwdIter->second.GetSumL();
            slidingSumTT += fwdIter->second.GetSumTT();
            slidingSumLT += fwdIter->second.GetSumLT();
            slidingSumLL += fwdIter->second.GetSumLL();
            slidingNPoints += fwdIter->second.GetNPoints();
        }

        const int bwdLayer(iLayer - layerFitHalfWindow - 1);
        LayerFitContributionMap::const_iterator bwdIter = layerFitContributionMap.find(bwdLayer);

        if (layerFitContributionMap.end() != bwdIter)
        {
            slidingSumT -= bwdIter->second.GetSumT();
            slidingSumL -= bwdIter->second.GetSumL();
            slidingSumTT -= bwdIter->second.GetSumTT();
            slidingSumLT -= bwdIter->second.GetSumLT();
            slidingSumLL -= bwdIter->second.GetSumLL();
            slidingNPoints -= bwdIter->second.GetNPoints();
        }

        // require three points for meaningful results
        if (slidingNPoints <= 2)
            continue;

        // only fill the result map if there is an entry in the contribution map
        if (layerFitContributionMap.end() == layerFitContributionMap.find(iLayer))
            continue;

        const double denominator(slidingSumLL - slidingSumL * slidingSumL / static_cast<double>(slidingNPoints));

        if (std::fabs(denominator) < std::numeric_limits<float>::epsilon())
            continue;

        const double gradient((slidingSumLT - slidingSumL * slidingSumT / static_cast<double>(slidingNPoints)) / denominator);
        const double intercept((slidingSumLL * slidingSumT / static_cast<double>(slidingNPoints) - slidingSumL * slidingSumLT / static_cast<double>(slidingNPoints)) / denominator);
        double variance((slidingSumTT - 2. * intercept * slidingSumT - 2. * gradient * slidingSumLT + intercept * intercept * static_cast<double>(slidingNPoints) + 2. * gradient * intercept * slidingSumL + gradient * gradient * slidingSumLL) / (1. + gradient * gradient));

        if (variance < -std::numeric_limits<float>::epsilon())
            continue;

        if (variance < std::numeric_limits<float>::epsilon())
            variance = 0.;

        const double rms(std::sqrt(variance / static_cast<double>(slidingNPoints)));
        const double l(this->GetL(iLayer));
        const double fitT(intercept + gradient * l);

        const LayerFitResult layerFitResult(l, fitT, gradient, rms);
        (void) m_layerFitResultMap.insert(LayerFitResultMap::value_type(iLayer, layerFitResult));
    }

    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::FindSlidingFitSegments()
{
    if (!m_fitSegmentList.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    unsigned int nSustainedSteps(0);
    float sustainedDirectionStartX(0.f), sustainedDirectionEndX(0.f);

    CartesianVector previousPosition(0.f, 0.f, 0.f);
    TransverseDirection previousDirection(UNKNOWN), sustainedDirection(UNKNOWN);

    LayerFitResultMap::const_iterator sustainedDirectionStartIter, sustainedDirectionEndIter;
    const LayerFitResultMap &layerFitResultMap(this->GetLayerFitResultMap());

    for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        CartesianVector position(0.f, 0.f, 0.f);
        this->GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);

        // TODO currentDirection could also be UNCHANGED_IN_X
        const TransverseDirection currentDirection(((position - previousPosition).GetX() > 0.f) ? POSITIVE_IN_X : NEGATIVE_IN_X);

        if (previousDirection == currentDirection)
        {
            ++nSustainedSteps;

            if (nSustainedSteps > 2)
            {
                sustainedDirection = currentDirection;
                sustainedDirectionEndIter = iter;
                sustainedDirectionEndX = position.GetX();
            }
        }
        else
        {
            if ((POSITIVE_IN_X == sustainedDirection) || (NEGATIVE_IN_X == sustainedDirection))
                m_fitSegmentList.push_back(FitSegment(sustainedDirectionStartIter->first, sustainedDirectionEndIter->first, sustainedDirectionStartX, sustainedDirectionEndX));

            nSustainedSteps = 0;
            sustainedDirection = UNKNOWN;
            sustainedDirectionStartIter = iter;
            sustainedDirectionStartX = position.GetX();
        }

        previousPosition = position;
        previousDirection = currentDirection;
    }

    if ((POSITIVE_IN_X == sustainedDirection) || (NEGATIVE_IN_X == sustainedDirection))
        m_fitSegmentList.push_back(FitSegment(sustainedDirectionStartIter->first, sustainedDirectionEndIter->first, sustainedDirectionStartX, sustainedDirectionEndX));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetMinAndMaxCoordinate(const bool isX, float &min, float &max) const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    min = std::numeric_limits<float>::max();
    max = -std::numeric_limits<float>::max();

    for (LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin(), iterEnd = m_layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        CartesianVector globalPosition(0.f, 0.f, 0.f);
        this->GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), globalPosition);
        const float coordinate(isX ? globalPosition.GetX() : globalPosition.GetZ());
        min = std::min(min, coordinate);
        max = std::max(max, coordinate);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalFitPosition(const LayerInterpolation &layerInterpolation) const
{
    const LayerFitResultMap::const_iterator firstLayerIter(layerInterpolation.GetStartLayerIter());
    const LayerFitResultMap::const_iterator secondLayerIter(layerInterpolation.GetEndLayerIter());

    const float firstWeight(layerInterpolation.GetStartLayerWeight());
    const float secondWeight(layerInterpolation.GetEndLayerWeight());

    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    CartesianVector firstLayerPosition(0.f, 0.f, 0.f);
    this->GetGlobalPosition(firstLayerIter->second.GetL(), firstLayerIter->second.GetFitT(), firstLayerPosition);

    if (firstLayerIter == secondLayerIter)
        return firstLayerPosition;

    CartesianVector secondLayerPosition(0.f, 0.f, 0.f);
    this->GetGlobalPosition(secondLayerIter->second.GetL(), secondLayerIter->second.GetFitT(), secondLayerPosition);

    if (firstWeight + secondWeight < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return ((firstLayerPosition * firstWeight + secondLayerPosition * secondWeight) * (1.f / (firstWeight + secondWeight)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalFitDirection(const LayerInterpolation &layerInterpolation) const
{
    const LayerFitResultMap::const_iterator firstLayerIter(layerInterpolation.GetStartLayerIter());
    const LayerFitResultMap::const_iterator secondLayerIter(layerInterpolation.GetEndLayerIter());

    const float firstWeight(layerInterpolation.GetStartLayerWeight());
    const float secondWeight(layerInterpolation.GetEndLayerWeight());

    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    CartesianVector firstLayerDirection(0.f,0.f,0.f);
    this->GetGlobalDirection(firstLayerIter->second.GetGradient(),firstLayerDirection);

    if (firstLayerIter == secondLayerIter)
        return firstLayerDirection;

    CartesianVector secondLayerDirection(0.f,0.f,0.f);
    this->GetGlobalDirection(secondLayerIter->second.GetGradient(),secondLayerDirection);

    if (firstWeight + secondWeight < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return ((firstLayerDirection * firstWeight + secondLayerDirection * secondWeight).GetUnitVector());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetFitRms(const LayerInterpolation &layerInterpolation) const
{
    const LayerFitResultMap::const_iterator firstLayerIter(layerInterpolation.GetStartLayerIter());
    const LayerFitResultMap::const_iterator secondLayerIter(layerInterpolation.GetEndLayerIter());

    const float firstWeight(layerInterpolation.GetStartLayerWeight());
    const float secondWeight(layerInterpolation.GetEndLayerWeight());

    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float firstLayerRms(firstLayerIter->second.GetRms());

    if (firstLayerIter == secondLayerIter)
        return firstLayerRms;

    const float secondLayerRms(secondLayerIter->second.GetRms());

    if (firstWeight + secondWeight < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return ((firstLayerRms * firstWeight + secondLayerRms * secondWeight) / (firstWeight + secondWeight));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::LongitudinalInterpolation(const float rL, LayerInterpolation &layerInterpolation) const
{
    double firstWeight(0.), secondWeight(0.);
    LayerFitResultMap::const_iterator firstLayerIter, secondLayerIter;

    const StatusCode statusCode(this->GetLongitudinalSurroundingLayers(rL, firstLayerIter, secondLayerIter));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    this->GetLongitudinalInterpolationWeights(rL, firstLayerIter, secondLayerIter, firstWeight, secondWeight);
    layerInterpolation = LayerInterpolation(firstLayerIter, secondLayerIter, firstWeight, secondWeight);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::TransverseInterpolation(const float x, const FitSegment &fitSegment, LayerInterpolation &layerInterpolation) const
{
    double firstWeight(0.), secondWeight(0.);
    LayerFitResultMap::const_iterator firstLayerIter, secondLayerIter;

    const StatusCode statusCode(this->GetTransverseSurroundingLayers(x, fitSegment.GetStartLayer(), fitSegment.GetEndLayer(), firstLayerIter, secondLayerIter));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    this->GetTransverseInterpolationWeights(x, firstLayerIter, secondLayerIter, firstWeight, secondWeight);
    layerInterpolation = LayerInterpolation(firstLayerIter, secondLayerIter, firstWeight, secondWeight);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::TransverseInterpolation(const float x, LayerInterpolationList &layerInterpolationList) const
{
    for (FitSegmentList::const_iterator iter = m_fitSegmentList.begin(), iterEnd = m_fitSegmentList.end(); iter != iterEnd; ++iter)
    {
        LayerInterpolation layerInterpolation;
        const StatusCode statusCode(this->TransverseInterpolation(x, *iter, layerInterpolation));

        if (STATUS_CODE_SUCCESS == statusCode)
        {
            layerInterpolationList.push_back(layerInterpolation);
        }
        else if (STATUS_CODE_NOT_FOUND != statusCode)
        {
            return statusCode;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetLongitudinalSurroundingLayers(const float rL, LayerFitResultMap::const_iterator &firstLayerIter,
    LayerFitResultMap::const_iterator &secondLayerIter) const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // Get minimum, maximum and input layers
    const int minLayer(m_layerFitResultMap.begin()->first), maxLayer(m_layerFitResultMap.rbegin()->first);
    const int thisLayer(this->GetLayer(rL));

    // Allow special case of single-layer sliding fit result
    if (minLayer == thisLayer && thisLayer == maxLayer)
    {
        firstLayerIter = m_layerFitResultMap.find(minLayer);
        secondLayerIter = m_layerFitResultMap.find(maxLayer);
        return STATUS_CODE_SUCCESS;
    }

    // For multi-layer sliding fit result, set start layer and bail out if out of range
    const int startLayer((thisLayer >= maxLayer) ? thisLayer - 1 : thisLayer);

    if ((startLayer < minLayer) || (startLayer >= maxLayer))
        return STATUS_CODE_NOT_FOUND;

    // First layer iterator
    firstLayerIter = m_layerFitResultMap.end();

    for (int iLayer = startLayer; iLayer >= minLayer; --iLayer)
    {
        firstLayerIter = m_layerFitResultMap.find(iLayer);

        if (m_layerFitResultMap.end() != firstLayerIter)
            break;
    }

    if (m_layerFitResultMap.end() == firstLayerIter)
        return STATUS_CODE_NOT_FOUND;

    // Second layer iterator
    secondLayerIter = m_layerFitResultMap.end();

    for (int iLayer = startLayer + 1; iLayer <= maxLayer; ++iLayer)
    {
        secondLayerIter = m_layerFitResultMap.find(iLayer);

        if (m_layerFitResultMap.end() != secondLayerIter)
            break;
    }

    if (m_layerFitResultMap.end() == secondLayerIter)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitResult::GetTransverseSurroundingLayers(const float x, const int minLayer, const int maxLayer,
    LayerFitResultMap::const_iterator &firstLayerIter, LayerFitResultMap::const_iterator &secondLayerIter) const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_iterator minLayerIter = m_layerFitResultMap.find(minLayer);
    if (m_layerFitResultMap.end() == minLayerIter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    LayerFitResultMap::const_iterator maxLayerIter = m_layerFitResultMap.find(maxLayer);
    if (m_layerFitResultMap.end() == maxLayerIter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    CartesianVector minPosition(0.f, 0.f, 0.f), maxPosition(0.f, 0.f, 0.f);
    this->GetGlobalPosition(minLayerIter->second.GetL(), minLayerIter->second.GetFitT(), minPosition);
    this->GetGlobalPosition(maxLayerIter->second.GetL(), maxLayerIter->second.GetFitT(), maxPosition);

    if ((std::fabs(maxPosition.GetX() - minPosition.GetX()) < std::numeric_limits<float>::epsilon()))
        return STATUS_CODE_NOT_FOUND;

    // Find start layer
    const float minL(minLayerIter->second.GetL());
    const float maxL(maxLayerIter->second.GetL());
    const float startL(minL + (maxL - minL) * (x - minPosition.GetX()) / (maxPosition.GetX() - minPosition.GetX()));
    const int startLayer(std::max(minLayer, std::min(maxLayer, this->GetLayer(startL))));

    // Find nearest layer iterator to start layer
    LayerFitResultMap::const_iterator startLayerIter = m_layerFitResultMap.end();
    CartesianVector startLayerPosition(0.f, 0.f, 0.f);

    for (int iLayer = startLayer; iLayer <= maxLayer; ++iLayer)
    {
        startLayerIter = m_layerFitResultMap.find(iLayer);

        if (m_layerFitResultMap.end() != startLayerIter)
            break;
    }

    if (m_layerFitResultMap.end() == startLayerIter)
        return STATUS_CODE_NOT_FOUND;

    this->GetGlobalPosition(startLayerIter->second.GetL(), startLayerIter->second.GetFitT(), startLayerPosition);

    const bool startIsAhead((startLayerPosition.GetX() - x) > std::numeric_limits<float>::epsilon());
    const bool increasesWithLayers(maxPosition.GetX() > minPosition.GetX());
    const int increment = ((startIsAhead == increasesWithLayers) ? -1 : +1);

    // Find surrounding layer iterators
    // (Second layer iterator comes immediately after the fit has crossed the target X coordinate
    //  and first layer iterator comes immediately before the second layer iterator).
    firstLayerIter = m_layerFitResultMap.end();
    secondLayerIter = m_layerFitResultMap.end();

    CartesianVector firstLayerPosition(0.f, 0.f, 0.f);
    CartesianVector secondLayerPosition(0.f, 0.f, 0.f);

    for (int iLayer = startLayerIter->first; (iLayer >= minLayer) && (iLayer <= maxLayer); iLayer += increment)
    {
        LayerFitResultMap::const_iterator tempIter = m_layerFitResultMap.find(iLayer);
        if (m_layerFitResultMap.end() == tempIter)
            continue;

        firstLayerIter = secondLayerIter;
        firstLayerPosition = secondLayerPosition;
        secondLayerIter = tempIter;

        this->GetGlobalPosition(secondLayerIter->second.GetL(), secondLayerIter->second.GetFitT(), secondLayerPosition);
        const bool isAhead(secondLayerPosition.GetX() > x);

        if (startIsAhead != isAhead)
            break;

        firstLayerIter = m_layerFitResultMap.end();
    }

    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetLongitudinalInterpolationWeights(const float rL, const LayerFitResultMap::const_iterator &firstLayerIter,
    const LayerFitResultMap::const_iterator &secondLayerIter, double &firstWeight, double &secondWeight) const
{
    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const double deltaL(rL - firstLayerIter->second.GetL());
    const double deltaLLayers(secondLayerIter->second.GetL() - firstLayerIter->second.GetL());

    if (std::fabs(deltaLLayers) > std::numeric_limits<float>::epsilon())
    {
        firstWeight = 1. - deltaL / deltaLLayers;
        secondWeight = deltaL / deltaLLayers;
    }
    else
    {
        firstWeight = 0.5;
        secondWeight = 0.5;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetTransverseInterpolationWeights(const float x, const LayerFitResultMap::const_iterator &firstLayerIter,
    const LayerFitResultMap::const_iterator &secondLayerIter, double &firstWeight, double &secondWeight) const
{
    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector firstLayerPosition(0.f, 0.f, 0.f);
    CartesianVector secondLayerPosition(0.f, 0.f, 0.f);

    this->GetGlobalPosition(firstLayerIter->second.GetL(), firstLayerIter->second.GetFitT(), firstLayerPosition);
    this->GetGlobalPosition(secondLayerIter->second.GetL(), secondLayerIter->second.GetFitT(), secondLayerPosition);

    const double deltaP(x - firstLayerPosition.GetX());
    const double deltaPLayers(secondLayerPosition.GetX() - firstLayerPosition.GetX());

    if (std::fabs(deltaPLayers) > std::numeric_limits<float>::epsilon())
    {
        firstWeight = 1. - deltaP / deltaPLayers;
        secondWeight = deltaP / deltaPLayers;
    }
    else
    {
        firstWeight = 0.5;
        secondWeight = 0.5;
    }
}

} // namespace lar_content
