/**
 *  @file   larpandoracontent/LArHelpers/LArHitWidthHelper.cc
 *
 *  @brief  Implementation of the lar hit width helper class.
 *
 *  $Log: $
 */
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

LArHitWidthHelper::ConstituentHit::ConstituentHit(const CartesianVector &positionVector, const float hitWidth, const Cluster *const pParentClusterAddress) :
    m_positionVector(positionVector),
    m_hitWidth(hitWidth),
    m_pParentClusterAddress(pParentClusterAddress)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHitWidthHelper::ConstituentHit::SortByDistanceToPoint::operator() (const ConstituentHit &lhs, const ConstituentHit &rhs)
{
    const CartesianVector &lhsPosition(lhs.GetPositionVector());
    const CartesianVector &rhsPosition(rhs.GetPositionVector());

    return (m_referencePoint.GetDistanceSquared(lhsPosition) < m_referencePoint.GetDistanceSquared(rhsPosition));   
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
    
LArHitWidthHelper::ClusterParameters::ClusterParameters(const Cluster *const pCluster, const float maxConstituentHitWidth, const bool isUniformHits,
        const float hitWidthScalingFactor) :
    m_pCluster(pCluster),
    m_numCaloHits(pCluster->GetNCaloHits()), 
    m_constituentHitVector(LArHitWidthHelper::GetConstituentHits(pCluster, maxConstituentHitWidth, hitWidthScalingFactor, isUniformHits)),
    m_totalWeight(LArHitWidthHelper::GetTotalClusterWeight(m_constituentHitVector)),
    m_lowerXExtrema(LArHitWidthHelper::GetExtremalCoordinatesLowerX(m_constituentHitVector)), 
    m_higherXExtrema(LArHitWidthHelper::GetExtremalCoordinatesHigherX(m_constituentHitVector))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
LArHitWidthHelper::ClusterParameters::ClusterParameters(const Cluster *const pCluster, const LArHitWidthHelper::ConstituentHitVector &constituentHitVector) :
    m_pCluster(pCluster),
    m_numCaloHits(pCluster->GetNCaloHits()), 
    m_constituentHitVector(constituentHitVector),
    m_totalWeight(LArHitWidthHelper::GetTotalClusterWeight(m_constituentHitVector)),
    m_lowerXExtrema(LArHitWidthHelper::GetExtremalCoordinatesLowerX(m_constituentHitVector)), 
    m_higherXExtrema(LArHitWidthHelper::GetExtremalCoordinatesHigherX(m_constituentHitVector))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHitWidthHelper::ClusterParameters::ClusterParameters(const Cluster *const pCluster, const unsigned int numCaloHits, const float totalWeight,
        const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, const CartesianVector &lowerXExtrema, const CartesianVector &higherXExtrema) :
    m_pCluster(pCluster),
    m_numCaloHits(numCaloHits), 
    m_constituentHitVector(constituentHitVector),
    m_totalWeight(totalWeight),
    m_lowerXExtrema(lowerXExtrema), 
    m_higherXExtrema(higherXExtrema)
{
}
  
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHitWidthHelper::SortByHigherXExtrema::SortByHigherXExtrema(const ClusterToParametersMap &clusterToParametersMap) :
    m_clusterToParametersMap(clusterToParametersMap)
{
    if (m_clusterToParametersMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHitWidthHelper::SortByHigherXExtrema::operator() (const Cluster *const pLhs, const Cluster *const pRhs)
{
    const LArHitWidthHelper::ClusterParameters& lhsClusterParameters(LArHitWidthHelper::GetClusterParameters(pLhs, m_clusterToParametersMap));
    const LArHitWidthHelper::ClusterParameters& rhsClusterParameters(LArHitWidthHelper::GetClusterParameters(pRhs, m_clusterToParametersMap));

    return (lhsClusterParameters.GetHigherXExtrema().GetX() < rhsClusterParameters.GetHigherXExtrema().GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

const LArHitWidthHelper::ClusterParameters& LArHitWidthHelper::GetClusterParameters(const Cluster *const pCluster, const ClusterToParametersMap &clusterToParametersMap)
{
    if (clusterToParametersMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    const auto clusterParametersIter(clusterToParametersMap.find(pCluster));

    if (clusterParametersIter == clusterToParametersMap.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return clusterParametersIter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHitWidthHelper::ConstituentHitVector LArHitWidthHelper::GetConstituentHits(const Cluster *const pCluster, const float maxConstituentHitWidth,
    const float hitWidthScalingFactor, const bool isUniform)
{
    if (maxConstituentHitWidth < std::numeric_limits<float>::epsilon())
    {
        std::cout << "LArHitWidthHelper::GetConstituentHits - Negative or equivalent to zero constitent hit width not allowed" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }
    
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    ConstituentHitVector constituentHitVector;
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const float hitWidth = pCaloHit->GetCellSize1() * hitWidthScalingFactor;
            const unsigned int numberOfConstituentHits = std::ceil(hitWidth / maxConstituentHitWidth);
            if (isUniform)
            {
                LArHitWidthHelper::SplitHitIntoConstituents(pCaloHit, pCluster, numberOfConstituentHits, maxConstituentHitWidth, constituentHitVector);
            }
            else
            {
                const float constituentHitWidth = hitWidth / numberOfConstituentHits;
                LArHitWidthHelper::SplitHitIntoConstituents(pCaloHit, pCluster, numberOfConstituentHits, constituentHitWidth, constituentHitVector);
            }
        }
    }
    
    return constituentHitVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHitWidthHelper::SplitHitIntoConstituents(const CaloHit *const pCaloHit, const Cluster *const pCluster, const unsigned int numberOfConstituentHits,
    const float constituentHitWidth, LArHitWidthHelper::ConstituentHitVector &constituentHitVector)
{
    const CartesianVector &hitCenter(pCaloHit->GetPositionVector());
    const bool isOdd(numberOfConstituentHits % 2 == 1);
    float xDistanceFromCenter(0.f);

    // find constituent hit centers by moving out from the original hit center position
    unsigned int loopIterations(std::ceil(numberOfConstituentHits / 2.0));
    for (unsigned int i = 0; i < loopIterations; ++i)
    {
        if (i == 0)
        {
            if (isOdd)
            {
                constituentHitVector.push_back(ConstituentHit(hitCenter, constituentHitWidth, pCluster));
                continue;
            }
            else
            {
                xDistanceFromCenter += constituentHitWidth / 2;
            }
        }
        else
        {
            xDistanceFromCenter += constituentHitWidth;
        }

        CartesianVector positivePosition(hitCenter + CartesianVector(xDistanceFromCenter, 0.f, 0.f)), negativePosition(hitCenter - CartesianVector(xDistanceFromCenter, 0.f, 0.f));
        
        constituentHitVector.push_back(ConstituentHit(positivePosition, constituentHitWidth, pCluster));
        constituentHitVector.push_back(ConstituentHit(negativePosition, constituentHitWidth, pCluster));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHitWidthHelper::GetTotalClusterWeight(const ConstituentHitVector &constituentHitVector)
{
    float clusterWeight(0.f); 
    for (const ConstituentHit &constituentHit : constituentHitVector)
        clusterWeight += constituentHit.GetHitWidth();
    
    return clusterWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHitWidthHelper::GetOriginalTotalClusterWeight(const Cluster *const pCluster)
{
    float clusterWeight(0.f);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
            clusterWeight += pCaloHit->GetCellSize1();
    }

    return clusterWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianPointVector LArHitWidthHelper::GetConstituentHitPositionVector(const ConstituentHitVector &constituentHitVector)
{
    CartesianPointVector constituentHitPositionVector;

    for (const ConstituentHit &constituentHit : constituentHitVector)
        constituentHitPositionVector.push_back(constituentHit.GetPositionVector());

    return constituentHitPositionVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArHitWidthHelper::GetExtremalCoordinatesLowerX(const ConstituentHitVector &constituentHitVector) 
{
    CartesianVector lowerXCoordinate(0.f, 0.f, 0.f), higherXCoordinate(0.f, 0.f, 0.f);
    LArHitWidthHelper::GetExtremalCoordinatesX(constituentHitVector, lowerXCoordinate, higherXCoordinate);

    return lowerXCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArHitWidthHelper::GetExtremalCoordinatesHigherX(const ConstituentHitVector &constituentHitVector) 
{
    CartesianVector lowerXCoordinate(0.f, 0.f, 0.f), higherXCoordinate(0.f, 0.f, 0.f);
    LArHitWidthHelper::GetExtremalCoordinatesX(constituentHitVector, lowerXCoordinate, higherXCoordinate);

    return higherXCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHitWidthHelper::GetExtremalCoordinatesX(const ConstituentHitVector &constituentHitVector, CartesianVector &lowerXCoordinate,
    CartesianVector &higherXCoordinate)
{
    const CartesianPointVector &constituentHitPositionVector(GetConstituentHitPositionVector(constituentHitVector));

    CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(constituentHitPositionVector, innerCoordinate, outerCoordinate);

    // set the lower/higher XCoordinate (in the event of a tie, use z)
    const float deltaX(outerCoordinate.GetX() - innerCoordinate.GetX());
    const float deltaZ(outerCoordinate.GetZ() - innerCoordinate.GetZ());

    if ((deltaX > 0.f) || ((std::fabs(deltaX) < std::numeric_limits<float>::epsilon()) && (deltaZ > 0.f)))
    {
        lowerXCoordinate = innerCoordinate;
        higherXCoordinate = outerCoordinate;
    }
    else
    {
        lowerXCoordinate = outerCoordinate;
        higherXCoordinate = innerCoordinate;
    }
}

} // namespace lar_content
