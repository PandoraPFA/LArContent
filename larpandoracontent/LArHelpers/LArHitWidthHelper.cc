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

LArHitWidthHelper::ConstituentHit::ConstituentHit(const CartesianVector &positionVector, const float hitWidth, const Cluster *const parentClusterAddress) :
    m_positionVector(positionVector),
    m_hitWidth(hitWidth),
    m_parentClusterAddress(parentClusterAddress)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster* LArHitWidthHelper::ConstituentHit::GetParentClusterAddress() const
{
    if (!m_parentClusterAddress)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_parentClusterAddress;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHitWidthHelper::ConstituentHit::SortByDistanceToPoint::operator() (const ConstituentHit &lhs, const ConstituentHit &rhs)
{
    CartesianVector lhsPosition(lhs.GetPositionVector());
    CartesianVector rhsPosition(rhs.GetPositionVector());

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

const pandora::Cluster* LArHitWidthHelper::ClusterParameters::GetClusterAddress() const
{
    if (!m_pCluster)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
  
    return m_pCluster;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHitWidthHelper::SortByHigherXExtrema2::SortByHigherXExtrema2(const ClusterToParametersMap &clusterToParametersMap) :
    m_clusterToParametersMap(clusterToParametersMap)
{
    if (m_clusterToParametersMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHitWidthHelper::SortByHigherXExtrema2::operator() (const Cluster *const pLhs, const Cluster *const pRhs)
{
    const LArHitWidthHelper::ClusterParameters& lhsClusterParameters(LArHitWidthHelper::GetClusterParameters(pLhs, m_clusterToParametersMap));
    const LArHitWidthHelper::ClusterParameters& rhsClusterParameters(LArHitWidthHelper::GetClusterParameters(pRhs, m_clusterToParametersMap));

    return (lhsClusterParameters.GetHigherXExtrema().GetX() < rhsClusterParameters.GetHigherXExtrema().GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
  
LArHitWidthHelper::ClusterToParametersMapStore* LArHitWidthHelper::ClusterToParametersMapStore::m_instance{nullptr};

//------------------------------------------------------------------------------------------------------------------------------------------

LArHitWidthHelper::ClusterToParametersMapStore* LArHitWidthHelper::ClusterToParametersMapStore::Instance()
{
    if (!m_instance)
        m_instance = new ClusterToParametersMapStore();
    
    return m_instance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHitWidthHelper::ClusterToParametersMap& LArHitWidthHelper::ClusterToParametersMapStore::GetMap()
{
    return m_clusterToParametersMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

const LArHitWidthHelper::ClusterParameters& LArHitWidthHelper::GetClusterParameters(const Cluster *const pCluster)
{
    LArHitWidthHelper::ClusterToParametersMapStore* pClusterToParametersMapStore(LArHitWidthHelper::ClusterToParametersMapStore::Instance());
    const LArHitWidthHelper::ClusterToParametersMap &clusterToParametersMap(pClusterToParametersMapStore->GetMap());

    if (clusterToParametersMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    const auto clusterParametersIter(clusterToParametersMap.find(pCluster));

    if (clusterParametersIter == clusterToParametersMap.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return clusterParametersIter->second;
}

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
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    ConstituentHitVector constituentHitVector;
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const float hitWidth = pCaloHit->GetCellSize1() * hitWidthScalingFactor;
            const unsigned int numberOfConstituentHits = ceil(hitWidth / maxConstituentHitWidth);
            if (isUniform)
            {
                SplitHitIntoConstituents(pCaloHit, pCluster, numberOfConstituentHits, maxConstituentHitWidth, constituentHitVector);
            }
            else
            {
                const float constituentHitWidth = hitWidth / numberOfConstituentHits;
                SplitHitIntoConstituents(pCaloHit, pCluster, numberOfConstituentHits, constituentHitWidth, constituentHitVector);
            }
        }
    }
    
    return constituentHitVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHitWidthHelper::SplitHitIntoConstituents(const CaloHit *const pCaloHit, const Cluster *const pCluster, const unsigned int numberOfConstituentHits, const float constituentHitWidth,
        LArHitWidthHelper::ConstituentHitVector &constituentHitVector)
{
    // begin at the central hit position
    CartesianVector positionAlongHit(pCaloHit->GetPositionVector());
    if (numberOfConstituentHits % 2 == 1)
    {
        // if odd number of constituents, include the central position
        constituentHitVector.push_back(ConstituentHit(positionAlongHit, constituentHitWidth, pCluster));

        // find remaining positions either side of central position
	    for (unsigned int i = 0; i < ceil((numberOfConstituentHits - 1) / 2); ++i)
	    {
            positionAlongHit += CartesianVector(constituentHitWidth, 0.f, 0.f);
            constituentHitVector.push_back(ConstituentHit(positionAlongHit, constituentHitWidth, pCluster));
        }

        // reset position to centre
        positionAlongHit = pCaloHit->GetPositionVector();
	    for (unsigned int i = 0; i < ceil((numberOfConstituentHits - 1) / 2); ++i)
        {
            positionAlongHit -= CartesianVector(constituentHitWidth, 0.f, 0.f);
            constituentHitVector.push_back(ConstituentHit(positionAlongHit, constituentHitWidth, pCluster));
        }
    }
    else
    {
        // if even number of constituents, begin with a numberOfConstituents / 2 x offset from central position
        for (unsigned int i = 0; i < ceil(numberOfConstituentHits / 2); ++i)
        {
            positionAlongHit += (i == 0) ? CartesianVector(constituentHitWidth / 2, 0.f, 0.f) : CartesianVector(constituentHitWidth, 0.f, 0.f);                
            constituentHitVector.push_back(ConstituentHit(positionAlongHit, constituentHitWidth, pCluster));
        }

        positionAlongHit = pCaloHit->GetPositionVector();
	    for (unsigned int i =0; i < ceil(numberOfConstituentHits / 2); ++i)
	    {
            positionAlongHit -= (i == 0) ? CartesianVector(constituentHitWidth / 2, 0.f, 0.f) : CartesianVector(constituentHitWidth, 0.f, 0.f);
            constituentHitVector.push_back(ConstituentHit(positionAlongHit, constituentHitWidth, pCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHitWidthHelper::GetTotalClusterWeight(const ConstituentHitVector &constituentHitVector)
{
    float hitWeight(0.f); 
    for (const ConstituentHit &constituentHit : constituentHitVector)
        hitWeight += constituentHit.GetHitWidth();
    
    return hitWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHitWidthHelper::GetOriginalTotalClusterWeight(const Cluster *const pCluster)
{
    float hitWeight(0.f);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
            hitWeight += pCaloHit->GetCellSize1();
    }

    return hitWeight;
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

bool LArHitWidthHelper::SortByHigherXExtrema(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const LArHitWidthHelper::ClusterParameters lhsClusterParameters(LArHitWidthHelper::GetClusterParameters(pLhs));
    const LArHitWidthHelper::ClusterParameters rhsClusterParameters(LArHitWidthHelper::GetClusterParameters(pRhs));

    return (lhsClusterParameters.GetHigherXExtrema().GetX() < rhsClusterParameters.GetHigherXExtrema().GetX());
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

void LArHitWidthHelper::GetExtremalCoordinatesX(const ConstituentHitVector &constituentHitVector, CartesianVector &lowerXCoordinate, CartesianVector &higherXCoordinate)
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
