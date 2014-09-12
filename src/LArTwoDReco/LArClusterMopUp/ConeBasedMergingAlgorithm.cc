/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the cone based merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArPlugins/LArTransformationPlugin.h"

#include "LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void ConeBasedMergingAlgorithm::ClusterMopUp(const ClusterList &pfoClusters, const ClusterList &remnantClusters,
    const ClusterToListNameMap &clusterToListNameMap) const
{
    ClusterAssociationMap clusterAssociationMap;
    const float slidingFitPitch(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());

    for (ClusterList::const_iterator pIter = pfoClusters.begin(), pIterEnd = pfoClusters.end(); pIter != pIterEnd; ++pIter)
    {
        try
        {
            Cluster *pPfoCluster(*pIter);
            const ConeParameters coneParameters(pPfoCluster, m_slidingFitWindow, slidingFitPitch, m_coneAngleCentile);

            for (ClusterList::const_iterator rIter = remnantClusters.begin(), rIterEnd = remnantClusters.end(); rIter != rIterEnd; ++rIter)
            {
                Cluster *pRemnantCluster(*rIter);
                const float boundedFraction(this->GetBoundedFraction(pRemnantCluster, coneParameters));

                if (boundedFraction < m_minBoundedFraction)
                    continue;

                AssociationDetails &associationDetails(clusterAssociationMap[pRemnantCluster]);

                if (!associationDetails.insert(AssociationDetails::value_type(pPfoCluster, boundedFraction)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }

    this->MakeClusterMerges(clusterAssociationMap, clusterToListNameMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeBasedMergingAlgorithm::GetBoundedFraction(const Cluster *const pCluster, const ConeParameters &coneParameters) const
{
    unsigned int nMatchedHits(0);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;
            const CartesianVector &positionVector(pCaloHit->GetPositionVector());

            if (coneParameters.GetPositionCosHalfAngle(positionVector) < coneParameters.GetConeCosHalfAngle())
                continue;

            if (coneParameters.GetConeAxisProjection(positionVector) > m_maxConeLengthMultiplier * coneParameters.GetConeLength())
                continue;

            ++nMatchedHits;
        }
    }

    return (static_cast<float>(nMatchedHits) / static_cast<float>(pCluster->GetNCaloHits()));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeBasedMergingAlgorithm::ConeParameters::ConeParameters(Cluster *pCluster, const unsigned int slidingFitWindow, const float slidingFitPitch,
        const float coneAngleCentile) :
    m_pCluster(pCluster),
    m_direction(0.f, 0.f, 0.f),
    m_apex(0.f, 0.f, 0.f),
    m_baseCentre(0.f, 0.f, 0.f),
    m_coneLength(0.f),
    m_coneCosHalfAngle(0.f),
    m_isForward(true)
{
    const TwoDSlidingFitResult fitResult(pCluster, slidingFitWindow, slidingFitPitch);
    this->GetDirectionEstimate(fitResult, m_isForward, m_direction);
    const CartesianVector basePosition(m_isForward ? fitResult.GetGlobalMaxLayerPosition() : fitResult.GetGlobalMinLayerPosition());

    m_apex = m_isForward ? fitResult.GetGlobalMinLayerPosition() : fitResult.GetGlobalMaxLayerPosition();
    m_baseCentre = m_apex + m_direction * (basePosition - m_apex).GetDotProduct(m_direction);
    m_coneLength = (m_baseCentre - m_apex).GetMagnitude();
    m_coneCosHalfAngle = this->GetCosHalfAngleEstimate(pCluster, m_direction, m_apex, coneAngleCentile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeBasedMergingAlgorithm::ConeParameters::GetDirectionEstimate(const TwoDSlidingFitResult &fitResult, bool &isForward, CartesianVector &direction) const
{
    HitsPerLayerMap hitsPerLayerMap;
    const OrderedCaloHitList &orderedCaloHitList(fitResult.GetCluster()->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            float rL(0.f), rT(0.f);
            fitResult.GetLocalPosition((*hitIter)->GetPositionVector(), rL, rT);
            hitsPerLayerMap[fitResult.GetLayer(rL)] += 1;
        }
    }

    CartesianVector directionSum(0.f, 0.f, 0.f);
    int layerCounter(0), hitCounter(0), layerSum(0), layerHitWeightedSum(0);
    const LayerFitResultMap &layerFitResultMap(fitResult.GetLayerFitResultMap());

    for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const int layer(iter->first);
            const unsigned int hitsInLayer(hitsPerLayerMap[layer]);

            CartesianVector fitDirection(0.f, 0.f, 0.f);
            fitResult.GetGlobalDirection(iter->second.GetGradient(), fitDirection);

            if (0 < hitsInLayer)
            {
                layerCounter += 1;
                hitCounter += hitsInLayer;

                layerSum += layer;
                layerHitWeightedSum += layer * hitsInLayer;
                directionSum += fitDirection * static_cast<float>(hitsInLayer);
            }
        }
        catch (StatusCodeException &) {}
    }

    if ((0 == hitCounter) || (0 == layerCounter))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float meanLayer(static_cast<float>(layerSum) / static_cast<float>(layerCounter));
    const float meanHitWeightedLayer(static_cast<float>(layerHitWeightedSum) / static_cast<float>(hitCounter));
    isForward = (meanHitWeightedLayer > meanLayer);

    const CartesianVector meanDirection((directionSum * (1.f / static_cast<float>(hitCounter))).GetUnitVector());
    direction = isForward ? meanDirection : meanDirection * -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeBasedMergingAlgorithm::ConeParameters::GetCosHalfAngleEstimate(const Cluster *const pCluster, const CartesianVector &direction,
    const CartesianVector &apex, const float coneAngleCentile) const
{
    FloatVector cosHalfAngleValues;
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
            cosHalfAngleValues.push_back(direction.GetCosOpeningAngle((*hitIter)->GetPositionVector() - apex));
    }

    std::sort(cosHalfAngleValues.begin(), cosHalfAngleValues.end());

    if (cosHalfAngleValues.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const unsigned int cosHalfAngleBin(coneAngleCentile * cosHalfAngleValues.size());
    return cosHalfAngleValues.at(cosHalfAngleBin);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeBasedMergingAlgorithm::ConeParameters::GetPositionCosHalfAngle(const CartesianVector &position) const
{
    return (m_direction.GetCosOpeningAngle(position - m_apex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeBasedMergingAlgorithm::ConeParameters::GetConeAxisProjection(const CartesianVector &position) const
{
    return (m_direction.GetDotProduct(position - m_apex));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeBasedMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "SlidingFitWindow", m_slidingFitWindow));

    m_coneAngleCentile = 0.25f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ConeAngleCentile", m_coneAngleCentile));

    m_maxConeLengthMultiplier = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxConeLengthMultiplier", m_maxConeLengthMultiplier));

    m_minBoundedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinBoundedFraction", m_minBoundedFraction));

    return ClusterMopUpAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
