/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the longitudinal association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

LongitudinalAssociationAlgorithm::LongitudinalAssociationAlgorithm() :
    m_minClusterLayers(4),
    m_maxGapLayers(7),
    m_fitLayers(30),
    m_maxGapDistanceSquared(10.f),
    m_minCosRelativeAngle(0.985f),
    m_maxTransverseDisplacement(2.f),
    m_maxLongitudinalDisplacement(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < m_minClusterLayers)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // ATTN This method assumes that clusters have been sorted by layer
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterIEnd = clusterVector.end(); iterI != iterIEnd; ++iterI)
    {
        Cluster *pInnerCluster = *iterI;

        for (ClusterVector::const_iterator iterJ = iterI, iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            Cluster *pOuterCluster = *iterJ;

            if (pInnerCluster == pOuterCluster)
                continue;

            if (!this->AreClustersAssociated(pInnerCluster, pOuterCluster))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    const float currentEnergy(pCurrentCluster->GetHadronicEnergy());
    const float testEnergy(pTestCluster->GetHadronicEnergy());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && (testEnergy > currentEnergy))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && (testEnergy > currentEnergy))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalAssociationAlgorithm::AreClustersAssociated(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
{
    if (pOuterCluster->GetInnerPseudoLayer() < pInnerCluster->GetInnerPseudoLayer())
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // TODO Remove hardcoded numbers
    if ((pOuterCluster->GetInnerPseudoLayer() < pInnerCluster->GetInnerPseudoLayer() + 3) ||
        (pInnerCluster->GetOuterPseudoLayer() + 3 > pOuterCluster->GetOuterPseudoLayer()))
    {
        return false;
    }

    if ((pInnerCluster->GetOuterPseudoLayer() > pOuterCluster->GetInnerPseudoLayer() + 1) ||
        (pOuterCluster->GetInnerPseudoLayer() > pInnerCluster->GetOuterPseudoLayer() + m_maxGapLayers))
    {
        return false;
    }

    if ((2 * pInnerCluster->GetOuterPseudoLayer() < pOuterCluster->GetInnerPseudoLayer() + pInnerCluster->GetInnerPseudoLayer()) ||
        (pInnerCluster->GetOuterPseudoLayer() + pOuterCluster->GetOuterPseudoLayer() < 2 * pOuterCluster->GetInnerPseudoLayer()))
    {
        return false;
    }

    const CartesianVector innerEndCentroid(pInnerCluster->GetCentroid(pInnerCluster->GetOuterPseudoLayer()));
    const CartesianVector outerStartCentroid(pOuterCluster->GetCentroid(pOuterCluster->GetInnerPseudoLayer()));

    if ((innerEndCentroid-outerStartCentroid).GetMagnitudeSquared() > m_maxGapDistanceSquared)
        return false;

    CaloHit *pOuterLayerHit = *(pInnerCluster->GetOrderedCaloHitList().rbegin()->second->begin());
    const float hitSizeX(pOuterLayerHit->GetCellLengthScale());
    const float hitSizeZ(pOuterLayerHit->GetCellThickness());

    ClusterFitResult innerEndFit, outerStartFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterFitHelper::FitEnd(pInnerCluster, m_fitLayers, innerEndFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterFitHelper::FitStart(pOuterCluster, m_fitLayers, outerStartFit));

    if (this->AreClustersAssociated(innerEndCentroid, outerStartCentroid, hitSizeX, hitSizeZ, innerEndFit, outerStartFit))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalAssociationAlgorithm::AreClustersAssociated(const CartesianVector &innerClusterEnd, const CartesianVector &outerClusterStart,
    const float hitSizeX, const float hitSizeZ, const ClusterFitResult &innerFit, const ClusterFitResult &outerFit) const
{
    if (!innerFit.IsFitSuccessful() || !outerFit.IsFitSuccessful())
        return false;

    if (innerFit.GetDirection().GetCosOpeningAngle(outerFit.GetDirection()) < m_minCosRelativeAngle)
        return false;

    const CartesianVector innerEndFit1(innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(innerClusterEnd - innerFit.GetIntercept())));
    const CartesianVector innerEndFit2(outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(innerClusterEnd - outerFit.GetIntercept())));

    const CartesianVector outerStartFit1(outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(outerClusterStart - outerFit.GetIntercept())));
    const CartesianVector outerStartFit2(innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(outerClusterStart - innerFit.GetIntercept())));

    const CartesianVector clusterSeparation(outerClusterStart - innerClusterEnd);

    if ((std::fabs(clusterSeparation.GetX() / hitSizeX) < m_maxTransverseDisplacement) &&
        (std::fabs(clusterSeparation.GetZ() / hitSizeZ) < m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedSeparation(outerStartFit1 - innerEndFit1);

    if ((std::fabs(fittedSeparation.GetX() / hitSizeX) < m_maxTransverseDisplacement) &&
        (std::fabs(fittedSeparation.GetZ() / hitSizeZ) < m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedInnerSeparation(innerEndFit2 - innerEndFit1);

    if ((std::fabs(fittedInnerSeparation.GetX() / hitSizeX) < m_maxTransverseDisplacement) &&
        (std::fabs(fittedInnerSeparation.GetZ() / hitSizeZ) < m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedOuterSeparation(outerStartFit2 - outerStartFit1);

    if ((std::fabs(fittedOuterSeparation.GetX() / hitSizeX) < m_maxTransverseDisplacement) &&
        (std::fabs(fittedOuterSeparation.GetZ() / hitSizeZ) < m_maxLongitudinalDisplacement))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGapLayers", m_maxGapLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FitLayers", m_fitLayers));

    float maxGapDistance = std::sqrt(m_maxGapDistanceSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGapDistance", maxGapDistance));
    m_maxGapDistanceSquared = maxGapDistance * maxGapDistance;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
