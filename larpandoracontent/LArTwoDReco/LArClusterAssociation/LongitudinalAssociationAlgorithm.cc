/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the longitudinal association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"

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
    m_maxLongitudinalDisplacement(2.f),
    m_hitSizeZ(0.3f),
    m_hitSizeX(0.5f),
    m_view(TPC_3D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    if (!pClusterList->empty())
        m_view = LArClusterHelper::GetClusterHitType(pClusterList->front());

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

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
        const Cluster *const pInnerCluster = *iterI;

        for (ClusterVector::const_iterator iterJ = iterI, iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            const Cluster *const pOuterCluster = *iterJ;

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

bool LongitudinalAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
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

    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), m_view)};
    const float maxGapDistanceSquaredAdjusted{ratio * ratio * m_maxGapDistanceSquared};

    if ((innerEndCentroid - outerStartCentroid).GetMagnitudeSquared() > maxGapDistanceSquaredAdjusted)
        return false;

    ClusterFitResult innerEndFit, outerStartFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterFitHelper::FitEnd(pInnerCluster, m_fitLayers, innerEndFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterFitHelper::FitStart(pOuterCluster, m_fitLayers, outerStartFit));

    if (this->AreClustersAssociated(innerEndCentroid, outerStartCentroid, innerEndFit, outerStartFit))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalAssociationAlgorithm::AreClustersAssociated(const CartesianVector &innerClusterEnd,
    const CartesianVector &outerClusterStart, const ClusterFitResult &innerFit, const ClusterFitResult &outerFit) const
{
    if (!innerFit.IsFitSuccessful() || !outerFit.IsFitSuccessful())
        return false;

    if (innerFit.GetDirection().GetCosOpeningAngle(outerFit.GetDirection()) < m_minCosRelativeAngle)
        return false;

    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), m_view)};
    const float maxTransverseDisplacementAdjusted{ratio * m_maxTransverseDisplacement};
    const float maxLongitudinalDisplacementAdjusted{ratio * m_maxLongitudinalDisplacement};

    const CartesianVector innerEndFit1(
        innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(innerClusterEnd - innerFit.GetIntercept())));
    const CartesianVector innerEndFit2(
        outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(innerClusterEnd - outerFit.GetIntercept())));

    const CartesianVector outerStartFit1(
        outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(outerClusterStart - outerFit.GetIntercept())));
    const CartesianVector outerStartFit2(
        innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(outerClusterStart - innerFit.GetIntercept())));

    const CartesianVector clusterSeparation(outerClusterStart - innerClusterEnd);

    if ((std::fabs(clusterSeparation.GetX()) < m_hitSizeX * maxTransverseDisplacementAdjusted) &&
        (std::fabs(clusterSeparation.GetZ()) < m_hitSizeZ * maxLongitudinalDisplacementAdjusted))
        return true;

    const CartesianVector fittedSeparation(outerStartFit1 - innerEndFit1);

    if ((std::fabs(fittedSeparation.GetX()) < m_hitSizeX * maxTransverseDisplacementAdjusted) &&
        (std::fabs(fittedSeparation.GetZ()) < m_hitSizeZ * maxLongitudinalDisplacementAdjusted))
        return true;

    const CartesianVector fittedInnerSeparation(innerEndFit2 - innerEndFit1);

    if ((std::fabs(fittedInnerSeparation.GetX()) < m_hitSizeX * maxTransverseDisplacementAdjusted) &&
        (std::fabs(fittedInnerSeparation.GetZ()) < m_hitSizeZ * maxLongitudinalDisplacementAdjusted))
        return true;

    const CartesianVector fittedOuterSeparation(outerStartFit2 - outerStartFit1);

    if ((std::fabs(fittedOuterSeparation.GetX()) < m_hitSizeX * maxTransverseDisplacementAdjusted) &&
        (std::fabs(fittedOuterSeparation.GetZ()) < m_hitSizeZ * maxLongitudinalDisplacementAdjusted))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapLayers", m_maxGapLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FitLayers", m_fitLayers));

    float maxGapDistance = std::sqrt(m_maxGapDistanceSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapDistance", maxGapDistance));
    m_maxGapDistanceSquared = maxGapDistance * maxGapDistance;

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitSizeZ", m_hitSizeZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitSizeX", m_hitSizeX));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
