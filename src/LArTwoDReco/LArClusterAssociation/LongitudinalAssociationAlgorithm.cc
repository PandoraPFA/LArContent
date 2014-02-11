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

namespace lar
{

void LongitudinalAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < 4)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // ATTN: This method assumes that clusters have been sorted by layer
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

    if ((pOuterCluster->GetInnerPseudoLayer() < pInnerCluster->GetInnerPseudoLayer() + 3) ||
        (pInnerCluster->GetOuterPseudoLayer() + 3 > pOuterCluster->GetOuterPseudoLayer()))
    {
        return false;
    }

    if ((pInnerCluster->GetOuterPseudoLayer() > pOuterCluster->GetInnerPseudoLayer() + 1) ||
        (pOuterCluster->GetInnerPseudoLayer() > pInnerCluster->GetOuterPseudoLayer() + 7))
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

    if ((innerEndCentroid-outerStartCentroid).GetMagnitudeSquared() > 10.f)
        return false;

    CaloHit *pOuterLayerHit = *(pInnerCluster->GetOrderedCaloHitList().rbegin()->second->begin());
    const float hitSizeX(pOuterLayerHit->GetCellLengthScale());
    const float hitSizeZ(pOuterLayerHit->GetCellThickness());

    ClusterHelper::ClusterFitResult innerEndFit, outerStartFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pInnerCluster, 30, innerEndFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pOuterCluster, 30, outerStartFit));

    if (this->AreClustersAssociated(innerEndCentroid, outerStartCentroid, hitSizeX, hitSizeZ, innerEndFit, outerStartFit))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalAssociationAlgorithm::AreClustersAssociated(const CartesianVector &innerClusterEnd, const CartesianVector &outerClusterStart,
    const float hitSizeX, const float hitSizeZ, const ClusterHelper::ClusterFitResult &innerFit, const ClusterHelper::ClusterFitResult &outerFit) const
{
    if (!innerFit.IsFitSuccessful() || !outerFit.IsFitSuccessful())
        return false;

    if (innerFit.GetDirection().GetCosOpeningAngle(outerFit.GetDirection()) < 0.985f)
        return false;

    const CartesianVector innerEndFit1(innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(innerClusterEnd - innerFit.GetIntercept())));
    const CartesianVector innerEndFit2(outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(innerClusterEnd - outerFit.GetIntercept())));

    const CartesianVector outerStartFit1(outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(outerClusterStart - outerFit.GetIntercept())));
    const CartesianVector outerStartFit2(innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(outerClusterStart - innerFit.GetIntercept())));

    const CartesianVector clusterSeparation(outerClusterStart - innerClusterEnd);

    if ((std::fabs(clusterSeparation.GetX() / hitSizeX) < 2.f) && (std::fabs(clusterSeparation.GetZ() / hitSizeZ) < 2.f))
        return true;

    const CartesianVector fittedSeparation(outerStartFit1 - innerEndFit1);

    if ((std::fabs(fittedSeparation.GetX() / hitSizeX) < 2.f) && (std::fabs(fittedSeparation.GetZ() / hitSizeZ) < 2.f))
        return true;

    const CartesianVector fittedInnerSeparation(innerEndFit2 - innerEndFit1);

    if ((std::fabs(fittedInnerSeparation.GetX() / hitSizeX) < 2.f) && (std::fabs(fittedInnerSeparation.GetZ() / hitSizeZ) < 2.f))
        return true;

    const CartesianVector fittedOuterSeparation(outerStartFit2 - outerStartFit1);

    if ((std::fabs(fittedOuterSeparation.GetX() / hitSizeX) < 2.f) && (std::fabs(fittedOuterSeparation.GetZ() / hitSizeZ) < 2.f))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
