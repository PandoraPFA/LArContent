/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/ThreeDAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the 3D association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ThreeDAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeDAssociationAlgorithm::ThreeDAssociationAlgorithm() :
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

void ThreeDAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    if (!pClusterList->empty())
        m_view = LArClusterHelper::GetClusterHitType(pClusterList->front());

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        float xmin, xmax;
        pCluster->GetClusterSpanX(xmin, xmax);

        float zmin, zmax;
        pCluster->GetClusterSpanZ(xmin, xmax, zmin, zmax);

        float ymin = +std::numeric_limits<float>::max();
        float ymax = -std::numeric_limits<float>::max();

        const OrderedCaloHitList *pCaloHitList = &pCluster->GetOrderedCaloHitList();
        for( OrderedCaloHitList::const_iterator caloListIter = pCaloHitList->begin();  caloListIter != pCaloHitList->end(); ++caloListIter ){
            for( CaloHitList::const_iterator hIter = caloListIter->second->begin(); hIter != caloListIter->second->end(); ++hIter ){
                const CaloHit *const pCaloHit = *hIter;

                const CartesianVector &hit(pCaloHit->GetPositionVector());

                if( hit.GetX() < xmin || hit.GetX() > xmax || 
                    hit.GetZ() < zmin || hit.GetZ() > zmax )
                    continue;

                ymin = std::min(hit.GetY(), ymin);
                ymax = std::max(hit.GetY(), ymax);
            }
        }

        float clusterLength = std::sqrt(
            std::pow(xmax - xmin, 2) + 
            std::pow(ymax - ymin, 2) + 
            std::pow(zmax - zmin, 2)
        );

        float min_clusterLength = 3 * LArGeometryHelper::GetWirePitch(this->GetPandora(), m_view); 
        if( clusterLength < min_clusterLength ){
            ClusterList cc;
            cc.push_back(pCluster);

            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&cc, "cluster", RED);
            std::cout << "SKIPPED! Pseudo layer: " << 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            continue;
        }

        /*
	    if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < m_minClusterLayers){
            ClusterList cc;
            cc.push_back(pCluster);

            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&cc, "cluster", RED);
            std::cout << "SKIPPED! Pseudo layer: " << 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            continue;
	    }

        ClusterList cc;
        cc.push_back(pCluster);

        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&cc, "cluster", RED);
        std::cout << "Passed!" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
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

bool ThreeDAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
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

bool ThreeDAssociationAlgorithm::AreClustersAssociated(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
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

    /*
    ClusterList innerClusters, outerClusters;
    innerClusters.push_back(pInnerCluster);
    outerClusters.push_back(pOuterCluster);

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&innerClusters, "InnerClusters", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&outerClusters, "OuterClusters", BLUE);
    */

    if (this->AreClustersAssociated(innerEndCentroid, outerStartCentroid, innerEndFit, outerStartFit)){
        //std::cout <<"CLUSTERS ARE ASSOCIATED!!!" << std::endl;
	//PandoraMonitoringApi::ViewEvent(this->GetPandora());

        return true;
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDAssociationAlgorithm::AreClustersAssociated(const CartesianVector &innerClusterEnd,
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

    const CartesianVector outerStartFit1(outerFit.GetIntercept() +
        outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(outerClusterStart - outerFit.GetIntercept())));
    const CartesianVector outerStartFit2(innerFit.GetIntercept() +
        innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(outerClusterStart - innerFit.GetIntercept())));

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

StatusCode ThreeDAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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
