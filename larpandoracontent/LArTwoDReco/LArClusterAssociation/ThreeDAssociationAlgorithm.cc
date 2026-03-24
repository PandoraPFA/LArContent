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

using namespace std::literals;
using namespace pandora;

namespace lar_content
{

ThreeDAssociationAlgorithm::ThreeDAssociationAlgorithm() :
    m_minClusterLayers(4),
    m_maxGapLayers(7),
    m_fitLayers(30),
    m_maxGapDistanceSquared(2.f),
    m_minCosRelativeAngle(0.94f),
    m_minClusterLength(1.f),
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

        float clusterLength = LArClusterHelper::GetLength(pCluster);

        if( clusterLength < m_minClusterLength )
            continue;

        this->PCAFit(pCluster);

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

                /*
            if( !(iterJ - clusterVector.begin() == 27) || (iterI - clusterVector.begin() == 27) ){
                continue;
                this->VisualizeClusters(pInnerCluster, pOuterCluster, "Consider deeezz"sv);
            } else {
                continue;
            }
                */

            if (!this->AreClustersAssociated(pInnerCluster, pOuterCluster))
                continue;

            //this->VisualizeClusters(pInnerCluster, pOuterCluster, "Accepted"sv );

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

    float closestDistance = LArClusterHelper::GetClosestDistance(pInnerCluster, pOuterCluster);

    if( closestDistance > 4 * m_minClusterLength )
        return false;

    const CartesianVector innerDirection = clusterToPCAMap.at(pInnerCluster).direction;
    const CartesianVector outerDirection = clusterToPCAMap.at(pOuterCluster).direction;
    const CartesianVector innerCentroid = clusterToPCAMap.at(pInnerCluster).centroid;
    const CartesianVector outerCentroid = clusterToPCAMap.at(pOuterCluster).centroid;

    const float openingAngle{ innerDirection.GetCosOpeningAngle(outerDirection) };

    if( openingAngle < m_minCosRelativeAngle )
        return false;

    const auto innerCaloHitList = clusterToHitListMap.at( pInnerCluster );

    constexpr float float_limit = std::numeric_limits<float>::max();

    float innerMinAvgDir{ +float_limit };
    float innerMaxAvgDir{ -float_limit };
    float innerMinDir{ +float_limit };
    float innerMaxDir{ -float_limit };

    // Project the start and end of the first and second clusters on the average direction
    // and their respective principal axis
    CartesianVector averageDirection{ ((innerDirection + outerDirection) * 0.5).GetUnitVector() };

    for( const CaloHit *const pCaloHit : innerCaloHitList ){ 
        innerMinAvgDir = std::min(innerMinAvgDir, averageDirection.GetDotProduct(pCaloHit->GetPositionVector()));
        innerMaxAvgDir = std::max(innerMaxAvgDir, averageDirection.GetDotProduct(pCaloHit->GetPositionVector()));

        innerMinDir = std::min(
            innerMinDir, 
            innerDirection.GetDotProduct(pCaloHit->GetPositionVector() - innerCentroid)
        );
        innerMaxDir = std::max(
            innerMaxDir, 
            innerDirection.GetDotProduct(pCaloHit->GetPositionVector() - innerCentroid)
        );
    }

    const auto outerCaloHitList = clusterToHitListMap.at( pOuterCluster );

    float outerMinAvgDir{ +float_limit };
    float outerMaxAvgDir{ -float_limit };
    float outerMinDir{ +float_limit };
    float outerMaxDir{ -float_limit };

    for( const CaloHit *const pCaloHit : outerCaloHitList ){
        outerMinAvgDir = std::min(outerMinAvgDir, averageDirection.GetDotProduct(pCaloHit->GetPositionVector()));
        outerMaxAvgDir = std::max(outerMaxAvgDir, averageDirection.GetDotProduct(pCaloHit->GetPositionVector()));

        outerMinDir = std::min(
            outerMinDir, 
            outerDirection.GetDotProduct(pCaloHit->GetPositionVector() - outerCentroid)
        );
        outerMaxDir = std::max(
            outerMaxDir, 
            outerDirection.GetDotProduct(pCaloHit->GetPositionVector() - outerCentroid)
        );
    }

    if( innerMinAvgDir > innerMaxAvgDir || outerMinAvgDir > outerMaxAvgDir || 
        innerMinDir > innerMaxDir || outerMinDir > outerMaxDir )
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // Overlapping clusters

    /* 
    * C1: *----------*
    * C2:       *--------* 
    */
    if( (outerMinAvgDir + 0.4f > innerMinAvgDir && outerMinAvgDir + 0.4f < innerMaxAvgDir) || 
        (innerMinAvgDir + 0.4f > outerMinAvgDir && innerMinAvgDir + 0.4f < outerMinAvgDir) )
        return false;

    /*
    * C1:       *--------*
    * C2:  *--------*
    */
    if( (outerMaxAvgDir - 0.4f > innerMinAvgDir && outerMaxAvgDir - 0.4f < innerMaxAvgDir) || 
        (innerMaxAvgDir - 0.4f > outerMinAvgDir && innerMaxAvgDir - 0.4f < outerMaxAvgDir) )
        return false;

    CartesianVector outerClusterStart(1.0f, 1.0f, 1.0f);
    CartesianVector innerClusterEnd(1.0f, 1.0f, 1.0f);

    if( innerMaxAvgDir > outerMaxAvgDir ){
        // inner -> outer; outer -> inner
        outerClusterStart = innerCentroid + innerDirection * innerMinDir;
        innerClusterEnd = outerCentroid + outerDirection * outerMaxDir;
    } else {
        outerClusterStart = outerCentroid + outerDirection * outerMinDir;
        innerClusterEnd = innerCentroid + innerDirection * innerMaxDir;
    }

    if( (outerClusterStart - innerClusterEnd).GetMagnitudeSquared() > m_maxGapDistanceSquared ){
        //this->VisualizeClusters(pInnerCluster, pOuterCluster, "Rejected"sv );
        
        return false;
    }
    
    //this->VisualizeClusters(pInnerCluster, pOuterCluster, "Accepted"sv );

    /*
    float endCluster1{ -std::numeric_limits<float>::max() };
    float endCluster2{ +std::numeric_limits<float>::max() };

    // Compute intercept and endpoints
    for( const CaloHit *const pCaloHit : innerCaloHitList ){
        endCluster1 = 
            std::max(endCluster1, innerDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - innerCentroid));
        endCluster2 = 
            std::min(endCluster2, outerDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - outerCentroid));
    }

    const CartesianVector innerClusterEnd1 = innerCentroid + innerDirection*endCluster1;
    const CartesianVector innerClusterEnd2 = outerCentroid + outerDirection*endCluster2;

    float startCluster1{ +std::numeric_limits<float>::max() };
    float startCluster2{ -std::numeric_limits<float>::max() };

    // Compute intercept and endpoints
    for( const CaloHit *const pCaloHit : outerCaloHitList ){
        startCluster1 = 
            std::min(startCluster1, outerDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - outerCentroid));
        startCluster2 = 
            std::max(startCluster2, innerDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - innerCentroid));
    }

    const CartesianVector outerClusterStart1 = outerCentroid + outerDirection*startCluster1;
    const CartesianVector outerClusterStart2 = innerCentroid + innerDirection*startCluster2;
    */

    return true;
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

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDAssociationAlgorithm::PCAFit( const pandora::Cluster *const pCluster ) const
{
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);   

    CartesianVector centroid (1.0f, 1.0f, 1.0f);
    CartesianVector direction(1.0f, 1.0f, 1.0f);

    this->PCAFit(clusterHits, centroid, direction);

    clusterToHitListMap.insert(
        std::make_pair(pCluster, clusterHits)
    );

    clusterToPCAMap.insert(
        std::make_pair(pCluster, PCAAttr{centroid, direction})
    );

}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDAssociationAlgorithm::PCAFit( const CaloHitList& mergedClusterCaloHit,  CartesianVector &centroid, CartesianVector &direction ) const
{
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;

    LArPcaHelper::RunPca(mergedClusterCaloHit, centroid, eigenValues, eigenVecs);

    // Compute direction: Eigen vector that has a positive z-component -> By convention
    direction = eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f;
}


void ThreeDAssociationAlgorithm::VisualizeClusters(const pandora::Cluster *const pInner, const pandora::Cluster *const pOuter, 
    const std::string_view str) const
{
    ClusterList innerClusters, outerClusters;
    innerClusters.push_back(pInner);
    outerClusters.push_back(pOuter);

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&innerClusters, "InnerClusters", AUTOITER);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&outerClusters, "OuterClusters", AUTOITER);

    std::cout << str << "\n";

	PandoraMonitoringApi::ViewEvent(this->GetPandora());
}


void ThreeDAssociationAlgorithm::PrintCluster( const pandora::Cluster *const pCluster ) const{
    std::cout << "  Nhits: " << pCluster->GetNCaloHits() << "\n";
    std::cout << "  Length: " << LArClusterHelper::GetLength(pCluster) << "\n";
}

} // namespace lar_content
