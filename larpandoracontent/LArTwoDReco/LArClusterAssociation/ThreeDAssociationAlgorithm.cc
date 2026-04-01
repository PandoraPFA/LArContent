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
    m_maxGapDistanceSquared(5.f),
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

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerPosition);
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

    const float closestDistance = LArClusterHelper::GetClosestDistance(pInnerCluster, pOuterCluster);

    if( closestDistance > 5 * m_minClusterLength )
        return false;

    // TODO: Need a better way to hand this as computing extremal coordinates will be expensive
    CartesianVector innerClusterInnerCoordinate(0.f, 0.f, 0.f);
    CartesianVector innerClusterOuterCoordinate(0.f, 0.f, 0.f);
    CartesianVector outerClusterInnerCoordinate(0.f, 0.f, 0.f);
    CartesianVector outerClusterOuterCoordinate(0.f, 0.f, 0.f);

    LArClusterHelper::GetExtremalCoordinates(pInnerCluster, innerClusterInnerCoordinate, innerClusterOuterCoordinate);
    LArClusterHelper::GetExtremalCoordinates(pOuterCluster, outerClusterInnerCoordinate, outerClusterOuterCoordinate);

    if( !LArClusterHelper::SortCoordinatesByPosition(innerClusterInnerCoordinate, outerClusterInnerCoordinate) )
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    const CartesianVector innerDirection = clusterToPCAMap.at(pInnerCluster).direction;
    const CartesianVector outerDirection = clusterToPCAMap.at(pOuterCluster).direction;
    const CartesianVector innerCentroid  = clusterToPCAMap.at(pInnerCluster).centroid;
    const CartesianVector outerCentroid  = clusterToPCAMap.at(pOuterCluster).centroid;

    const float openingAngle{ innerDirection.GetCosOpeningAngle(outerDirection) };

    // Bypass opening angle cut if a short cluster is the extension of a longer one
    if( (closestDistance > m_minClusterLength) ||
        ((LArClusterHelper::GetLength(pInnerCluster) < 5*m_minClusterLength) ==
         (LArClusterHelper::GetLength(pOuterCluster) < 5*m_minClusterLength) )){

        if( openingAngle < m_minCosRelativeAngle )
            return false;
    }

    const CartesianVector innerEndFit1{ innerCentroid +  innerDirection * innerDirection.GetDotProduct(innerClusterOuterCoordinate - innerCentroid) }; 
    const CartesianVector innerEndFit2{ outerCentroid +  outerDirection * outerDirection.GetDotProduct(innerClusterOuterCoordinate - outerCentroid) }; 

    const CartesianVector outerStartFit1{ outerCentroid +  outerDirection * outerDirection.GetDotProduct(outerClusterInnerCoordinate - outerCentroid) }; 
    const CartesianVector outerEndFit1{ outerCentroid +  outerDirection * outerDirection.GetDotProduct(outerClusterOuterCoordinate - outerCentroid) }; 
    const CartesianVector outerStartFit2{ innerCentroid +  innerDirection * innerDirection.GetDotProduct(outerClusterInnerCoordinate - innerCentroid) }; 
    const CartesianVector outerEndFit2{ innerCentroid +  innerDirection * innerDirection.GetDotProduct(outerClusterOuterCoordinate - innerCentroid) }; 

    // Check for overlapping clusters
    if( ((innerEndFit1 - innerCentroid).GetMagnitudeSquared() > (outerStartFit2 - innerCentroid).GetMagnitudeSquared() + 2.56f) ){
        std::cout << (innerEndFit1 - innerCentroid).GetMagnitudeSquared() << "/n";
        std::cout << (outerStartFit2 - innerCentroid).GetMagnitudeSquared() << "/n";

        //this->VisualizeClusters(pInnerCluster, pOuterCluster, "Rejected"sv );

        return false;
    }


        /*
    if( (innerEndFit1 - innerCentroid).GetMagnitudeSquared() > (outerStartFit2 - innerCentroid).GetMagnitudeSquared() + 1.f ){
        //(innerEndFit2 - outerCentroid).GetMagnitudeSquared() < (outerStartFit1 - outerCentroid).GetMagnitudeSquared() + 0.16f ){
        return false;
    }
    */

    //const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), m_view)};

    if( (innerClusterOuterCoordinate - outerClusterInnerCoordinate).GetMagnitudeSquared() < m_maxGapDistanceSquared )
        return true;

    if( (innerEndFit1 - outerStartFit2).GetMagnitudeSquared() < m_maxGapDistanceSquared )
        return true;

    if( (innerEndFit2 - outerStartFit1).GetMagnitudeSquared() < m_maxGapDistanceSquared )
        return true;

    return false;

    /*
    const auto innerCaloHitList = clusterToHitListMap.at( pInnerCluster );

    constexpr float float_limit = std::numeric_limits<float>::max();

    float innerMinAvgDir{ +float_limit };
    float innerMaxAvgDir{ -float_limi

    float innerMinDir{ +float_limit };
    float innerMaxDir{ -float_limit };

    // Project the start and end of the first and second clusters on the average direction
    // akjasdfasdfnd their respective principal axis
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

        if( (outerMinAvgDir + 0.4f > innerMinAvgDir && outerMinAvgDir + 0.4f < innerMaxAvgDir) || 
        (innerMinAvgDir + 0.4f > outerMinAvgDir && innerMinAvgDir + 0.4f < outerMinAvgDir) )
        return false;

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

    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), m_view)};

    if( (outerClusterStart - innerClusterEnd).GetMagnitudeSquared() > ratio * ratio *  m_maxGapDistanceSquared ){
        //this->VisualizeClusters(pInnerCluster, pOuterCluster, "Rejected"sv );
        
        return false;
    }
    
    return true;
    */
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

    clusterToHitListMap.insert(std::make_pair(pCluster, clusterHits));

    clusterToPCAMap.insert(std::make_pair(pCluster, PCAAttr{centroid, direction}));
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

} // namespace lar_content
