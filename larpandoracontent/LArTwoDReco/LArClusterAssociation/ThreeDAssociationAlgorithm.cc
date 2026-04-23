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
    m_minClusterLength(1.1f),
    m_minCosRelativeAngle(0.94f),
    m_maxGapDistanceSquared(10.f),
    m_view(TPC_3D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    if (!pClusterList->empty())
        m_view = LArClusterHelper::GetClusterHitType(pClusterList->front());

    clusterAttrMap.clear();

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;
        float clusterLength = LArClusterHelper::GetLength(pCluster);

        if( clusterLength < m_minClusterLength )
            continue;

        this->PCAFit(pCluster);

        clusterVector.push_back(pCluster);

        CartesianVector innerCoordinate(0.f, 0.f, 0.f);
        CartesianVector outerCoordinate(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);

        clusterAttrMap[pCluster].length = clusterLength;
        clusterAttrMap[pCluster].innerCoordinate = innerCoordinate;
        clusterAttrMap[pCluster].outerCoordinate = outerCoordinate;
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
    CartesianVector currentClusterInnerCoordinate{ clusterAttrMap.at(pCurrentCluster).innerCoordinate };
    CartesianVector currentClusterOuterCoordinate{ clusterAttrMap.at(pCurrentCluster).outerCoordinate };

    CartesianVector testClusterInnerCoordinate{ clusterAttrMap.at(pTestCluster).innerCoordinate };
    CartesianVector testClusterOuterCoordinate{ clusterAttrMap.at(pTestCluster).outerCoordinate };

    const CartesianVector currentPosition{ isForward ? currentClusterOuterCoordinate : currentClusterInnerCoordinate };
    const CartesianVector testPosition{ isForward ? testClusterOuterCoordinate : testClusterInnerCoordinate };

    if( isForward && ((LArClusterHelper::SortCoordinatesByPosition(currentPosition, testPosition) || 
                      ((currentPosition == testPosition) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster)))))
        return true;

    if( !isForward && ((LArClusterHelper::SortCoordinatesByPosition(testPosition, currentPosition) || 
                      ((currentPosition == testPosition) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster)))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDAssociationAlgorithm::AreClustersAssociated(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
{

    const float closestDistance = LArClusterHelper::GetClosestDistance(pInnerCluster, pOuterCluster);

    if( closestDistance > m_maxGapDistanceSquared )
        return false;

    const CartesianVector innerClusterInnerCoordinate{ clusterAttrMap.at(pInnerCluster).innerCoordinate };
    const CartesianVector innerClusterOuterCoordinate{ clusterAttrMap.at(pInnerCluster).outerCoordinate };
    const CartesianVector outerClusterInnerCoordinate{ clusterAttrMap.at(pOuterCluster).innerCoordinate };
    const CartesianVector outerClusterOuterCoordinate{ clusterAttrMap.at(pOuterCluster).outerCoordinate };

    if( !LArClusterHelper::SortCoordinatesByPosition(innerClusterInnerCoordinate, outerClusterInnerCoordinate) )
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    const LArPcaHelper::EigenVectors innerVectors{ clusterAttrMap.at(pInnerCluster).eigenVectors };
    const LArPcaHelper::EigenVectors outerVectors{ clusterAttrMap.at(pOuterCluster).eigenVectors };

    const CartesianVector innerPrimaryAxis  { innerVectors.at(0).GetZ() > 0.f ? innerVectors.at(0) : innerVectors.at(0) * -1.f };
    const CartesianVector innerSecondaryAxis{ innerVectors.at(1).GetZ() > 0.f ? innerVectors.at(1) : innerVectors.at(1) * -1.f };
    const CartesianVector innerTertiaryAxis { innerVectors.at(2).GetZ() > 0.f ? innerVectors.at(2) : innerVectors.at(2) * -1.f };

    const CartesianVector outerPrimaryAxis  { outerVectors.at(0).GetZ() > 0.f ? outerVectors.at(0) : outerVectors.at(0) * -1.f };
    const CartesianVector outerSecondaryAxis{ outerVectors.at(1).GetZ() > 0.f ? outerVectors.at(1) : outerVectors.at(1) * -1.f };
    const CartesianVector outerTertiaryAxis { outerVectors.at(2).GetZ() > 0.f ? outerVectors.at(2) : outerVectors.at(2) * -1.f };

    const CartesianVector innerCentroid{ clusterAttrMap.at(pInnerCluster).centroid };
    const CartesianVector outerCentroid{ clusterAttrMap.at(pOuterCluster).centroid };

    const float openingAngle{ innerPrimaryAxis.GetCosOpeningAngle(outerPrimaryAxis) };

    constexpr float shortClusterLength = 2.50f;

    // Bypass opening angle cut if a short cluster is the extension of a longer one
    if( (closestDistance > 0.5f) ||
        ((clusterAttrMap.at(pInnerCluster).length <= shortClusterLength ) == (clusterAttrMap.at(pOuterCluster).length <= shortClusterLength))){

        if( openingAngle < m_minCosRelativeAngle )
            return false;
    }

    bool isInnerLongest{ LArClusterHelper::SortByNHits(pInnerCluster, pOuterCluster) && 
        (clusterAttrMap.at(pInnerCluster).length > clusterAttrMap.at(pOuterCluster).length) };
    
    const CartesianVector centroid{ isInnerLongest ? innerCentroid : outerCentroid };

    const CartesianVector primaryAxis  { isInnerLongest ? innerPrimaryAxis   : outerPrimaryAxis   };
    const CartesianVector secondaryAxis{ isInnerLongest ? innerSecondaryAxis : outerSecondaryAxis };
    const CartesianVector tertiaryAxis { isInnerLongest ? innerTertiaryAxis  : outerTertiaryAxis  };

    const CartesianVector innerEndPrimaryAxis  { centroid + primaryAxis * primaryAxis.GetDotProduct(innerClusterOuterCoordinate - centroid) };
    const CartesianVector outerStartPrimaryAxis{ centroid + primaryAxis * primaryAxis.GetDotProduct(outerClusterInnerCoordinate - centroid) };

    const CartesianVector innerEndSecondaryAxis  { centroid + secondaryAxis * secondaryAxis.GetDotProduct(innerClusterOuterCoordinate - centroid) };
    const CartesianVector outerStartSecondaryAxis{ centroid + secondaryAxis * secondaryAxis.GetDotProduct(outerClusterInnerCoordinate - centroid) };
    
    const CartesianVector innerEndTertiaryAxis  { centroid + tertiaryAxis * tertiaryAxis.GetDotProduct(innerClusterOuterCoordinate - centroid) };
    const CartesianVector outerStartTertiaryAxis{ centroid + tertiaryAxis * tertiaryAxis.GetDotProduct(outerClusterInnerCoordinate - centroid) };

    if( (innerEndTertiaryAxis - outerStartTertiaryAxis).GetMagnitudeSquared() > shortClusterLength )
        return false;

    if( (innerEndSecondaryAxis - outerStartSecondaryAxis).GetMagnitudeSquared() > shortClusterLength )
        return false;

    if( (innerEndPrimaryAxis - outerStartPrimaryAxis).GetMagnitudeSquared() < m_maxGapDistanceSquared )
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MaxGapDistanceSquared", m_maxGapDistanceSquared));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDAssociationAlgorithm::PCAFit( const pandora::Cluster *const pCluster ) const
{
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);   

    CartesianVector centroid (0.f, 0.f, 0.f);
    CartesianVector direction(0.f, 0.f, 0.f);

    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;

    LArPcaHelper::RunPca(clusterHits, centroid, eigenValues, eigenVectors);

    clusterAttrMap[pCluster].centroid = centroid;
    clusterAttrMap[pCluster].eigenValues = eigenValues;
    clusterAttrMap[pCluster].eigenVectors = eigenVectors;
}

} // namespace lar_content
