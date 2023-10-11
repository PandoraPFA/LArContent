/**
 *  @file   larpandoracontent/LArCheating/CheatingLongitudinalAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the longitudinal association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArCheating/CheatingLongitudinalAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingLongitudinalAssociationAlgorithm::CheatingLongitudinalAssociationAlgorithm() :
    m_minClusterLayers(4),
    m_maxGapLayers(7),
    m_fitLayers(30),
    m_maxGapDistanceSquared(10.f),
    m_minCosRelativeAngle(0.985f),
    m_maxTransverseDisplacement(2.f),
    m_maxLongitudinalDisplacement(2.f),
    m_hitSizeZ(0.3f),
    m_hitSizeX(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingLongitudinalAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
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

void CheatingLongitudinalAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
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

//-----------------------------------------------------------------------------------------------------------------------------------------
StatusCode CheatingLongitudinalAssociationAlgorithm::GetPrincipalAxisVector() const
{ 
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS,!=,PandoraContentApi::GetCurrentList(*this,pClusterList));
    const std::string m_treename{"mc"};
   
    for(const Cluster *const pCluster : *pClusterList)
    {
        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);   
 
        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecs;
        LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        LArPcaHelper::RunPca(clusterCaloHitList, centroid, eigenValues, eigenVecs);

        if(!eigenVecs.empty())
        {
            const float pcaX{eigenVecs.front().GetX()};
            const float pcaY{eigenVecs.front().GetY()};
            const float pcaZ{eigenVecs.front().GetZ()};	
 
            for(const CartesianVector &axis : eigenVecs)	
       	    {
                std::cout << axis  << std::endl;
            }
     
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "PcaX", pcaX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "PcaY", pcaY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "PcaZ", pcaZ));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    return STATUS_CODE_SUCCESS;    
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingLongitudinalAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
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

bool CheatingLongitudinalAssociationAlgorithm::AreClustersAssociated(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
{
    const MCParticle *const pMCOuterParticle(MCParticleHelper::GetMainMCParticle(pOuterCluster));
    const MCParticle *const pMCInnerParticle(MCParticleHelper::GetMainMCParticle(pInnerCluster)); 
    std::cout << "Cluster, Inner: ";
    std::cout << pMCInnerParticle->GetParticleId() << std::endl;
    std::cout << "Cluster, Outer: " << pMCOuterParticle->GetParticleId() << std::endl;
     
    if (pMCOuterParticle == pMCInnerParticle)
        return true;

    return false;
}
//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingLongitudinalAssociationAlgorithm::AreClustersAssociated(const CartesianVector &innerClusterEnd,
    const CartesianVector &outerClusterStart, const ClusterFitResult &innerFit, const ClusterFitResult &outerFit) const
{
    if (!innerFit.IsFitSuccessful() || !outerFit.IsFitSuccessful())
        return false;

    if (innerFit.GetDirection().GetCosOpeningAngle(outerFit.GetDirection()) < m_minCosRelativeAngle)
        return false;

    const CartesianVector innerEndFit1(
        innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(innerClusterEnd - innerFit.GetIntercept())));
    const CartesianVector innerEndFit2(
        outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(innerClusterEnd - outerFit.GetIntercept())));

    const CartesianVector outerStartFit1(
        outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(outerClusterStart - outerFit.GetIntercept())));
    const CartesianVector outerStartFit2(
        innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(outerClusterStart - innerFit.GetIntercept())));

    const CartesianVector clusterSeparation(outerClusterStart - innerClusterEnd);

    if ((std::fabs(clusterSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) &&
        (std::fabs(clusterSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedSeparation(outerStartFit1 - innerEndFit1);

    if ((std::fabs(fittedSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) &&
        (std::fabs(fittedSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedInnerSeparation(innerEndFit2 - innerEndFit1);

    if ((std::fabs(fittedInnerSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) &&
        (std::fabs(fittedInnerSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedOuterSeparation(outerStartFit2 - outerStartFit1);

    if ((std::fabs(fittedOuterSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) &&
        (std::fabs(fittedOuterSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingLongitudinalAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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
