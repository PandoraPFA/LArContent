/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional hit creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ThreeDHitCreationAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            ParticleFlowObject *pPfo = *iter;

            Cluster *pClusterU(NULL), *pClusterV(NULL), *pClusterW(NULL);
            this->GetClusters(pPfo, pClusterU, pClusterV, pClusterW);

            TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
            LArClusterHelper::LArTwoDSlidingFit(pClusterU, m_slidingFitWindow, slidingFitResultU);
            LArClusterHelper::LArTwoDSlidingFit(pClusterV, m_slidingFitWindow, slidingFitResultV);
            LArClusterHelper::LArTwoDSlidingFit(pClusterW, m_slidingFitWindow, slidingFitResultW);

            CaloHitList caloHitListU, caloHitListV, caloHitListW;
            pClusterU->GetOrderedCaloHitList().GetCaloHitList(caloHitListU);
            pClusterV->GetOrderedCaloHitList().GetCaloHitList(caloHitListV);
            pClusterW->GetOrderedCaloHitList().GetCaloHitList(caloHitListW);

            CaloHitList newThreeDHits, ommittedTwoDHits;
            this->CreateThreeDHits(caloHitListU, slidingFitResultV, slidingFitResultW, newThreeDHits, ommittedTwoDHits);
            this->CreateThreeDHits(caloHitListV, slidingFitResultU, slidingFitResultW, newThreeDHits, ommittedTwoDHits);
            this->CreateThreeDHits(caloHitListW, slidingFitResultU, slidingFitResultV, newThreeDHits, ommittedTwoDHits);

            Cluster *pCluster3D(NULL);
            this->CreateThreeDCluster(newThreeDHits, pCluster3D);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster3D));

            CaloHitList extrapolatedHits;
            this->CreateExtrapolatedHits(ommittedTwoDHits, pCluster3D, extrapolatedHits);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster3D, &extrapolatedHits));
        }
        catch (StatusCodeException &)
        {
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::GetClusters(const ParticleFlowObject *const pPfo, Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW) const
{
    pClusterU = NULL; pClusterV = NULL; pClusterW = NULL;
    const ClusterList &clusterList(pPfo->GetClusterList());

    if (3 != clusterList.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);
        const HitType hitType(LArThreeDHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V == hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        Cluster *&pTargetCluster((TPC_VIEW_U == hitType) ? pClusterU : (TPC_VIEW_V == hitType) ? pClusterV : pClusterW);
        pTargetCluster = pCluster;
    }

    if ((NULL == pClusterU) || (NULL == pClusterV) || (NULL == pClusterW))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateThreeDHits(const CaloHitList &inputTwoDHits, const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    CaloHitList &newThreeDHits, CaloHitList &ommittedTwoDHits) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateThreeDCluster(const CaloHitList &caloHitList, Cluster *&pCluster) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateExtrapolatedHits(const CaloHitList &ommittedHits, Cluster *pCluster, CaloHitList &extrapolatedHits) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));

    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
