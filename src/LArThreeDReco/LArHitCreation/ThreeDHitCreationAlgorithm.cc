/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional hit creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

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
//            Pfo *pPfo = *iter;
//
//            Cluster *pClusterU(NULL), *pClusterV(NULL), *pClusterW(NULL);
//            this->GetClusters(pPfo, pClusterU, pClusterV, pClusterW);
//
//            TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
//            LArClusterHelper::LArTwoDSlidingFit(pClusterU, m_slidingFitWindow, slidingFitResultU);
//            LArClusterHelper::LArTwoDSlidingFit(pClusterV, m_slidingFitWindow, slidingFitResultV);
//            LArClusterHelper::LArTwoDSlidingFit(pClusterW, m_slidingFitWindow, slidingFitResultW);
//
//            CaloHitList caloHitListU, caloHitListV, caloHitListW;
//            pClusterU->GetOrderedCaloHitList().GetCaloHitList(caloHitListU);
//            pClusterV->GetOrderedCaloHitList().GetCaloHitList(caloHitListV);
//            pClusterW->GetOrderedCaloHitList().GetCaloHitList(caloHitListW);
//
//            // TODO, need a way of obtaining address of newly created hits
//            CaloHitList newThreeDHits, ommittedTwoDHits;
//            this->CreateThreeDHits(caloHitListU, slidingFitResultV, slidingFitResultW, newThreeDHits, ommittedTwoDHits);
//            this->CreateThreeDHits(caloHitListU, slidingFitResultV, slidingFitResultW, newThreeDHits, ommittedTwoDHits);
//            this->CreateThreeDHits(caloHitListU, slidingFitResultV, slidingFitResultW, newThreeDHits, ommittedTwoDHits);
//
//            Cluster *pCluster3D(NULL);
//            this->CreateCluster(pPfo, newThreeDHits, pCluster3D);
//
//            CaloHitList extrapolatedHits;
//            this->GetExtrapolatedHits(pCluster3D, extrapolatedHits);
//
//            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster3D, &extrapolatedHits));
//
//            // Sort out lists
        }
        catch (StatusCodeException &)
        {
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
