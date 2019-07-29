/**
 *  @file   larpandoracontent/LArDeepLearning/DeepLearningTrackShowerIdAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower id algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArDeepLearning/DeepLearningTrackShowerIdAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DeepLearningTrackShowerIdAlgorithm::DeepLearningTrackShowerIdAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningTrackShowerIdAlgorithm::Run()
{
//    const CaloHitList *pCaloHitList = nullptr;
//    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
//    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
