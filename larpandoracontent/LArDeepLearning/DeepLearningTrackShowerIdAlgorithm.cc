/**
 *  @file   larpandoracontent/LArDeepLearning/DeepLearningTrackShowerIdAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower id algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include <torch/script.h>

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
    std::shared_ptr<torch::jit::script::Module> module = torch::jit::load(m_modelFileName);

    assert(module != nullptr);
    std::cout << "ok\n";

//    const CaloHitList *pCaloHitList = nullptr;
//    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
//    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", m_modelFileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
