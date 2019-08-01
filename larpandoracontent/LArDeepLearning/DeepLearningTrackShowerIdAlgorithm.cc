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

DeepLearningTrackShowerIdAlgorithm::DeepLearningTrackShowerIdAlgorithm() :
    m_xMin(-700),
    m_xMax(700),
    m_zMin(-400),
    m_zMax(1000),
    m_nBins(1024),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningTrackShowerIdAlgorithm::Run()
{
    // Load the model.pt file.
    std::shared_ptr<torch::jit::script::Module> pModule(nullptr);

    try
    {
        pModule = torch::jit::load(m_modelFileName);
    }
    catch (const c10::Error &e)
    {
        std::cout << "Error loading the PyTorch module" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    if (m_visualize)
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        const float xSpan(m_xMax - m_xMin), zSpan(m_zMax - m_zMin);

        typedef std::map<const CaloHit*, std::pair<int, int>> CaloHitToBinMap;
        CaloHitToBinMap caloHitToBinMap;

        // Start with RGB picture of black pixels.  Four indices: first default size 1, second index is RGB indices, third is xBin,
        // fourth is zBin
        torch::Tensor input = torch::zeros({1, 3, m_nBins, m_nBins});
        auto accessor = input.accessor<float, 4>();

        // Create a map of calo hits to x/z bin values.  Set the output track shower id of the pixel using the RGB values at the pixel
        // containing the calo hit
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());

            const int xBin(std::floor((x-m_xMin)*m_nBins/xSpan));
            const int zBin(std::floor((z-m_zMin)*m_nBins/zSpan));

            // ATTN: Set pixels containing a calo hit to white
            if (xBin >= 0 && xBin <= m_nBins && zBin >= 0 && zBin <= m_nBins)
            {
                caloHitToBinMap.insert(std::make_pair(pCaloHit, std::make_pair(xBin, zBin)));
                accessor[0][0][xBin][zBin] = 1;
                accessor[0][1][xBin][zBin] = 1;
                accessor[0][2][xBin][zBin] = 1;
            }
        }

        // Pass as input the input Tensor containing the calo hit picture
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(input);

        // Run the input through the trained model and get the output accessor
        at::Tensor output = pModule->forward(inputs).toTensor();
        auto outputAccessor = output.accessor<float, 4>();

        // Colour in the shower and track bits (and other) in a visual display for first performance inspection
        CaloHitList showerHits, trackHits, other;

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            if (caloHitToBinMap.find(pCaloHit) == caloHitToBinMap.end())
            {
                other.push_back(pCaloHit);
                continue;
            }

            const int xBin(caloHitToBinMap.at(pCaloHit).first);
            const int zBin(caloHitToBinMap.at(pCaloHit).second);

            // Is the R value bigger than the B value.  In training the target picture was coloured such that showers were red and tracks blue
            const bool isTrack(outputAccessor[0][0][xBin][zBin] > outputAccessor[0][2][xBin][zBin] ? false : true);
            object_creation::CaloHit::Metadata metadata;

            if (isTrack)
            {
                trackHits.push_back(pCaloHit);
            }
            else
            {
                showerHits.push_back(pCaloHit);
            }
        }

        if (m_visualize)
        {
            const std::string trackListName("TrackHits_" + listName);
            const std::string showerListName("ShowerHits_" + listName);
            const std::string otherListName("OtherHits_" + listName);
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHits, trackListName, BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHits, showerListName, RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &other, otherListName, BLACK));
        }
    }

    if (m_visualize)
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", m_modelFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Visualize", m_visualize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageXMin", m_xMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageXMax", m_xMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMin", m_zMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMax", m_zMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumberOfBins", m_nBins));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
