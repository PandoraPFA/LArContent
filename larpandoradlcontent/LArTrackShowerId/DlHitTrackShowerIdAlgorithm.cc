/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower id algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <chrono>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlHitTrackShowerIdAlgorithm::DlHitTrackShowerIdAlgorithm() :
    m_imageHeight(256), m_imageWidth(256), m_tileSize(128.f), m_visualize(false), m_useTrainingMode(false), m_trainingOutputFile("")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlHitTrackShowerIdAlgorithm::~DlHitTrackShowerIdAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::Run()
{
    if (m_useTrainingMode)
        return this->Train();
    else
        return this->Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::Train()
{
    const int SHOWER{1}, TRACK{2};
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

        const HitType view{pCaloHitList->front()->GetHitType()};

        if (!(view == TPC_VIEW_U || view == TPC_VIEW_V || view == TPC_VIEW_W))
            return STATUS_CODE_NOT_ALLOWED;

        std::string trainingOutputFileName(m_trainingOutputFile);

        if (view == TPC_VIEW_U)
            trainingOutputFileName += "_CaloHitListU.csv";
        else if (view == TPC_VIEW_V)
            trainingOutputFileName += "_CaloHitListV.csv";
        else if (view == TPC_VIEW_W)
            trainingOutputFileName += "_CaloHitListW.csv";

        LArMCParticleHelper::PrimaryParameters parameters;
        // Only care about reconstructability with respect to the current view, so skip good view check
        parameters.m_minHitsForGoodView = 0;
        // Turn off max photo propagation for now, only care about killing off daughters of neutrons
        parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

        LArMvaHelper::MvaFeatureVector featureVector;
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            int tag{TRACK};
            float inputEnergy{0.f};

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                // Throw away non-reconstructable hits
                if (targetMCParticleToHitsMap.find(pMCParticle) == targetMCParticleToHitsMap.end())
                    continue;
                if (LArMCParticleHelper::IsDescendentOf(pMCParticle, 2112))
                    continue;
                inputEnergy = pCaloHit->GetInputEnergy();
                if (inputEnergy < 0.f)
                    continue;

                const int pdg{std::abs(pMCParticle->GetParticleId())};
                if (pdg == 11 || pdg == 22)
                    tag = SHOWER;
                else
                    tag = TRACK;
            }
            catch (const StatusCodeException &)
            {
                continue;
            }

            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
            featureVector.push_back(static_cast<double>(tag));
            featureVector.push_back(static_cast<double>(inputEnergy));
        }
        // Add number of hits to end of vector than rotate (more efficient than direct insert at front)
        featureVector.push_back(static_cast<double>(featureVector.size() / 4));
        std::rotate(featureVector.rbegin(), featureVector.rbegin() + 1, featureVector.rend());

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(trainingOutputFileName, true, featureVector));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::Infer()
{
    const float eps{1.1920929e-7}; // Python float epsilon, used in image padding

    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }

    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        const HitType view{pCaloHitList->front()->GetHitType()};

        if (!(view == TPC_VIEW_U || view == TPC_VIEW_V || view == TPC_VIEW_W))
            return STATUS_CODE_NOT_ALLOWED;

        LArDLHelper::TorchModel &model{view == TPC_VIEW_U ? m_modelU : (view == TPC_VIEW_V ? m_modelV : m_modelW)};

        // Get bounds of hit region
        float xMin{};
        float xMax{};
        float zMin{};
        float zMax{};
        this->GetHitRegion(*pCaloHitList, xMin, xMax, zMin, zMax);
        const float xRange = (xMax + eps) - (xMin - eps);
        int nTilesX = static_cast<int>(std::ceil(xRange / m_tileSize));

        PixelToTileMap sparseMap;
        this->GetSparseTileMap(*pCaloHitList, xMin, zMin, nTilesX, sparseMap);
        const int nTiles = sparseMap.size();

        CaloHitList trackHits, showerHits, otherHits;
        // Process tile
        // ATTN: Be sure to reset all values to zero after each tile has been processed
        float **weights = new float *[m_imageHeight];
        for (int r = 0; r < m_imageHeight; ++r)
            weights[r] = new float[m_imageWidth]();
        for (int i = 0; i < nTiles; ++i)
        {
            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                const float x(pCaloHit->GetPositionVector().GetX());
                const float z(pCaloHit->GetPositionVector().GetZ());
                // Determine which tile the hit will be assigned to
                const int tileX = static_cast<int>(std::floor((x - xMin) / m_tileSize));
                const int tileZ = static_cast<int>(std::floor((z - zMin) / m_tileSize));
                const int tile = sparseMap.at(tileZ * nTilesX + tileX);
                if (tile == i)
                {
                    // Determine hit position within the tile
                    const float localX = std::fmod(x - xMin, m_tileSize);
                    const float localZ = std::fmod(z - zMin, m_tileSize);
                    // Determine hit pixel within the tile
                    const int pixelX = static_cast<int>(std::floor(localX * m_imageWidth / m_tileSize));
                    const int pixelZ = (m_imageHeight - 1) - static_cast<int>(std::floor(localZ * m_imageHeight / m_tileSize));
                    weights[pixelZ][pixelX] += pCaloHit->GetInputEnergy();
                }
            }

            // Find min and max charge to allow normalisation
            float chargeMin{std::numeric_limits<float>::max()}, chargeMax{-std::numeric_limits<float>::max()};
            for (int r = 0; r < m_imageHeight; ++r)
            {
                for (int c = 0; c < m_imageWidth; ++c)
                {
                    if (weights[r][c] > chargeMax)
                        chargeMax = weights[r][c];
                    if (weights[r][c] < chargeMin)
                        chargeMin = weights[r][c];
                }
            }
            float chargeRange{chargeMax - chargeMin};
            if (chargeRange <= 0.f)
                chargeRange = 1.f;

            // Populate accessor based on normalised weights
            CaloHitToPixelMap caloHitToPixelMap;
            LArDLHelper::TorchInput input;
            LArDLHelper::InitialiseInput({1, 1, m_imageHeight, m_imageWidth}, input);
            auto accessor = input.accessor<float, 4>();
            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                const float x(pCaloHit->GetPositionVector().GetX());
                const float z(pCaloHit->GetPositionVector().GetZ());
                // Determine which tile the hit will be assigned to
                const int tileX = static_cast<int>(std::floor((x - xMin) / m_tileSize));
                const int tileZ = static_cast<int>(std::floor((z - zMin) / m_tileSize));
                const int tile = sparseMap.at(tileZ * nTilesX + tileX);
                if (tile == i)
                {
                    // Determine hit position within the tile
                    const float localX = std::fmod(x - xMin, m_tileSize);
                    const float localZ = std::fmod(z - zMin, m_tileSize);
                    // Determine hit pixel within the tile
                    const int pixelX = static_cast<int>(std::floor(localX * m_imageWidth / m_tileSize));
                    const int pixelZ = (m_imageHeight - 1) - static_cast<int>(std::floor(localZ * m_imageHeight / m_tileSize));
                    accessor[0][0][pixelZ][pixelX] = (weights[pixelZ][pixelX] - chargeMin) / chargeRange;
                    caloHitToPixelMap.insert(std::make_pair(pCaloHit, std::make_tuple(tileZ, tileX, pixelZ, pixelX)));
                }
            }
            // Reset weights
            for (int r = 0; r < this->m_imageHeight; ++r)
                for (int c = 0; c < this->m_imageWidth; ++c)
                    weights[r][c] = 0.f;

            // Run the input through the trained model and get the output accessor
            LArDLHelper::TorchInputVector inputs;
            inputs.push_back(input);
            LArDLHelper::TorchOutput output;
            LArDLHelper::Forward(model, inputs, output);
            auto outputAccessor = output.accessor<float, 4>();

            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                auto found{caloHitToPixelMap.find(pCaloHit)};
                if (found == caloHitToPixelMap.end())
                    continue;
                auto pixelMap = found->second;
                const int tileZ(std::get<0>(pixelMap));
                const int tileX(std::get<1>(pixelMap));
                const int tile = sparseMap.at(tileZ * nTilesX + tileX);
                if (tile == i)
                { // Make sure we're looking at a hit in the correct tile
                    const int pixelZ(std::get<2>(pixelMap));
                    const int pixelX(std::get<3>(pixelMap));

                    // Apply softmax to loss to get actual probability
                    float probShower = exp(outputAccessor[0][1][pixelZ][pixelX]);
                    float probTrack = exp(outputAccessor[0][2][pixelZ][pixelX]);
                    float probNull = exp(outputAccessor[0][0][pixelZ][pixelX]);
                    if (probShower > probTrack && probShower > probNull)
                        showerHits.push_back(pCaloHit);
                    else if (probTrack > probShower && probTrack > probNull)
                        trackHits.push_back(pCaloHit);
                    else
                        otherHits.push_back(pCaloHit);
                    float recipSum = 1.f / (probShower + probTrack);
                    // Adjust probabilities to ignore null hits and update LArCaloHit
                    probShower *= recipSum;
                    probTrack *= recipSum;
                    LArCaloHit *pLArCaloHit{const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(pCaloHit))};
                    pLArCaloHit->SetShowerProbability(probShower);
                    pLArCaloHit->SetTrackProbability(probTrack);
                }
            }
        }
        for (int r = 0; r < this->m_imageHeight; ++r)
            delete[] weights[r];
        delete[] weights;

        if (m_visualize)
        {
            const std::string trackListName("TrackHits_" + listName);
            const std::string showerListName("ShowerHits_" + listName);
            const std::string otherListName("OtherHits_" + listName);
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHits, trackListName, BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHits, showerListName, RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &otherHits, otherListName, BLACK));
        }
    }

    if (m_visualize)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax)
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x(pCaloHit->GetPositionVector().GetX());
        const float z(pCaloHit->GetPositionVector().GetZ());
        if (x < xMin)
            xMin = x;
        if (x > xMax)
            xMax = x;
        if (z < zMin)
            zMin = z;
        if (z > zMax)
            zMax = z;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::GetSparseTileMap(
    const CaloHitList &caloHitList, const float xMin, const float zMin, const int nTilesX, PixelToTileMap &sparseMap)
{
    // Identify the tiles that actually contain hits
    std::map<int, bool> tilePopulationMap;
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x(pCaloHit->GetPositionVector().GetX());
        const float z(pCaloHit->GetPositionVector().GetZ());
        // Determine which tile the hit will be assigned to
        const int tileX = static_cast<int>(std::floor((x - xMin) / m_tileSize));
        const int tileZ = static_cast<int>(std::floor((z - zMin) / m_tileSize));
        const int tile = tileZ * nTilesX + tileX;
        tilePopulationMap.insert(std::make_pair(tile, true));
    }

    int nextTile = 0;
    for (auto element : tilePopulationMap)
    {
        if (element.second)
        {
            sparseMap.insert(std::make_pair(element.first, nextTile));
            ++nextTile;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        bool modelLoaded{false};
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameU", m_modelFileNameU));
        if (!m_modelFileNameU.empty())
        {
            m_modelFileNameU = LArFileHelper::FindFileInPath(m_modelFileNameU, "FW_SEARCH_PATH");
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_modelFileNameU, m_modelU));
            modelLoaded = true;
        }
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameV", m_modelFileNameV));
        if (!m_modelFileNameV.empty())
        {
            m_modelFileNameV = LArFileHelper::FindFileInPath(m_modelFileNameV, "FW_SEARCH_PATH");
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_modelFileNameV, m_modelV));
            modelLoaded = true;
        }
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameW", m_modelFileNameW));
        if (!m_modelFileNameW.empty())
        {
            m_modelFileNameW = LArFileHelper::FindFileInPath(m_modelFileNameW, "FW_SEARCH_PATH");
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_modelFileNameW, m_modelW));
            modelLoaded = true;
        }
        if (!modelLoaded)
        {
            std::cout << "Error: Inference requested, but no model files were successfully loaded" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_imageHeight));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_imageWidth));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TileSize", m_tileSize));
    if (m_imageHeight <= 0.f || m_imageWidth <= 0.f || m_tileSize <= 0.f)
    {
        std::cout << "Error: Invalid image size specification" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
