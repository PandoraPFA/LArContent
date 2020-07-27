/**
 *  @file   larpandoracontent/LArDeepLearning/DeepLearningTrackShowerIdAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower id algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArDeepLearning/DeepLearningTrackShowerIdAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <chrono>

using namespace pandora;

namespace lar_content
{

DeepLearningTrackShowerIdAlgorithm::DeepLearningTrackShowerIdAlgorithm() :
    m_xMin(-420),
    m_zMinU(-350),
    m_zMinV(0),
    m_zMinW(-25),
    m_nBins(512),
    m_visualize(false),
    m_useTrainingMode(false),
    m_profile(false),
    m_trainingOutputFile("")
{
    const float span(980);
    m_xMax = m_xMin + span;
    m_zMaxU = m_zMinU + span;
    m_zMaxV = m_zMinV + span;
    m_zMaxW = m_zMinW + span;
}

DeepLearningTrackShowerIdAlgorithm::~DeepLearningTrackShowerIdAlgorithm()
{
    if (m_profile)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "cpu_tree", "cpu.root", "UPDATE"));
        }
        catch(const StatusCodeException&)
        {
            std::cout << "DeepLearningTrackShowerIdAlgorithm: Unable to write tree to file" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool isNeutronDaughter(const MCParticle *const mcParticle)
{
    const MCParticle* particle = mcParticle;
    while(!particle->GetParentList().empty())
    {   
        if(particle->GetParentList().size() > 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const MCParticle* pParent = *(particle->GetParentList().begin());
        const int pdg = std::abs(pParent->GetParticleId());
        if (pdg == 2112)
        {
            return true;
        }   
        particle = pParent;
    }

    return false;
}


StatusCode DeepLearningTrackShowerIdAlgorithm::Run()
{
    if (m_useTrainingMode)
        return Train();
    else
        return Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "TrainingOutputFileName", m_trainingOutputFile));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", m_modelFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Visualize", m_visualize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Profile", m_profile));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageXMin", m_xMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageXMax", m_xMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMinU", m_zMinU));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMaxU", m_zMaxU));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMinV", m_zMinV));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMaxV", m_zMaxV));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMinW", m_zMinW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMaxW", m_zMaxW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumberOfBins", m_nBins));

    return STATUS_CODE_SUCCESS;
}

StatusCode DeepLearningTrackShowerIdAlgorithm::Train()
{
    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

        if (!isU && !isV && !isW) return STATUS_CODE_NOT_ALLOWED;

        std::string trainingOutputFileName(m_trainingOutputFile);
        LArMvaHelper::MvaFeatureVector featureVector;

        if (isU) trainingOutputFileName += "_CaloHitListU.txt";
        else if (isV) trainingOutputFileName += "_CaloHitListV.txt";
        else if (isW) trainingOutputFileName += "_CaloHitListW.txt";

        featureVector.push_back(static_cast<double>(pCaloHitList->size()));

        int discardedHits = 0;
        LArMCParticleHelper::PrimaryParameters parameters;
        // Only care about reconstructability with respect to the current view, so skip good view check
        parameters.m_minHitsForGoodView = 0;
        // Turn off max photo propagation for now, only care about killing off daughters of neutrons
        parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList,
                parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            int isReconstructable = 1;
            int pdg(-1);

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                // Throw away non-reconstructable hits
                if (targetMCParticleToHitsMap.find(pMCParticle) == targetMCParticleToHitsMap.end())
                {
                    ++discardedHits;
                    isReconstructable = 0;
                }
                pdg = pMCParticle->GetParticleId();
                if (isReconstructable)
                {
                    if(isNeutronDaughter(pMCParticle))
                    {
                        isReconstructable = 0;
                    }
                }
            }
            catch (...)
            {
                continue;
            }

            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetY()));
            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
            featureVector.push_back(static_cast<double>(pdg));
            featureVector.push_back(static_cast<double>(isReconstructable));
        }

        LArMvaHelper::ProduceTrainingExample(trainingOutputFileName, true, featureVector);
    }
    return STATUS_CODE_SUCCESS;
}

void GetHitRegion(const CaloHitList& caloHitList, float& xMin, float& xMax, float& zMin, float& zMax)
{
    xMin = std::numeric_limits<float>::max(); xMax = 0.f;
    zMin = std::numeric_limits<float>::max(); zMax = 0.f;
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x(pCaloHit->GetPositionVector().GetX());
        const float z(pCaloHit->GetPositionVector().GetZ());
        if (x < xMin) xMin = x;
        if (x > xMax) xMax = x;
        if (z < zMin) zMin = z;
        if (z > zMax) zMax = z;
    }
}

void GetSparseTileMap(const CaloHitList& caloHitList, const float xMin,
        const float zMin, const float tileSize, const int nTilesX, std::map<int, int>& sparseMap)
{
    // Identify the tiles that actually contain hits
    std::map<int, bool> tilePopulationMap;
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x(pCaloHit->GetPositionVector().GetX());
        const float z(pCaloHit->GetPositionVector().GetZ());
        // Determine which tile the hit will be assigned to
        const int tileX = static_cast<int>(std::floor((x - xMin) / tileSize));
        const int tileZ = static_cast<int>(std::floor((z - zMin) / tileSize));
        const int tile = tileZ * nTilesX + tileX;
        tilePopulationMap.insert(std::make_pair(tile, true));
    }

    int nextTile = 0;
    for(auto element : tilePopulationMap)
    {
        if (element.second)
        {
            sparseMap.insert(std::make_pair(element.first, nextTile));
            ++nextTile;
        }
    }
}

StatusCode DeepLearningTrackShowerIdAlgorithm::Infer()
{
    auto start = std::chrono::steady_clock::now();
    torch::jit::script::Module pModule;
    try
    {
        pModule = torch::jit::load(m_modelFileName);
    }
    catch (...)
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

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

        const float tileSize = 128.f;
        const int imageWidth = 256;
        const int imageHeight = 256;
        // Get bounds of hit region
        float xMin{}; float xMax{}; float zMin{}; float zMax{};
        GetHitRegion(*pCaloHitList, xMin, xMax, zMin, zMax);
        const float xRange = xMax - xMin;
        const float zRange = zMax - zMin;
        int nTilesX = static_cast<int>(std::ceil(xRange / tileSize));
        int nTilesZ = static_cast<int>(std::ceil(zRange / tileSize));
        // Need to add 1 to number of tiles for ranges exactly matching tile size and for zero ranges
        if (std::fmod(xRange, tileSize) == 0.f) ++nTilesX;
        if (std::fmod(zRange, tileSize) == 0.f) ++nTilesZ;

        std::map<int, int> sparseMap;
        GetSparseTileMap(*pCaloHitList, xMin, zMin, tileSize, nTilesX, sparseMap);
        const int nTiles = sparseMap.size();

        // Normalisation value from training set - should store and load this from file
        const float PIXEL_VALUE = (255.0 - 0.558497428894043) / 11.920382499694824;
        CaloHitList trackHits, showerHits, otherHits;
        // Populate tensor - may need to split the batch into separate single image batches
        //torch::Tensor input = torch::zeros({nTiles, 1, imageHeight, imageWidth});
        for (int i = 0; i < nTiles; ++i)
        {
            torch::Tensor input = torch::zeros({1, 1, imageHeight, imageWidth});
            typedef std::map<const CaloHit*, std::tuple<int, int, int, int>> CaloHitToPixelMap;
            CaloHitToPixelMap caloHitToPixelMap;
            auto accessor = input.accessor<float, 4>();
            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                const float x(pCaloHit->GetPositionVector().GetX());
                const float z(pCaloHit->GetPositionVector().GetZ());
                // Determine which tile the hit will be assigned to
                const int tileX = static_cast<int>(std::floor((x - xMin) / tileSize));
                const int tileZ = static_cast<int>(std::floor((z - zMin) / tileSize));
                const int tile = sparseMap.at(tileZ * nTilesX + tileX);
                if (tile == i)
                {
                    // Determine hit position within the tile
                    const float localX = std::fmod(x - xMin, tileSize);
                    const float localZ = std::fmod(z - zMin, tileSize);
                    // Determine hit pixel within the tile
                    const int pixelX = static_cast<int>(std::floor(localX * imageHeight / tileSize));
                    const int pixelZ = static_cast<int>(std::floor(localZ * imageWidth / tileSize));
                    //accessor[tile][0][pixelX][pixelZ] = PIXEL_VALUE;
                    accessor[0][0][pixelX][pixelZ] = PIXEL_VALUE;
                    caloHitToPixelMap.insert(std::make_pair(pCaloHit, std::make_tuple(tileX, tileZ, pixelX, pixelZ)));
                }
            }

            // Pass as input the input Tensor containing the calo hit picture
            std::vector<torch::jit::IValue> inputs;
            inputs.push_back(input);

            // Run the input through the trained model and get the output accessor
            at::Tensor output = pModule.forward(inputs).toTensor();
            auto outputAccessor = output.accessor<float, 4>();

            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                auto found = caloHitToPixelMap.find(pCaloHit);
                if (found == caloHitToPixelMap.end()) continue;
                auto pixelMap = found->second;
                //auto pixelMap = caloHitToPixelMap.at(pCaloHit);
                const int tileX(std::get<0>(pixelMap));
                const int tileZ(std::get<1>(pixelMap));
                const int tile = sparseMap.at(tileZ * nTilesX + tileX);
                if (tile == i)
                {
                    const int pixelX(std::get<2>(pixelMap));
                    const int pixelZ(std::get<3>(pixelMap));

                    object_creation::CaloHit::Metadata metadata;
                    // Apply softmax to loss to get actual probability
                    //float probShower = exp(outputAccessor[tile][1][pixelX][pixelZ]);
                    //float probTrack = exp(outputAccessor[tile][2][pixelX][pixelZ]);
                    //float probNull = exp(outputAccessor[tile][0][pixelX][pixelZ]);
                    float probShower = exp(outputAccessor[0][1][pixelX][pixelZ]);
                    float probTrack = exp(outputAccessor[0][2][pixelX][pixelZ]);
                    float probNull = exp(outputAccessor[0][0][pixelX][pixelZ]);
                    float recipSum = 1.f / (probShower + probTrack + probNull);
                    probShower *= recipSum;
                    probTrack *= recipSum;
                    probNull *= recipSum;
                    metadata.m_propertiesToAdd["Pshower"] = probShower;
                    metadata.m_propertiesToAdd["Ptrack"] = probTrack;
                    if (probShower > probTrack && probShower > probNull)
                        showerHits.push_back(pCaloHit);
                    else if (probTrack > probShower && probTrack > probNull)
                        trackHits.push_back(pCaloHit);
                    else
                        otherHits.push_back(pCaloHit);
                    const StatusCode &statusCode(PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
                    if (statusCode != STATUS_CODE_SUCCESS)
                        std::cout << "Cannot set calo hit meta data" << std::endl;
                }
            }
        }

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

    auto end = std::chrono::steady_clock::now();
    if (m_profile)
    {
        auto diff = end - start;
        float cpu(std::chrono::duration<double, std::milli>(diff).count());
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "cpu_tree", "cpu_time", cpu));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "cpu_tree"));
    }

    if (m_visualize)
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
