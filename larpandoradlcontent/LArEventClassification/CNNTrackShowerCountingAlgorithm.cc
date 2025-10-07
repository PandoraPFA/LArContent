/**
 *  @file   larpandoradlcontent/LArEventClassification/CNNTrackShowerCountingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track/shower counting algorithm
 *
 *  $Log: $
 */

#include <chrono>
#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"

#include "larpandoradlcontent/LArEventClassification/CNNTrackShowerCountingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

CNNTrackShowerCountingAlgorithm::CNNTrackShowerCountingAlgorithm() :
    m_trainingMode{false},
    m_trainingTreeName{""},
    m_trainingFileName{""},
    m_height{512},
    m_width{512},
    m_wiresPerPixel{1},
    m_driftStep{0.5f},
    m_goodMCTrackHits{5},
    m_goodMCShowerHits{5},
    m_mcHitWeightThreshold{0.5f},
    m_secondaryDistanceThreshold{2.5f},
    m_minHits{10},
    m_maxChargeThreshold{10.f},
    m_outputPfoListName{""}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

CNNTrackShowerCountingAlgorithm::~CNNTrackShowerCountingAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName, m_trainingFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "CNNTrackShowerCountingAlgorithm: Failed to save ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CNNTrackShowerCountingAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CNNTrackShowerCountingAlgorithm::PrepareTrainingSample()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const MCParticle *pMCNeutrino{nullptr};
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(pMCParticle))
        {
            pMCNeutrino = pMCParticle;
            break;
        }
    }
    if (pMCNeutrino == nullptr)
        return STATUS_CODE_FAILURE;

    const int nuPDG{pMCNeutrino->GetParticleId()};
    const float nuEnergy{pMCNeutrino->GetEnergy()};
    const CartesianVector nuTrueVertex{pMCNeutrino->GetVertex()};

    MCParticleList mcPrimaries;
    this->GetMCPrimaries(pMCParticleList, mcPrimaries);
    const InteractionDescriptor &intType{LArInteractionTypeHelper::GetInteractionDescriptor(mcPrimaries)};
    const bool isCC{intType.IsCC()};

    std::unordered_map<HitType, unsigned int> nTracksPerView;
    std::unordered_map<HitType, unsigned int> nShowersPerView;
    bool emptyView{false};
    for (const std::string &name : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, name, pCaloHitList));
        if (pCaloHitList->empty())
        {
            emptyView = true;
            break;
        }

        LArMCParticleHelper::MCContributionMap mcToHitsMap;
        this->GetVisibleParticles(pMCParticleList, pMCNeutrino, pCaloHitList, mcToHitsMap);

        if (mcToHitsMap.empty())
        {
            emptyView = true;
            break;
        }

        const HitType hitType{pCaloHitList->front()->GetHitType()};
        unsigned int nTracks{0}, nShowers{0};
        this->CountMCPrimaries(mcToHitsMap, nTracks, nShowers);
        nTracksPerView[hitType] = nTracks;
        nShowersPerView[hitType] = nShowers;
    }

    if (emptyView)
    {
        std::cout << "Skipping event as it does not have enough hits or associated primary particles to make a training sample" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuPDG", nuPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "isCC", static_cast<int>(isCC)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuEnergy", nuEnergy));

    const float nuTrueVertexX{nuTrueVertex.GetX()};
    const float nuTrueVertexY{nuTrueVertex.GetY()};
    const float nuTrueVertexZ{nuTrueVertex.GetZ()};
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuTrueVertexX", nuTrueVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuTrueVertexY", nuTrueVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuTrueVertexZ", nuTrueVertexZ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nTracksU", static_cast<int>(nTracksPerView[TPC_VIEW_U])));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nShowersU", static_cast<int>(nShowersPerView[TPC_VIEW_U])));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nTracksV", static_cast<int>(nTracksPerView[TPC_VIEW_V])));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nShowersV", static_cast<int>(nShowersPerView[TPC_VIEW_V])));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nTracksW", static_cast<int>(nTracksPerView[TPC_VIEW_W])));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nShowersW", static_cast<int>(nShowersPerView[TPC_VIEW_W])));

    ViewToVertexPositionMap nuRecoVertices;
    this->GetVertexPositions(nuRecoVertices);
    std::vector<int> pixel_view;
    std::vector<int> pixel_row;
    std::vector<int> pixel_column;
    std::vector<float> pixel_charge;

    // Create the pixel map and write out the features
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        VertexPosition &vertex(nuRecoVertices.at(pCaloHitList->front()->GetHitType()));
        PixelMap pixelMap;
        this->CreatePixelMap(pCaloHitList, pixelMap, vertex);
        if (pixelMap.size() < m_minHits)
            continue;

        const int view{static_cast<int>(pCaloHitList->front()->GetHitType())};
        for (auto const &pixel : pixelMap)
        {
            pixel_view.emplace_back(view);
            pixel_row.emplace_back(pixel.first.first);
            pixel_column.emplace_back(pixel.first.second);
            pixel_charge.emplace_back(pixel.second <= m_maxChargeThreshold ? pixel.second : m_maxChargeThreshold);
        }
    }

    const int nuRecoVertexDriftBinU(nuRecoVertices.at(TPC_VIEW_U).GetDriftBin());
    const int nuRecoVertexDriftBinV(nuRecoVertices.at(TPC_VIEW_V).GetDriftBin());
    const int nuRecoVertexDriftBinW(nuRecoVertices.at(TPC_VIEW_W).GetDriftBin());
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuRecoVertexDriftBinU", nuRecoVertexDriftBinU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuRecoVertexDriftBinV", nuRecoVertexDriftBinV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuRecoVertexDriftBinW", nuRecoVertexDriftBinW));

    const int nuRecoVertexWireBinU(nuRecoVertices.at(TPC_VIEW_U).GetWireBin());
    const int nuRecoVertexWireBinV(nuRecoVertices.at(TPC_VIEW_V).GetWireBin());
    const int nuRecoVertexWireBinW(nuRecoVertices.at(TPC_VIEW_W).GetWireBin());
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuRecoVertexWireBinU", nuRecoVertexWireBinU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuRecoVertexWireBinV", nuRecoVertexWireBinV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuRecoVertexWireBinW", nuRecoVertexWireBinW));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "pixelView", &pixel_view));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "pixelRow", &pixel_row));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "pixelColumn", &pixel_column));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "pixelCharge", &pixel_charge));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CNNTrackShowerCountingAlgorithm::Infer()
{
    ViewToVertexPositionMap vertices;
    this->GetVertexPositions(vertices);

    TrackShowerCountingResults result;
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        const HitType hitType(pCaloHitList->front()->GetHitType());

        VertexPosition &vertex(vertices.at(hitType));
        PixelMap pixelMap;
        this->CreatePixelMap(pCaloHitList, pixelMap, vertex);
        if (pixelMap.size() < m_minHits)
            continue;

        // Image first
        LArDLHelper::TorchInput imageInput;
        this->MakeNetworkInputFromPixelMap(pixelMap, imageInput);

        // Now for the reconstructed vertex
        LArDLHelper::TorchInput vtxInput;
        LArDLHelper::InitialiseInput({1, 2}, vtxInput);
        auto accessor = vtxInput.accessor<float, 2>();
        accessor[0][0] = vertices[hitType].GetDriftBin() / static_cast<float>(m_height);
        accessor[0][1] = vertices[hitType].GetWireBin() / static_cast<float>(m_width);

        LArDLHelper::TorchInputVector inputs;
        inputs.emplace_back(imageInput);
        inputs.emplace_back(vtxInput);
        LArDLHelper::TorchMultiOutput output;
        LArDLHelper::Forward(m_model, inputs, output);

        // Network has two outputs that are returned in Python as a list. Here that will be a tuple containing
        // two tensors, one for each of the outputs.
        const auto outputList = output.toList();

        // Due to the way PyTorch includes the softmax in the loss function we need to apply it here if we
        // want out scores to sum to one
        const auto nuOutput = torch::softmax(outputList.get(0).toTensor(), 1);
        const auto trkOutput = torch::softmax(outputList.get(1).toTensor(), 1);
        const auto shwOutput = torch::softmax(outputList.get(2).toTensor(), 1);

        result.AddScoresFromView(hitType, nuOutput, trkOutput, shwOutput);
    }
    return this->StorePredictions(result);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CNNTrackShowerCountingAlgorithm::StorePredictions(const TrackShowerCountingResults &result)
{
    // For now we will create a fake PFParticle and store results as metadata.
    // TODO: update this after changes to the SDK will allow us to store event-level information
    const PfoList *pDummyEventPfoList{nullptr};
    std::string dummyName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pDummyEventPfoList, dummyName));

    PandoraContentApi::ParticleFlowObject::Parameters dummyEventPfoParameters;
    dummyEventPfoParameters.m_particleId = NU_MU;
    dummyEventPfoParameters.m_charge = PdgTable::GetParticleCharge(dummyEventPfoParameters.m_particleId.Get());
    dummyEventPfoParameters.m_mass = PdgTable::GetParticleMass(dummyEventPfoParameters.m_particleId.Get());
    dummyEventPfoParameters.m_energy = 0.f;
    dummyEventPfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    dummyEventPfoParameters.m_propertiesToAdd["IsNeutrino"] = 1.f;
    const ParticleFlowObject *pDummyEventPfo{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, dummyEventPfoParameters, pDummyEventPfo));

    // Now we add some metadata to the Pfo
    object_creation::ParticleFlowObject::Metadata dummyEventPfoMetadata;

    const std::vector<std::string> nuClassNames{"predNuScoreNC", "predNuScoreCCNumu", "predNuScoreCCNue"};
    const std::vector<HitType> views{TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

    for (const HitType &view : views)
    {
        const std::string viewString{view == TPC_VIEW_U ? "U" : (view == TPC_VIEW_V ? "V" : "W")};

        const std::vector<float> nuScores{result.GetNuScoresFromView(view)};
        for (unsigned int n = 0; n < nuClassNames.size(); ++n)
        {
            const std::string &nuClass{nuClassNames.at(n)};
            dummyEventPfoMetadata.m_propertiesToAdd[nuClass + viewString] = nuScores.at(n);
        }

        const std::vector<float> trackScores{result.GetTrackScoresFromView(view)};
        for (unsigned int t = 0; t < trackScores.size(); ++t)
        {
            const std::string &className{"predTrackScore" + std::to_string(t)};
            dummyEventPfoMetadata.m_propertiesToAdd[className + viewString] = trackScores.at(t);
        }

        const std::vector<float> showerScores{result.GetShowerScoresFromView(view)};
        for (unsigned int s = 0; s < showerScores.size(); ++s)
        {
            const std::string &className{"predShowerScore" + std::to_string(s)};
            dummyEventPfoMetadata.m_propertiesToAdd[className + viewString] = showerScores.at(s);
        }
    }

    // Save the best class values too
    dummyEventPfoMetadata.m_propertiesToAdd["predNuClassU"] = result.GetNuClassPredictionFromView(TPC_VIEW_U);
    dummyEventPfoMetadata.m_propertiesToAdd["predNuClassV"] = result.GetNuClassPredictionFromView(TPC_VIEW_V);
    dummyEventPfoMetadata.m_propertiesToAdd["predNuClassW"] = result.GetNuClassPredictionFromView(TPC_VIEW_W);

    dummyEventPfoMetadata.m_propertiesToAdd["nPredTracksU"] = result.GetTrackClassPredictionFromView(TPC_VIEW_U);
    dummyEventPfoMetadata.m_propertiesToAdd["nPredTracksV"] = result.GetTrackClassPredictionFromView(TPC_VIEW_V);
    dummyEventPfoMetadata.m_propertiesToAdd["nPredTracksW"] = result.GetTrackClassPredictionFromView(TPC_VIEW_W);

    dummyEventPfoMetadata.m_propertiesToAdd["nPredShowersU"] = result.GetShowerClassPredictionFromView(TPC_VIEW_U);
    dummyEventPfoMetadata.m_propertiesToAdd["nPredShowersV"] = result.GetShowerClassPredictionFromView(TPC_VIEW_V);
    dummyEventPfoMetadata.m_propertiesToAdd["nPredShowersW"] = result.GetShowerClassPredictionFromView(TPC_VIEW_W);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pDummyEventPfo, dummyEventPfoMetadata));

    if (!pDummyEventPfoList->empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, dummyName, m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetVisibleParticles(const MCParticleList *const pMCParticleList, const MCParticle *const pMCNeutrino,
    const CaloHitList *const pCaloHitList, LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    // Find all MCParticles contributing > m_mcHitWeightThreshold to at least one CaloHit
    for (auto const *pCaloHit : *pCaloHitList)
    {
        for (auto const &mcHitWeightPair : pCaloHit->GetMCParticleWeightMap())
        {
            if (mcHitWeightPair.second > m_mcHitWeightThreshold)
                mcToHitsMap[mcHitWeightPair.first].emplace_back(pCaloHit);
        }
    }

    // We need to add any MC particles without hits into the map
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (!mcToHitsMap.count(pMCParticle))
            mcToHitsMap[pMCParticle] = CaloHitList();
    }

    // Apply the number of hits cut to primaries and set failed primaries to be removed
    std::vector<const MCParticle *> keysToRemove;
    for (auto const &mcToHitsPair : mcToHitsMap)
    {
        const MCParticle *pMCParticle{mcToHitsPair.first};
        const unsigned int nHitsCut{this->IsTracklike(pMCParticle) ? m_goodMCTrackHits : m_goodMCShowerHits};
        //const MCParticle *pParentParticle{*(pMCParticle->GetParentList().begin())};

        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle) && mcToHitsPair.second.size() < nHitsCut)
            keysToRemove.emplace_back(pMCParticle);
    }

    // Look for secondaries produced close to the vertex when their parent is not visible
    const HitType hitType{pCaloHitList->front()->GetHitType()};
    for (auto const &mcToHitsPair : mcToHitsMap)
    {
        const MCParticle *pMCParticle{mcToHitsPair.first};

        // Ignore primary particles this time
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))
            continue;

        // Ignore the incoming neutrino
        if (pMCParticle->GetParentList().empty())
        {
            keysToRemove.emplace_back(pMCParticle);
            continue;
        }

        // Look for secondaries produced by primaries that are considered not visible and see if they start close to the neutrino vertex position
        const MCParticle *pParentParticle{*(pMCParticle->GetParentList().begin())};
        const unsigned int nHitsCut{this->IsTracklike(pMCParticle) ? m_goodMCTrackHits : m_goodMCShowerHits};
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pParentParticle) &&
            std::find(keysToRemove.begin(), keysToRemove.end(), pParentParticle) != keysToRemove.end())
        {
            if (mcToHitsPair.second.size() < nHitsCut)
            {
                keysToRemove.emplace_back(pMCParticle);
                continue;
            }

            const CartesianVector thisVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCParticle->GetVertex(), hitType));
            const CartesianVector nuVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCNeutrino->GetVertex(), hitType));
            const float distance{std::sqrt(nuVertex.GetDistanceSquared(thisVertex))};

            if (distance > m_secondaryDistanceThreshold)
                keysToRemove.emplace_back(pMCParticle);
        }
        else
            keysToRemove.emplace_back(pMCParticle);
    }

    // We will be left with visible primaries and selected secondaries of invisible primaries
    for (auto const &key : keysToRemove)
        mcToHitsMap.erase(key);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CNNTrackShowerCountingAlgorithm::MakeNetworkInputFromPixelMap(const PixelMap &pixelMap, LArDLHelper::TorchInput &networkInput) const
{
    LArDLHelper::InitialiseInput({1, 1, m_height, m_width}, networkInput);
    auto accessor = networkInput.accessor<float, 4>();

    for (auto const &pixel : pixelMap)
    {
        const unsigned int row{pixel.first.first};
        const unsigned int col{pixel.first.second};
        const float q{pixel.second <= m_maxChargeThreshold ? pixel.second : m_maxChargeThreshold};
        // Rescale to be between 0 and 1.
        accessor[0][0][row][col] = q / m_maxChargeThreshold;
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::CreatePixelMap(const CaloHitList *const pCaloHitList, PixelMap &pixelMap, VertexPosition &vertex) const
{
    float xMin{std::numeric_limits<float>::max()}, xMax{-1.f * std::numeric_limits<float>::max()};
    float zMin{std::numeric_limits<float>::max()}, zMax{-1.f * std::numeric_limits<float>::max()};
    this->GetHitRegion(pCaloHitList, xMin, xMax, zMin, zMax);

    // Get the wire pitch for this view
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view{pCaloHitList->front()->GetHitType()};
    const float pitch{view == TPC_VIEW_U ? pTPC->GetWirePitchU() : (view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW())};
    const float wireStep{pitch * m_wiresPerPixel};

    const float zSpan{pitch * m_wiresPerPixel * m_width};
    const float xSpan{m_driftStep * m_height};

    // If the min / max values exceed the above spans then we need to crop
    const bool cropZ{(zMax - zMin) > zSpan ? true : false};
    const bool cropX{(xMax - xMin) > xSpan ? true : false};
    const float vtxPosX{vertex.GetPosition().GetX()};
    const float vtxPosZ{vertex.GetPosition().GetZ()};

    // Adjust the min and max depending on reconstructed vertex position for events where
    // the reconstructed vertex isn't on a CaloHit, e.g. NC pizero events
    xMin = std::min(xMin, vtxPosX);
    xMax = std::max(xMax, vtxPosX);
    zMin = std::min(zMin, vtxPosZ);
    zMax = std::max(zMax, vtxPosZ);

    constexpr float eps{0.0001f};
    if (cropX)
        this->GetCrop1D(pCaloHitList, vertex, xMin, xMax, xSpan, true);
    else
        vertex.AddVertexBin(true, static_cast<unsigned int>(std::floor(eps + (vtxPosX - xMin) / m_driftStep)));

    if (cropZ)
        this->GetCrop1D(pCaloHitList, vertex, zMin, zMax, zSpan, false);
    else
        vertex.AddVertexBin(false, static_cast<unsigned int>(std::floor(eps + (vtxPosZ - zMin) / wireStep)));

    // Now we have the cropped region we can make the pixel map
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const float driftPos{pCaloHit->GetPositionVector().GetX()};
        if (driftPos < xMin || driftPos > xMax)
            continue;
        const float wirePos{pCaloHit->GetPositionVector().GetZ()};
        if (wirePos < zMin || wirePos > zMax)
            continue;

        const unsigned int driftBin{static_cast<unsigned int>(std::floor(eps + (driftPos - xMin) / m_driftStep))};
        const unsigned int wireBin{static_cast<unsigned int>(std::floor(eps + (wirePos - zMin) / wireStep))};
        const Pixel pixel{driftBin, wireBin};
        if (!pixelMap.count(pixel))
            pixelMap[pixel] = pCaloHit->GetMipEquivalentEnergy();
        else
            pixelMap[pixel] += pCaloHit->GetMipEquivalentEnergy();
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetVertexPositions(ViewToVertexPositionMap &vertices) const
{
    const VertexList *pVertexList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    if (pVertexList->empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    if (pVertexList->front()->GetVertexType() != VERTEX_3D)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    vertices.emplace(
        TPC_VIEW_U, VertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_U)));
    vertices.emplace(
        TPC_VIEW_V, VertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_V)));
    vertices.emplace(
        TPC_VIEW_W, VertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_W)));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetCrop1D(const CaloHitList *const pCaloHitList, VertexPosition &vertexPosition, float &min,
    float &max, const float &span, const bool &isDrift) const
{
    const unsigned int imageDimension{isDrift ? m_height : m_width};
    const float step{span / imageDimension};
    // Make a vector of bins from min to max using the required pitch
    const unsigned int nBins{static_cast<unsigned int>(std::ceil((max - min) / step))};
    if (nBins <= imageDimension)
        return;
    std::vector<float> dimensionBins(nBins, 0.f);
    constexpr float eps{0.0001f};
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const CartesianVector pos{pCaloHit->GetPositionVector()};
        const float dim{(isDrift ? pos.GetX() : pos.GetZ())};
        unsigned int bin{static_cast<unsigned int>(std::floor(eps + (dim - min) / step))};
        // Protection for final bin upper edge
        if (std::fabs(dim - max) < std::numeric_limits<float>::epsilon())
            bin = nBins - 1;
        if (bin >= nBins)
        {
            std::cout << " - Error! Bin " << bin << " outside range 0 to " << nBins - 1 << std::endl;
            std::cout << "  - Underlying values: " << dim - max << ", " << dim << ", " << min << ", " << max << ", " << step << std::endl;
            continue;
        }
        dimensionBins[bin] += pCaloHit->GetMipEquivalentEnergy();
    }

    const float vtxPos{isDrift ? vertexPosition.GetPosition().GetX() : vertexPosition.GetPosition().GetZ()};
    const unsigned int vtxBin{static_cast<unsigned int>(std::floor(eps + (vtxPos - min) / step))};

    // Look for the required number of consecutive bins with the highest integrated charge
    std::unordered_map<unsigned int, float> startBinToCharge;
    for (unsigned int i = 0; i < (nBins - imageDimension); ++i)
    {
        const bool vtxInRange{(vtxBin >= i) && (vtxBin < (i + imageDimension))};
        startBinToCharge[i] = (vtxInRange ? std::accumulate(dimensionBins.begin() + i, dimensionBins.begin() + i + imageDimension, 0.f) : 0.f);
    }

    // Get the maximum element from the map
    const std::unordered_map<unsigned int, float>::iterator maxElementIter{
        std::max_element(startBinToCharge.begin(), startBinToCharge.end(), [](const auto &a, const auto &b) { return a.second < b.second; })};
    const unsigned int maxBin{maxElementIter->first};
    const float localMin{min + maxBin * step};
    const float localMax{localMin + imageDimension * step};

    vertexPosition.AddVertexBin(isDrift, vtxBin - maxBin);
    min = localMin;
    max = localMax;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetMCPrimaries(const MCParticleList *pMCParticleList, MCParticleList &mcPrimaries) const
{
    for (const MCParticle *const particle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(particle))
            mcPrimaries.push_back(particle);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool CNNTrackShowerCountingAlgorithm::IsTracklike(const MCParticle *const pMCParticle) const
{
    return !(this->IsShowerlike(pMCParticle));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool CNNTrackShowerCountingAlgorithm::IsShowerlike(const MCParticle *const pMCParticle) const
{
    const int pdg{pMCParticle->GetParticleId()};
    return (std::abs(pdg) != 11 && pdg != 22) ? false : true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::CountMCPrimaries(
    const LArMCParticleHelper::MCContributionMap &mcToHitsMap, unsigned int &nTracks, unsigned int &nShowers) const
{
    nTracks = 0;
    nShowers = 0;
    for (const auto &[pMCParticle, _] : mcToHitsMap)
    {
        if (this->IsShowerlike(pMCParticle))
            ++nShowers;
        else
            ++nTracks;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetHitRegion(const CaloHitList *const pCaloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = std::numeric_limits<float>::lowest();
    zMin = std::numeric_limits<float>::max();
    zMax = std::numeric_limits<float>::lowest();

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        xMin = std::min(x, xMin);
        xMax = std::max(x, xMax);
        zMin = std::min(z, zMin);
        zMax = std::max(z, zMax);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CNNTrackShowerCountingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_height));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_width));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WiresPerPixel", m_wiresPerPixel));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DriftStep", m_driftStep));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingTreeName", m_trainingTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
    }
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_model);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GoodMCTrackHits", m_goodMCTrackHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GoodMCShowerHits", m_goodMCShowerHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCHitWeightThreshold", m_mcHitWeightThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SecondaryDistanceThreshold", m_secondaryDistanceThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHits", m_minHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxChargeThreshold", m_maxChargeThreshold));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));
    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

CNNTrackShowerCountingAlgorithm::VertexPosition::VertexPosition() :
    m_position{0.f, 0.f, 0.f},
    m_driftBin{0},
    m_wireBin{0}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

CNNTrackShowerCountingAlgorithm::VertexPosition::VertexPosition(const CartesianVector &pos) :
    m_position{pos},
    m_driftBin{0},
    m_wireBin{0}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::VertexPosition::AddVertexBin(const bool isDrift, const unsigned int bin)
{
    if (isDrift)
        m_driftBin = bin;
    else
        m_wireBin = bin;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

CartesianVector CNNTrackShowerCountingAlgorithm::VertexPosition::GetPosition() const
{
    return m_position;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CNNTrackShowerCountingAlgorithm::VertexPosition::GetDriftBin() const
{
    return m_driftBin;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CNNTrackShowerCountingAlgorithm::VertexPosition::GetWireBin() const
{
    return m_wireBin;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::TensorToVector(const torch::Tensor &scores) const
{
    // These are technically 2D as they have a batch dimension
    const auto scoreAccessor = scores.accessor<float, 2>();
    std::vector<float> scoreVector;
    for (unsigned int element = 0; element < scores.size(1); ++element)
        scoreVector.emplace_back(scoreAccessor[0][element]);

    return scoreVector;
}

void CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::AddScoresFromView(
    const HitType &view, const torch::Tensor &nuScores, const torch::Tensor &trackScores, const torch::Tensor &showerScores)
{
    if (m_nuScores.count(view) || m_trackScores.count(view) || m_showerScores.count(view))
    {
        std::cerr << "TrackShowerCountingResults::AddScoresFromView: scores for view " << view << " already exist. Doing nothing." << std::endl;
        return;
    }
    m_nuScores[view] = this->TensorToVector(nuScores);
    m_trackScores[view] = this->TensorToVector(trackScores);
    m_showerScores[view] = this->TensorToVector(showerScores);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::GetNuScoresFromView(const pandora::HitType &view) const
{
    return m_nuScores.at(view);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::GetNuClassPredictionFromView(const pandora::HitType &view) const
{
    const std::vector<float> &nuScores{m_nuScores.at(view)};
    return std::distance(nuScores.begin(), std::max_element(nuScores.begin(), nuScores.end()));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::GetTrackScoresFromView(const pandora::HitType &view) const
{
    return m_trackScores.at(view);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::GetTrackClassPredictionFromView(const pandora::HitType &view) const
{
    const std::vector<float> &trkScores{m_trackScores.at(view)};
    return std::distance(trkScores.begin(), std::max_element(trkScores.begin(), trkScores.end()));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::GetShowerScoresFromView(const pandora::HitType &view) const
{
    return m_showerScores.at(view);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

unsigned int CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::GetShowerClassPredictionFromView(const pandora::HitType &view) const
{
    const std::vector<float> &shwScores{m_showerScores.at(view)};
    return std::distance(shwScores.begin(), std::max_element(shwScores.begin(), shwScores.end()));
}

} // namespace lar_dl_content
