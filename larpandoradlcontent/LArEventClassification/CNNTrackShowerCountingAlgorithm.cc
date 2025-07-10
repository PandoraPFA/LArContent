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
//    m_useVertexForCrops{true},
    m_useSimpleTruthLabels{false},
    m_goodMCTrackHits{5},
    m_goodMCShowerHits{5},
    m_mcHitWeightThreshold{0.5f},
    m_secondaryDistanceThreshold{2.5f},
    m_minHits{10},
    m_outputPfoListName{""}
{}

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
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const MCParticle *pMCNeutrino{nullptr};
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(pMCParticle))
            pMCNeutrino = pMCParticle;
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

    std::map<HitType, unsigned int> nTracksPerView;
    std::map<HitType, unsigned int> nShowersPerView;
    bool emptyView{false};
    for (const std::string &name : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, name, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

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
        std::cout << " - Got " << nTracks << " tracks and " << nShowers << " showers in view " << static_cast<int>(hitType) << std::endl;
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

    if (nTracksPerView.count(TPC_VIEW_U))
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nTracksU", static_cast<int>(nTracksPerView[TPC_VIEW_U])));
    if (nShowersPerView.count(TPC_VIEW_U))
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nShowersU", static_cast<int>(nShowersPerView[TPC_VIEW_U])));
    if (nTracksPerView.count(TPC_VIEW_V))
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nTracksV", static_cast<int>(nTracksPerView[TPC_VIEW_V])));
    if (nShowersPerView.count(TPC_VIEW_V))
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nShowersV", static_cast<int>(nShowersPerView[TPC_VIEW_V])));
    if (nTracksPerView.count(TPC_VIEW_W))
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nTracksW", static_cast<int>(nTracksPerView[TPC_VIEW_W])));
    if (nShowersPerView.count(TPC_VIEW_W))
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

//        std::cout << "Creating pixel map for hit collection " << listname << " with nHits = " << pCaloHitList->size() << std::endl;
        VertexPosition &vertex(nuRecoVertices.at(pCaloHitList->front()->GetHitType()));
        PixelMap pixelMap;
        this->CreatePixelMap(pCaloHitList, pixelMap, vertex);
        if (pixelMap.size() < m_minHits)
            continue;

//        std::cout << "Writing pixel map with " << pixelMap.size() << " non-zero entries to the tree" << std::endl;
        const int view{static_cast<int>(pCaloHitList->front()->GetHitType())};
        for (auto const &pixel : pixelMap)
        {
            pixel_view.emplace_back(view);
            pixel_row.emplace_back(pixel.first.first);
            pixel_column.emplace_back(pixel.first.second);
            pixel_charge.emplace_back(pixel.second);
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
    TrackShowerCountingResults result;
    ViewToVertexPositionMap vertices;
    this->GetVertexPositions(vertices);

    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        VertexPosition &vertex(vertices.at(pCaloHitList->front()->GetHitType()));
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
        accessor[0][0] = vertices[pCaloHitList->front()->GetHitType()].GetDriftBin() / static_cast<float>(m_height);
        accessor[0][1] = vertices[pCaloHitList->front()->GetHitType()].GetWireBin() / static_cast<float>(m_width);

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
        const auto trkOutput = torch::softmax(outputList.get(0).toTensor(), 1);
        const auto shwOutput = torch::softmax(outputList.get(1).toTensor(), 1);
        
        result.AddScoresFromView(pCaloHitList->front()->GetHitType(), trkOutput, shwOutput);
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

    dummyEventPfoMetadata.m_propertiesToAdd["nTracksU"] = result.GetTrackClassPredictionFromView(TPC_VIEW_U);
    dummyEventPfoMetadata.m_propertiesToAdd["nTracksV"] = result.GetTrackClassPredictionFromView(TPC_VIEW_V);
    dummyEventPfoMetadata.m_propertiesToAdd["nTracksW"] = result.GetTrackClassPredictionFromView(TPC_VIEW_W);

    dummyEventPfoMetadata.m_propertiesToAdd["nShowersU"] = result.GetShowerClassPredictionFromView(TPC_VIEW_U);
    dummyEventPfoMetadata.m_propertiesToAdd["nShowersV"] = result.GetShowerClassPredictionFromView(TPC_VIEW_V);
    dummyEventPfoMetadata.m_propertiesToAdd["nShowersW"] = result.GetShowerClassPredictionFromView(TPC_VIEW_W);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pDummyEventPfo, dummyEventPfoMetadata));

    if (!pDummyEventPfoList->empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, dummyName, m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetVisibleParticles(const MCParticleList *const pMCParticleList, const MCParticle *const pMCNeutrino, const CaloHitList *const pCaloHitList, LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    if (m_useSimpleTruthLabels)
    {
        LArMCParticleHelper::PrimaryParameters parameters;
        parameters.m_minPrimaryGoodHits = m_goodMCTrackHits - 1; // No choice but to just use one of these. TODO: remove the simple option...
        parameters.m_minHitsForGoodView = m_goodMCTrackHits - 1;
        parameters.m_minPrimaryGoodViews = 1;
        parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();

        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToHitsMap);
    }
    else
    {
        // Find all MCParticles contributing > m_mcHitWeightThreshold to at least one CaloHit
        for (auto const *pCaloHit : *pCaloHitList)
        {
            for (auto const &mcHitWeightPair : pCaloHit->GetMCParticleWeightMap())
            {
                const MCParticle *pMCParticle{mcHitWeightPair.first};

                if (mcHitWeightPair.second > m_mcHitWeightThreshold)
                    mcToHitsMap[pMCParticle].emplace_back(pCaloHit);
            }
        }
    
        // We need to add any MC particles without hits into the map
        for (const MCParticle *const pMCParticle : *pMCParticleList)
        {
            if (!mcToHitsMap.count(pMCParticle))
                mcToHitsMap[pMCParticle] = CaloHitList();
        }
    }

    // Apply the number of hits cut to primaries and set failed primaries to be removed
    std::vector<const MCParticle*> keysToRemove;
    for (auto const &mcToHitsPair : mcToHitsMap)
    {
        const MCParticle *pMCParticle{mcToHitsPair.first};
        const unsigned int nHitsCut{this->IsTracklike(pMCParticle) ? m_goodMCTrackHits : m_goodMCShowerHits};
        //const MCParticle *pParentParticle{*(pMCParticle->GetParentList().begin())};

        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))
        {
              //std::cout << "Particle " << pMCParticle << " of type " << pMCParticle->GetParticleId() << " has " << mcToHitsPair.second.size() << " hits has parent " << pParentParticle << " of type " << pParentParticle->GetParticleId() << std::endl;
            if (mcToHitsPair.second.size() < nHitsCut)
                keysToRemove.emplace_back(pMCParticle);
        }
    }

    // Look for secondaries produced close to the vertex when their parent is not visible
    const HitType hitType{pCaloHitList->front()->GetHitType()};
    for (auto const &mcToHitsPair : mcToHitsMap)
    {
        const MCParticle *pMCParticle{mcToHitsPair.first};

        // Ignore the incoming neutrino
        if (pMCParticle->GetParentList().empty())
        {
            keysToRemove.emplace_back(pMCParticle);
            continue;
        }

        // Ignore primary particles this time
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))
            continue;

        const MCParticle *pParentParticle{*(pMCParticle->GetParentList().begin())};
        const unsigned int nHitsCut{this->IsTracklike(pMCParticle) ? m_goodMCTrackHits : m_goodMCShowerHits};
        // Look for secondaries produced by primaries that are considered not visible and see if they start close to the neutrino vertex position
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pParentParticle) && std::find(keysToRemove.begin(), keysToRemove.end(), pParentParticle) != keysToRemove.end())
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
//            else
//                std::cout << "Keeping particle of type " << pMCParticle->GetParticleId() << " with primary parent type " << pParentParticle->GetParticleId() << " that failed the selection criteria starts at a distance of " << distance << " from the neutrino interaction vertex" << std::endl;
        }
        else
            keysToRemove.emplace_back(pMCParticle);
    }

    for (auto const &key : keysToRemove)
        mcToHitsMap.erase(key);

//    for (auto const &mcToHitsPair : mcToHitsMap)
//    {
//        const MCParticle *pMCParticle{mcToHitsPair.first};
//        const MCParticle *pParentParticle{*(pMCParticle->GetParentList().begin())};
//        
//        const CartesianVector thisVtx(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCParticle->GetVertex(), hitType));
//        const CartesianVector thisEnd(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCParticle->GetEndpoint(), hitType));
//        const float particleLength{std::sqrt(thisVtx.GetDistanceSquared(thisEnd))};
//        std::cout << "Particle " << pMCParticle << " of type " << pMCParticle->GetParticleId() << " has " << mcToHitsPair.second.size() << " hits has parent " << pParentParticle << " of type " << pParentParticle->GetParticleId() << " and length " << particleLength << std::endl;
//    }

    // Calculate the angle between all remaining tracks and showers
    std::map<std::pair<const MCParticle *const, const MCParticle *const>, float> particleAngles;
    for (auto const &mcToHitsPair1 : mcToHitsMap)
    {
        for (auto const &mcToHitsPair2 : mcToHitsMap)
        {
            const MCParticle *const pMCParticle1{mcToHitsPair1.first};
            const MCParticle *const pMCParticle2{mcToHitsPair2.first};
            if (pMCParticle1 == pMCParticle2)
                continue;
            if (particleAngles.count(std::make_pair(pMCParticle2, pMCParticle1)))
                continue;
            const CartesianVector momentum1(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCParticle1->GetMomentum(), hitType)); 
            const CartesianVector momentum2(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCParticle2->GetMomentum(), hitType));
            const float openingAngle(momentum1.GetOpeningAngle(momentum2) * 180.f / M_PI);
            particleAngles[std::make_pair(mcToHitsPair1.first, mcToHitsPair2.first)] = openingAngle;
//            std::cout << "Angle between " << mcToHitsPair1.first << " and " << mcToHitsPair2.first << " = " << openingAngle << std::endl;
        }
    }

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
        float q{pixel.second};

        // Cap the charge at 10 ADC units and then rescale to be between 0 and 1.
        if (q > 10)
            q = 10;
        accessor[0][0][row][col] = q / 10.f;
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::CreatePixelMap(const CaloHitList *const pCaloHitList, PixelMap &pixelMap, VertexPosition &vertex) const
{
//    ViewToVertexPositionMap vertices;
//    if (m_useVertexForCrops)
//    this->GetVertexPositions(vertices);

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
//    std::cout << "Pixel Map Creation... do we need to crop? " << cropX << ", " << cropZ << std::endl;
//    CartesianVector vtx{CartesianVector(0.f, 0.f, 0.f)};
//    if (m_useVertexForCrops)

    float vtxPosX{vertex.GetPosition().GetX()};
    float vtxPosZ{vertex.GetPosition().GetZ()};

    // Adjust the min and max depending on reconstructed vertex position for events where
    // the reconstructed vertex isn't on a CaloHit, e.g. NC pizero events
    if (vtxPosX < xMin)
        xMin = vtxPosX;
    if (vtxPosX > xMax)
        xMax = vtxPosX;
    if (vtxPosZ < zMin)
        zMin = vtxPosZ;
    if (vtxPosZ > zMax)
        zMax = vtxPosZ;

    const float eps{0.0001f};
    if (cropX)
        this->GetCrop1D(pCaloHitList, vertex, xMin, xMax, xSpan, true);
    else
        vertex.AddVertexBin(true, static_cast<unsigned int>(std::floor(eps + (vtxPosX - xMin) / m_driftStep)));

    if (cropZ)
        this->GetCrop1D(pCaloHitList, vertex, zMin, zMax, zSpan, false);
    else
        vertex.AddVertexBin(false, static_cast<unsigned int>(std::floor(eps + (vtxPosZ - zMin) / wireStep)));


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
        Pixel pixel{driftBin, wireBin};
        if (!pixelMap.count(pixel))
            pixelMap[pixel] = pCaloHit->GetMipEquivalentEnergy();
        else
            pixelMap[pixel] += pCaloHit->GetMipEquivalentEnergy(); 
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetVertexPositions(ViewToVertexPositionMap &vertices) const
{
//    if (!m_useVertexForCrops)
//        return;

    const VertexList *pVertexList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    if (pVertexList->empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    if (pVertexList->front()->GetVertexType() != VERTEX_3D)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    vertices.emplace(TPC_VIEW_U, VertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_U)));
    vertices.emplace(TPC_VIEW_V, VertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_V)));
    vertices.emplace(TPC_VIEW_W, VertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_W)));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetCrop1D(const CaloHitList *const pCaloHitList, VertexPosition &vertexPosition,
    float &min, float &max, const float &span, const bool &isDrift) const
{
    // For some NC events we can expect the vertex to be outside the range from the hits
    const float vtxPos{isDrift ? vertexPosition.GetPosition().GetX() : vertexPosition.GetPosition().GetZ()};
    if (vtxPos < min)
        min = vtxPos;
    if (vtxPos > max)
        max = vtxPos;

    const unsigned int imageDimension{isDrift ? m_height : m_width};
    const float step{span / imageDimension};
    // Make a vector of bins from min to max using the required pitch
    const unsigned int nBins{static_cast<unsigned int>(std::ceil((max - min) / step))};
    if (nBins <= imageDimension)
        return;
    std::vector<float> dimensionBins(nBins, 0.f);
//    std::cout << "Min " << min << " and Max " << max << " with step " << step << " gives " << dimensionBins.size() << " bins" << std::endl;
    const float eps{0.0001f};
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const CartesianVector pos{pCaloHit->GetPositionVector()};
        const float dim{(isDrift ? pos.GetX() : pos.GetZ())};
        unsigned int bin{static_cast<unsigned int>(std::floor(eps + (dim - min) / step))};
        if (bin == nBins)
            bin = nBins - 1;
        if (bin >= nBins)
        {
            std::cout << " - Error! Bin " << bin << " outside range 0 to " << nBins - 1 << std::endl;
            std::cout << "  - Underlying values: " << dim - max << ", " << dim << ", " << min << ", " << max << ", " << step << std::endl;
            continue;
        }
        dimensionBins[bin] += pCaloHit->GetMipEquivalentEnergy();
    }

    const unsigned int vtxBin{static_cast<unsigned int>(std::floor(eps + (vtxPos - min) / step))};
    // Look for the required number of consecutive bins with the highest integrated charge
    std::map<unsigned int, float> startBinToCharge;
    //const float totalCharge{std::accumulate(dimensionBins.begin(), dimensionBins.end(), 0.f)};
    //std::cout << "Vertex: " << vtxPos << ", " << vtxBin << " :: " << nBins << std::endl;
    for (unsigned int i = 0; i < (nBins - imageDimension); ++i)
    {
        const bool vtxInRange{(vtxBin >= i) && (vtxBin < (i + imageDimension))};
        startBinToCharge[i] = (vtxInRange ? std::accumulate(dimensionBins.begin() + i, dimensionBins.begin() + i + imageDimension, 0.f) : 0.f);
    }

    // Get the maximum element from the map
    std::map<unsigned int, float>::iterator maxElementIter{std::max_element(
        startBinToCharge.begin(), startBinToCharge.end(), [](const auto &a, const auto &b){return a.second < b.second;})};
    const unsigned int maxBin{maxElementIter->first};
//    std::cout << maxBin << ", " << maxElementIter->second << " :: " << maxBinLoop << ", " << maxValue << std::endl;
    const float localMin{min + maxBin * step};
    const float localMax{localMin + imageDimension * step};

    vertexPosition.AddVertexBin(isDrift, vtxBin - maxBin);

//    std::cout << "Crop " << isDrift << ": " << min << ", " << vtxPos << ", " << vtxBin << ", " << maxBin << ", " << max << std::endl;

    min = localMin;
    max = localMax;

	  if (vtxPos < localMin || vtxPos > localMax)
	      std::cout << "CNNTrackShowerCountingAlgorithm::GetCrop1D: reconstructed vertex outside cropped region! " << vtxPos << ", " << localMin << ", " << localMax << std::endl;
//    std::cout << "Best bin " << maxBin << " with charge sum " << startBinToCharge[maxBin] << " gives minimum " << min << " and max " << max << std::endl;
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

void CNNTrackShowerCountingAlgorithm::CountMCPrimaries(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, unsigned int &nTracks, unsigned int &nShowers) const
{
    nTracks = 0;
    nShowers = 0;
    for (const auto &[mc, hits] : mcToHitsMap)
    {
        if (this->IsShowerlike(mc))
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

    
//    PANDORA_RETURN_RESULT_IF_AND_IF(
//        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseVertexForCrops", m_useVertexForCrops));
//    if (m_useVertexForCrops)
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseSimpleTruthLabels", m_useSimpleTruthLabels));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GoodMCTrackHits", m_goodMCTrackHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GoodMCShowerHits", m_goodMCShowerHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCHitWeightThreshold", m_mcHitWeightThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SecondaryDistanceThreshold", m_secondaryDistanceThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHits", m_minHits));

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

void CNNTrackShowerCountingAlgorithm::VertexPosition::AddVertexBin(const bool &isDrift, const unsigned int &bin)
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

void CNNTrackShowerCountingAlgorithm::TrackShowerCountingResults::AddScoresFromView(const HitType &view, const torch::Tensor &trackScores, const torch::Tensor &showerScores)
{
    if (m_trackScores.count(view) || m_showerScores.count(view))
    {
        std::cerr << "TrackShowerCountingResults::AddScoresFromView: scores for view " << view << " already exist. Doing nothing." << std::endl;
        return;
    }

    // The tensors have a batch dimension so are technically 2D
    auto trackAccessor = trackScores.accessor<float, 2>();
    std::vector<float> trkScoreVector;
    for (unsigned int element = 0; element < trackScores.size(1); ++element)
        trkScoreVector.emplace_back(trackAccessor[0][element]);
    m_trackScores[view] = trkScoreVector;

    auto showerAccessor = showerScores.accessor<float, 2>();
    std::vector<float> shwScoreVector;
    for (unsigned int element = 0; element < showerScores.size(1); ++element)
        shwScoreVector.emplace_back(showerAccessor[0][element]);
    m_showerScores[view] = shwScoreVector;
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

}

