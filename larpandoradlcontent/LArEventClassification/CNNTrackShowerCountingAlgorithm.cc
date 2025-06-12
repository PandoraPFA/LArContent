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
    m_validationMode{false},
    m_trainingTreeName{""},
    m_trainingFileName{""},
    m_height{512},
    m_width{512},
    m_wiresPerPixel{1},
    m_driftStep{0.5f},
    m_useVertexForCrops{true},
    m_goodMCPrimaryHits{5},
    m_minHits{10}
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
    StatusCode result{STATUS_CODE_FAILURE};

    // Slightly more complex that naively expected due to validation requiring both truth and predictions
    if (m_trainingMode || m_validationMode)
    {
        result = this->PrepareTrainingSample();
        if (STATUS_CODE_FAILURE == result)
            return result;
    }

    if (!m_trainingMode)
        result = this->Infer();

    return result;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CNNTrackShowerCountingAlgorithm::PrepareTrainingSample()
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    // We need to define the number of tracks and showers per view, so we need to do this per view
    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = m_goodMCPrimaryHits;
    parameters.m_minHitsForGoodView = m_goodMCPrimaryHits;
    parameters.m_minPrimaryGoodViews = 1;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();

    bool isCC{false};
    int nuPDG{0};
    std::map<HitType, unsigned int> nTracksPerView;
    std::map<HitType, unsigned int> nShowersPerView;
    float nuEnergy{0.f};
    CartesianVector nuVertex{0.f, 0.f, 0.f};
    unsigned int emptyViews{0};
    for (const std::string &name : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, name, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        LArMCParticleHelper::MCContributionMap mcToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(
            pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToHitsMap);
//        std::cout << " - Got " << mcToHitsMap.size() << " contributing MC particles of " << pMCParticleList->size() << " for hit list " << name << " with " << pCaloHitList->size() << " hits" << std::endl;

        if (mcToHitsMap.empty())
        {
            ++emptyViews;
            continue;
        }

        MCParticleList hierarchy;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CompleteMCHierarchy(mcToHitsMap, hierarchy));
        MCParticleList primaries;

        // Only fill the event level information once
        if (nuPDG == 0)
        {
            MCParticleList mcPrimaries;
            this->GetMCPrimaries(hierarchy, mcPrimaries);
            const InteractionDescriptor &intType{LArInteractionTypeHelper::GetInteractionDescriptor(mcPrimaries)};
            isCC = intType.IsCC();
            const MCParticle *const pMCNeutrino{hierarchy.front()};
            if (!LArMCParticleHelper::IsNeutrino(pMCNeutrino))
                return STATUS_CODE_FAILURE;
            nuPDG = pMCNeutrino->GetParticleId();
            nuEnergy = pMCNeutrino->GetEnergy();
            nuVertex = pMCNeutrino->GetVertex();
//            std::cout << " - Got " << (isCC ? "CC " : "NC ") << nuPDG << " event" << std::endl;
        }

        const HitType hitType{pCaloHitList->front()->GetHitType()};
        unsigned int nTracks{0}, nShowers{0};
        this->CountMCPrimaries(mcToHitsMap, nTracks, nShowers);
//        std::cout << " - Got " << nTracks << " tracks and " << nShowers << " showers in view " << static_cast<int>(hitType) << std::endl;
        nTracksPerView[hitType] = nTracks;
        nShowersPerView[hitType] = nShowers;
    }

    if (m_validationMode)
    {
        // Print the true values and stop
        std::cout << "View: " << TPC_VIEW_U << ", true number of tracks =  " << nTracksPerView.at(TPC_VIEW_U) << " and true number of showers = " << nShowersPerView.at(TPC_VIEW_U) << std::endl;
        std::cout << "View: " << TPC_VIEW_V << ", true number of tracks =  " << nTracksPerView.at(TPC_VIEW_V) << " and true number of showers = " << nShowersPerView.at(TPC_VIEW_V) << std::endl;
        std::cout << "View: " << TPC_VIEW_W << ", true number of tracks =  " << nTracksPerView.at(TPC_VIEW_W) << " and true number of showers = " << nShowersPerView.at(TPC_VIEW_W) << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    if (emptyViews)
    {
        std::cout << "Skipping event as it does not have enough hits or associated primary particles to make a training sample" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    const float nuVertexX{nuVertex.GetX()};
    const float nuVertexY{nuVertex.GetY()};
    const float nuVertexZ{nuVertex.GetZ()};

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuPDG", nuPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "isCC", static_cast<int>(isCC)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuEnergy", nuEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuVertexX", nuVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuVertexY", nuVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "nuVertexZ", nuVertexZ));
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
        PixelMap pixelMap;
        this->CreatePixelMap(pCaloHitList, pixelMap);
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

    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        PixelMap pixelMap;
        this->CreatePixelMap(pCaloHitList, pixelMap);
        if (pixelMap.size() < m_minHits)
            continue;

        LArDLHelper::TorchInput input;
        this->MakeNetworkInputFromPixelMap(pixelMap, input);
        LArDLHelper::TorchInputVector inputs;
        inputs.emplace_back(input);
        LArDLHelper::TorchMultiOutput output;
        LArDLHelper::Forward(m_model, inputs, output);

        // Network has two outputs that are returned in Python as a list. Here that will be a tuple containing
        // two tensors, one for each of the outputs.
        const auto outputList = output.toList();

        // Extract the tensors from the tuple
        // Due to the way PyTorch includes the softmax in the loss function we need to apply it here if we
        // want out scores to sum to one
//        const auto trkOutput = torch::softmax(outputTuple->elements()[0].toTensor(), 1);
//        const auto shwOutput = torch::softmax(outputTuple->elements()[1].toTensor(), 1);
        const auto trkOutput = torch::softmax(outputList.get(0).toTensor(), 1);
        const auto shwOutput = torch::softmax(outputList.get(1).toTensor(), 1);
        
        result.AddScoresFromView(pCaloHitList->front()->GetHitType(), trkOutput, shwOutput);
    }

    return STATUS_CODE_SUCCESS;
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

void CNNTrackShowerCountingAlgorithm::CreatePixelMap(const CaloHitList *const pCaloHitList, PixelMap &pixelMap) const
{
    ViewToVertexPositionMap vertices;
    if (m_useVertexForCrops)
        this->GetVertexPositions(vertices);

    float xMin{std::numeric_limits<float>::max()}, xMax{-1.f * std::numeric_limits<float>::max()};
    float zMin{std::numeric_limits<float>::max()}, zMax{-1.f * std::numeric_limits<float>::max()};
    this->GetHitRegion(pCaloHitList, xMin, xMax, zMin, zMax);

    // Get the wire pitch for this view
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view{pCaloHitList->front()->GetHitType()};
    const float pitch{view == TPC_VIEW_U ? pTPC->GetWirePitchU() : (view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW())};

    const float zSpan{pitch * m_wiresPerPixel * m_width};
    const float xSpan{m_driftStep * m_height};
    // If the min / max values exceed the above spans then we need to crop
    const bool cropZ{(zMax - zMin) > zSpan ? true : false};
    const bool cropX{(xMax - xMin) > xSpan ? true : false};
//    std::cout << "Pixel Map Creation... do we need to crop? " << cropX << ", " << cropZ << std::endl;
    CartesianVector vtx{CartesianVector(0.f, 0.f, 0.f)};
    if (m_useVertexForCrops)
        vtx = vertices.at(view);
    if (cropX)
        this->GetCrop1D(pCaloHitList, vtx, xMin, xMax, xSpan, true);
    if (cropZ)
        this->GetCrop1D(pCaloHitList, vtx, zMin, zMax, zSpan, false);

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const float driftPos{pCaloHit->GetPositionVector().GetX()};
        if (driftPos < xMin || driftPos > xMax)
            continue;
        const float wirePos{pCaloHit->GetPositionVector().GetZ()};
        if (wirePos < zMin || wirePos > zMax)
            continue;

        const float eps{0.0001f};
        const unsigned int driftBin{static_cast<unsigned int>(std::floor(eps + (driftPos - xMin) / m_driftStep))};
        const float wireStep{pitch * m_wiresPerPixel};
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
    if (!m_useVertexForCrops)
        return;

    const VertexList *pVertexList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    if (pVertexList->empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    if (pVertexList->front()->GetVertexType() != VERTEX_3D)
        return;

    vertices.emplace(TPC_VIEW_U, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_U));
    vertices.emplace(TPC_VIEW_V, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_V));
    vertices.emplace(TPC_VIEW_W, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_W));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetCrop1D(const CaloHitList *const pCaloHitList, const CartesianVector &vertexPosition,
    float &min, float &max, const float &span, const bool &isDrift) const
{
    // For some NC events we can expect the vertex to be outside the range from the hits
    const float vtxPos{m_useVertexForCrops ? (isDrift ? vertexPosition.GetX() : vertexPosition.GetZ()) : min};
    if (m_useVertexForCrops)
    {
        if (vtxPos < min)
            min = vtxPos;
        if (vtxPos > max)
            max = vtxPos;
    }

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
        const unsigned int bin{static_cast<unsigned int>(std::floor(eps + (dim - min) / step))};
        if (bin > dimensionBins.size() - 1)
        {
            std::cout << " - Error! Bin " << bin << " outside range 0 to " << dimensionBins.size() - 1 << std::endl;
            std::cout << "  - Underlying values: " << dim << ", " << min << ", " << max << ", " << step << std::endl;
            continue;
        }
        dimensionBins[bin] += pCaloHit->GetMipEquivalentEnergy();
    }

    const unsigned int vtxBin{static_cast<unsigned int>(std::floor(eps + (vtxPos - min) / step)) - 1};

    // Look for the required number of consecutive bins with the highest integrated charge
    std::map<unsigned int, float> startBinToCharge;
    const float totalCharge{std::accumulate(dimensionBins.begin(), dimensionBins.end(), 0.f)};
    //std::cout << "Vertex: " << vtxPos << ", " << vtxBin << " :: " << nBins << std::endl;
    for (unsigned int i = 0; i < (nBins - imageDimension); ++i)
    {
        startBinToCharge[i] = std::accumulate(dimensionBins.begin() + i, dimensionBins.begin() + i + imageDimension, 0.f);
        // If we are using the reco vertex then set the sum to zero if it doesn't include the vertex
        if (m_useVertexForCrops)
        {
            if (vtxBin < i || vtxBin > (i + imageDimension - 1))
                startBinToCharge[i] = 0;
        }
        else
        {
            // Give some upstream bias
            startBinToCharge[i] -= (totalCharge / dimensionBins.size()) * i;
        }
    }

    // Get the maximum element from the map
    std::map<unsigned int, float>::iterator maxElementIter{std::max_element(
        startBinToCharge.begin(), startBinToCharge.end(), [](const auto &a, const auto &b){return a.second < b.second;})};
    const unsigned int maxBin{maxElementIter->first};
//    std::cout << maxBin << ", " << maxElementIter->second << " :: " << maxBinLoop << ", " << maxValue << std::endl;
    const float localMin{min + maxBin * step};
    const float localMax{localMin + imageDimension * step};
    if (m_useVertexForCrops)
    {
        const float vtxPos{isDrift ? vertexPosition.GetX() : vertexPosition.GetZ()}; 
	      if (vtxPos < localMin || vtxPos > localMax)
	          std::cout << "CNNTrackShowerCountingAlgorithm::GetCrop1D: reconstructed vertex outside cropped region! " << vtxPos << ", " << localMin << ", " << localMax << std::endl;
    }

    min = localMin;
    max = localMax;

//    std::cout << "Best bin " << maxBin << " with charge sum " << startBinToCharge[maxBin] << " gives minimum " << min << " and max " << max << std::endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CNNTrackShowerCountingAlgorithm::CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCParticleList &mcHierarchy) const
{
    try
    {
        for (const auto &[mc, hits] : mcToHitsMap)
        {
            (void)hits;
            mcHierarchy.push_back(mc);
            LArMCParticleHelper::GetAllAncestorMCParticles(mc, mcHierarchy);
        }
    }
    catch (const StatusCodeException &e)
    {
        return e.GetStatusCode();
    }

    // Move the neutrino to the front of the list
    auto pivot =
        std::find_if(mcHierarchy.begin(), mcHierarchy.end(), [](const MCParticle *mc) -> bool { return LArMCParticleHelper::IsNeutrino(mc); });
    (void)pivot;
    if (pivot != mcHierarchy.end())
        std::rotate(mcHierarchy.begin(), pivot, std::next(pivot));
    else
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::GetMCPrimaries(const MCParticleList &mcHierarchy, MCParticleList &mcPrimaries) const
{
    for (const MCParticle *const particle : mcHierarchy)
    {
        if (LArMCParticleHelper::IsPrimary(particle))
            mcPrimaries.push_back(particle);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CNNTrackShowerCountingAlgorithm::CountMCPrimaries(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, unsigned int &nTracks, unsigned int &nShowers) const
{
    nTracks = 0;
    nShowers = 0;
    for (const auto &[mc, hits] : mcToHitsMap)
    {
        // Check if it made some hits or not
        if (hits.size() < m_goodMCPrimaryHits)
            continue;

        const int pdg{mc->GetParticleId()};
        const bool isShw{(std::abs(pdg) != 11 && pdg != 22) ? false : true};

//        std::cout << "Found Good MC Particle of type " << pdg << " (" << isShw << ") with " << hits.size() << " hits" << std::endl;

        if (isShw)
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidationMode", m_validationMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_height));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_width));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WiresPerPixel", m_wiresPerPixel));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DriftStep", m_driftStep));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseVertexForCrops", m_useVertexForCrops));
    if (m_useVertexForCrops)
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GoodMCPrimaryHits", m_goodMCPrimaryHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHits", m_minHits));

    return STATUS_CODE_SUCCESS;
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

    const size_t nTrkPred(std::distance(trkScoreVector.begin(), std::max_element(trkScoreVector.begin(), trkScoreVector.end())));
    const size_t nShwPred(std::distance(shwScoreVector.begin(), std::max_element(shwScoreVector.begin(), shwScoreVector.end())));

    std::cout << "View: " << view << ", predicted number of tracks =  " << nTrkPred << " and prediced number of showers = " << nShwPred << std::endl;
}

}

