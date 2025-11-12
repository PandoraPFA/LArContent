/**
 *  @file   larpandoradlcontent/LArVertex/DlVertexingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning vertexing algorithm.
 *
 *  $Log: $
 */

#include <chrono>
#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoradlcontent/LArHelpers/LArCanvasHelper.h"

#include "larpandoradlcontent/LArVertex/DlVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlVertexingAlgorithm::DlVertexingAlgorithm() :
    m_event{-1},
    m_visualise{false},
    m_writeTree{false},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
{
}

DlVertexingAlgorithm::~DlVertexingAlgorithm()
{
    if (m_writeTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "VertexAssessmentAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlVertexingAlgorithm::PrepareTrainingSample()
{
    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(pCaloHitList2D, pMCParticleList, mcToHitsMap);
    MCParticleList hierarchy;
    LArMCParticleHelper::CompleteMCHierarchy(mcToHitsMap, hierarchy);

    // Get boundaries for hits and make x dimension common
    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!(isU || isV || isW))
            return STATUS_CODE_NOT_ALLOWED;

        CartesianPointVector vertices;
        for (const MCParticle *mc : hierarchy)
        {
            if (LArMCParticleHelper::IsNeutrino(mc))
                vertices.push_back(mc->GetVertex());
        }
        if (vertices.empty())
            continue;
        const CartesianVector &vertex{vertices.front()};
        const std::string trainingFilename{m_trainingOutputFile + "_" + listname + ".csv"};
        const unsigned long nVertices{1};
        unsigned long nHits{0};
        const unsigned int nuance{LArMCParticleHelper::GetNuanceCode(hierarchy.front())};

        // Vertices
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const double xVtx{vertex.GetX()};
        double zVtx{0.};
        if (isW)
            zVtx = transform->YZtoW(vertex.GetY(), vertex.GetZ());
        else if (isV)
            zVtx = transform->YZtoV(vertex.GetY(), vertex.GetZ());
        else
            zVtx = transform->YZtoU(vertex.GetY(), vertex.GetZ());

        // Calo hits
        double xMin{driftMin}, xMax{driftMax}, zMin{wireMin[view]}, zMax{wireMax[view]};

        // Only train on events where the vertex resides within the image - with a small tolerance
        if (!(xVtx > (xMin - 1.f) && xVtx < (xMax + 1.f) && zVtx > (zMin - 1.f) && zVtx < (zMax + 1.f)))
            continue;

        LArMvaHelper::MvaFeatureVector featureVector;
        featureVector.emplace_back(static_cast<double>(nuance));
        featureVector.emplace_back(static_cast<double>(nVertices));
        featureVector.emplace_back(xVtx);
        featureVector.emplace_back(zVtx);
        // Retain the hit region
        featureVector.emplace_back(xMin);
        featureVector.emplace_back(xMax);
        featureVector.emplace_back(zMin);
        featureVector.emplace_back(zMax);

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()}, adc{pCaloHit->GetMipEquivalentEnergy()};
            // If on a refinement pass, drop hits outside the region of interest
            if (m_pass > 1 && (x < xMin || x > xMax || z < zMin || z > zMax))
                continue;
            featureVector.emplace_back(static_cast<double>(x));
            featureVector.emplace_back(static_cast<double>(z));
            featureVector.emplace_back(static_cast<double>(adc));
            ++nHits;
        }
        featureVector.insert(featureVector.begin() + 8, static_cast<double>(nHits));
        // Only write out the feature vector if there were enough hits in the region of interest
        if (nHits > 10)
            LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::Infer()
{
    if (m_pass == 1)
    {
        ++m_event;
    }
    else
    {
        // INFO: Check if there is a zoom in region for the second pass.
        const VertexList *pVertexList(nullptr);
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_inputVertexListName, pVertexList) || pVertexList == nullptr ||
            pVertexList->empty())
        {
            std::cout << "DLVertexing: Input vertex list is empty! Can't perform pass " << m_pass << std::endl;
            return STATUS_CODE_SUCCESS;
        }
    }

    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    CartesianPointVector vertexCandidatesU, vertexCandidatesV, vertexCandidatesW;
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        HitType view{pCaloHitList->front()->GetHitType()};
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

        LArDLHelper::TorchInput input;
        PixelVector pixelVector;
        this->MakeNetworkInputFromHits(*pCaloHitList, view, driftMin, driftMax, wireMin[view], wireMax[view], input, pixelVector);

        // Run the input through the trained model
        LArDLHelper::TorchInputVector inputs;
        inputs.push_back(input);
        LArDLHelper::TorchOutput output;
        if (isU)
            LArDLHelper::Forward(m_modelU, inputs, output);
        else if (isV)
            LArDLHelper::Forward(m_modelV, inputs, output);
        else
            LArDLHelper::Forward(m_modelW, inputs, output);

        int colOffset{0}, rowOffset{0}, canvasWidth{m_width}, canvasHeight{m_height};
        this->GetCanvasParameters(output, pixelVector, colOffset, rowOffset, canvasWidth, canvasHeight);

        float **canvas{new float *[canvasHeight]};
        for (int row = 0; row < canvasHeight; ++row)
            canvas[row] = new float[canvasWidth]{};

        // we want the maximum value in the num_classes dimension (1) for every pixel
        auto classes{torch::argmax(output, 1)};
        // the argmax result is a 1 x height x width tensor where each element is a class id
        auto classesAccessor{classes.accessor<int64_t, 3>()};
        const double scaleFactor{std::sqrt(m_height * m_height + m_width * m_width)};
        std::map<int, bool> haveSeenMap;
        for (const auto &[row, col] : pixelVector)
        {
            const auto cls{classesAccessor[0][row][col]};
            if (cls > 0 && cls < m_nClasses)
            {
                const int inner{static_cast<int>(std::round(std::ceil(scaleFactor * m_thresholds[cls - 1])))};
                const int outer{static_cast<int>(std::round(std::ceil(scaleFactor * m_thresholds[cls])))};
                LArCanvasHelper::DrawRing(canvas, row + rowOffset, col + colOffset, inner, outer, 1.f / (outer * outer - inner * inner));
            }
        }

        CartesianPointVector positionVector;
        this->MakeWirePlaneCoordinatesFromCanvas(
            canvas, canvasWidth, canvasHeight, colOffset, rowOffset, view, driftMin, driftMax, wireMin[view], wireMax[view], positionVector);
        if (isU)
            vertexCandidatesU.emplace_back(positionVector.front());
        else if (isV)
            vertexCandidatesV.emplace_back(positionVector.front());
        else
            vertexCandidatesW.emplace_back(positionVector.front());

#ifdef MONITORING
        if (m_visualise)
        {
            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            const MCParticleList *pMCParticleList{nullptr};
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

            CartesianVector trueVertex3D(0, 0, 0);
            if (LArMCParticleHelper::GetTrueVertex(pMCParticleList, trueVertex3D))
            {
                float x{0.f}, u{0.f}, v{0.f}, w{0.f};
                const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
                LArVertexHelper::GetTrueVertexPosition(trueVertex3D, transform, x, u, v, w);
                if (isU)
                {
                    const CartesianVector trueVertex(x, 0.f, u);
                    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueVertex, "U(true)", BLUE, 3));
                }
                else if (isV)
                {
                    const CartesianVector trueVertex(x, 0.f, v);
                    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueVertex, "V(true)", BLUE, 3));
                }
                else if (isW)
                {
                    const CartesianVector trueVertex(x, 0.f, w);
                    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueVertex, "W(true)", BLUE, 3));
                }
                for (const auto &pos : positionVector)
                {
                    std::string label{isU ? "U" : isV ? "V" : "W"};
                    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, label, RED, 3));
                }
            }
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
#endif

        for (int row = 0; row < canvasHeight; ++row)
            delete[] canvas[row];
        delete[] canvas;
    }

    int nEmptyLists{0};
    if (vertexCandidatesU.empty())
        ++nEmptyLists;
    if (vertexCandidatesV.empty())
        ++nEmptyLists;
    if (vertexCandidatesW.empty())
        ++nEmptyLists;
    std::vector<VertexTuple> vertexTuples;
    CartesianPointVector candidates3D;
    if (nEmptyLists == 0)
    {
        vertexTuples.emplace_back(VertexTuple(this->GetPandora(), vertexCandidatesU.front(), vertexCandidatesV.front(), vertexCandidatesW.front()));
    }
    else if (nEmptyLists == 1)
    {
        if (vertexCandidatesU.empty())
        { // V and W available
            vertexTuples.emplace_back(VertexTuple(this->GetPandora(), vertexCandidatesV.front(), vertexCandidatesW.front(), TPC_VIEW_V, TPC_VIEW_W));
        }
        else if (vertexCandidatesV.empty())
        { // U and W available
            vertexTuples.emplace_back(VertexTuple(this->GetPandora(), vertexCandidatesU.front(), vertexCandidatesW.front(), TPC_VIEW_U, TPC_VIEW_W));
        }
        else
        { // U and V available
            vertexTuples.emplace_back(VertexTuple(this->GetPandora(), vertexCandidatesU.front(), vertexCandidatesV.front(), TPC_VIEW_U, TPC_VIEW_V));
        }
    }
    else
    { // Not enough views to reconstruct a 3D vertex
        std::cout << "Insufficient 2D vertices to reconstruct a 3D vertex" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    if (m_visualise)
    {
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexTuples.front().GetPosition(), "candidate", GREEN, 1));
    }

    if (!vertexTuples.empty())
    {
        const CartesianVector &vertex{vertexTuples.front().GetPosition()};
        CartesianPointVector vertexCandidates;
        vertexCandidates.emplace_back(vertex);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->MakeCandidateVertexList(vertexCandidates));
    }
    else
    {
        std::cout << "Insufficient 2D vertices to reconstruct a 3D vertex" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeNetworkInputFromHits(const CaloHitList &caloHits, const HitType view, const float xMin,
    const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());

    // Determine the bin edges
    std::vector<double> xBinEdges(m_width + 1);
    std::vector<double> zBinEdges(m_height + 1);
    xBinEdges[0] = xMin - 0.5f * m_driftStep;
    const double dx = ((xMax + 0.5f * m_driftStep) - xBinEdges[0]) / m_width;
    for (int i = 1; i < m_width + 1; ++i)
        xBinEdges[i] = xBinEdges[i - 1] + dx;
    zBinEdges[0] = zMin - 0.5f * pitch;
    const double dz = ((zMax + 0.5f * pitch) - zBinEdges[0]) / m_height;
    for (int i = 1; i < m_height + 1; ++i)
        zBinEdges[i] = zBinEdges[i - 1] + dz;

    LArDLHelper::InitialiseInput({1, 1, m_height, m_width}, networkInput);
    auto accessor = networkInput.accessor<float, 4>();

    for (const CaloHit *pCaloHit : caloHits)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        if (m_pass > 1)
        {
            if (x < xMin || x > xMax || z < zMin || z > zMax)
                continue;
        }
        const float adc{pCaloHit->GetMipEquivalentEnergy()};
        const int pixelX{static_cast<int>(std::floor((x - xBinEdges[0]) / dx))};
        const int pixelZ{static_cast<int>(std::floor((z - zBinEdges[0]) / dz))};
        accessor[0][0][pixelZ][pixelX] += adc;
    }
    for (int row = 0; row < m_height; ++row)
    {
        for (int col = 0; col < m_width; ++col)
        {
            const float value{accessor[0][0][row][col]};
            if (value > 0)
                pixelVector.emplace_back(std::make_pair(row, col));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeWirePlaneCoordinatesFromCanvas(float **canvas, const int canvasWidth, const int canvasHeight,
    const int columnOffset, const int rowOffset, const HitType view, const float xMin, const float xMax, const float zMin, const float zMax,
    CartesianPointVector &positionVector) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());

    const double dx = ((xMax + 0.5f * m_driftStep) - (xMin - 0.5f * m_driftStep)) / m_width;
    const double dz = ((zMax + 0.5f * pitch) - (zMin - 0.5f * pitch)) / m_height;

    float best{-1.f};
    int rowBest{0}, colBest{0};
    for (int row = 0; row < canvasHeight; ++row)
        for (int col = 0; col < canvasWidth; ++col)
            if (canvas[row][col] > 0 && canvas[row][col] > best)
            {
                best = canvas[row][col];
                rowBest = row;
                colBest = col;
            }

    const float x{static_cast<float>((colBest - columnOffset) * dx + xMin)};
    const float z{static_cast<float>((rowBest - rowOffset) * dz + zMin)};

    CartesianVector pt(x, 0.f, z);
    positionVector.emplace_back(pt);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeCandidateVertexList(const CartesianPointVector &positions)
{
    const VertexList *pVertexList{nullptr};
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (const CartesianVector &position : positions)
    {
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pVertex(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));

    return STATUS_CODE_SUCCESS;
}

#ifdef MONITORING
void DlVertexingAlgorithm::PopulateRootTree(const std::vector<VertexTuple> &vertexTuples, const pandora::CartesianPointVector &vertexCandidatesU,
    const pandora::CartesianPointVector &vertexCandidatesV, const pandora::CartesianPointVector &vertexCandidatesW) const
{
    if (m_writeTree)
    {
        const MCParticleList *pMCParticleList{nullptr};
        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetCurrentList(*this, pMCParticleList))
        {
            if (pMCParticleList)
            {
                LArMCParticleHelper::MCContributionMap mcToHitsMap;
                MCParticleVector primaries;
                LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
                if (!primaries.empty())
                {
                    const MCParticle *primary{primaries.front()};
                    const MCParticleList &parents{primary->GetParentList()};
                    if (parents.size() == 1)
                    {
                        const MCParticle *trueNeutrino{parents.front()};
                        if (LArMCParticleHelper::IsNeutrino(trueNeutrino))
                        {
                            const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
                            const CartesianVector &trueVertex{primaries.front()->GetVertex()};
                            if (LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, m_volumeType))
                            {
                                const CartesianVector &recoVertex{vertexTuples.front().GetPosition()};
                                const float tx{trueVertex.GetX()};
                                const float tu{static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ()))};
                                const float tv{static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ()))};
                                const float tw{static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ()))};
                                const float rx_u{vertexCandidatesU.front().GetX()};
                                const float ru{vertexCandidatesU.front().GetZ()};
                                const float rx_v{vertexCandidatesV.front().GetX()};
                                const float rv{vertexCandidatesV.front().GetZ()};
                                const float rx_w{vertexCandidatesW.front().GetX()};
                                const float rw{vertexCandidatesW.front().GetZ()};
                                const float dr_u{std::sqrt((rx_u - tx) * (rx_u - tx) + (ru - tu) * (ru - tu))};
                                const float dr_v{std::sqrt((rx_v - tx) * (rx_v - tx) + (rv - tv) * (rv - tv))};
                                const float dr_w{std::sqrt((rx_w - tx) * (rx_w - tx) + (rw - tw) * (rw - tw))};
                                const CartesianVector &dv{recoVertex - trueVertex};
                                const float dr{dv.GetMagnitude()};
                                const float dx{dv.GetX()}, dy{dv.GetY()}, dz{dv.GetZ()};
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "event", m_event));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "pass", m_pass));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dr_u", dr_u));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dr_v", dr_v));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dr_w", dr_w));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dr", dr));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dx", dx));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dy", dy));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dz", dz));
                                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, DlVertexingBaseAlgorithm::ReadSettings(xmlHandle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));

    if (!m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
        if (m_writeTree)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
