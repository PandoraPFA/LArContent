/**
 *  @file   larpandoradlcontent/LArVertex/DlVertexingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning vertexing algorithm.
 *
 *  $Log: $
 */

#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoradlcontent/LArVertex/DlVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlVertexingAlgorithm::DlVertexingAlgorithm():
    m_trainingMode{false},
    m_trainingOutputFile{""},
    m_height{256},
    m_width{256},
    m_pixelShift{0.f},
    m_pixelScale{1.f},
    m_visualise{false}
{
}

DlVertexingAlgorithm::~DlVertexingAlgorithm()
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::Run()
{
    if (m_trainingMode)
        return this->Train();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlVertexingAlgorithm::Train()
{
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMCToHitsMap(mcToHitsMap));
    MCParticleList hierarchy;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CompleteMCHierarchy(mcToHitsMap, hierarchy));

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
            if (!LArMCParticleHelper::IsNeutrino(mc))
            {   // Secondary vertices - will have coicident start and endpoints, only add once
/*                if (LArMCParticleHelper::GetHierarchyTier(mc) == 1)
                {
                    const int pdg{std::abs(mc->GetParticleId())};
                    if (pdg == 11 || pdg == 22)
                    {   // Shower, don't retain endpoint
                        CartesianVector vertex{mc->GetVertex()};
                        if (std::find(vertices.begin(), vertices.end(), vertex) == vertices.end())
                            vertices.push_back(vertex);
                    }
                    else
                    {   // Track, retain endpoint
                        const CartesianVector &vertex{mc->GetVertex()};
                        const CartesianVector &endpoint{mc->GetEndpoint()};
                        if (std::find(vertices.begin(), vertices.end(), vertex) == vertices.end())
                            vertices.push_back(vertex);
                        if (std::find(vertices.begin(), vertices.end(), endpoint) == vertices.end())
                            vertices.push_back(endpoint);
                    }
                }*/
            }
            else
            {   // Primary vertex - unique, so add directly to vertex vector
                vertices.push_back(mc->GetVertex());
            }
        }

        const std::string trainingFilename{m_trainingOutputFile + "_" + listname + ".csv"};
        const unsigned long nVertices{vertices.size()};
        unsigned long nHits{0};
        const unsigned int nuance{LArMCParticleHelper::GetNuanceCode(hierarchy.front())};
        LArMvaHelper::MvaFeatureVector featureVector;

        featureVector.push_back(static_cast<double>(nuance));
        featureVector.push_back(static_cast<double>(nVertices));

        // Vertices
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        for (const CartesianVector vertex : vertices)
        {
            const double x{vertex.GetX()};
            featureVector.push_back(x);
            if (isW)
                featureVector.push_back(transform->YZtoW(vertex.GetY(), vertex.GetZ()));
            else if (isV)
                featureVector.push_back(transform->YZtoV(vertex.GetY(), vertex.GetZ()));
            else
                featureVector.push_back(transform->YZtoU(vertex.GetY(), vertex.GetZ()));
        }

        // Calo hits
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            try
            {
                const MCParticle *const pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                // Discard non-reconstructable hits
                if (mcToHitsMap.find(pMCParticle) == mcToHitsMap.end())
                    continue;

                const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()}, adc{pCaloHit->GetInputEnergy()};
                featureVector.push_back(static_cast<double>(x));
                featureVector.push_back(static_cast<double>(z));
                featureVector.push_back(static_cast<double>(adc));
                ++nHits;
            }
            catch (...)
            {
            }
        }
        featureVector.insert(featureVector.begin() + 2, static_cast<double>(nHits));
        LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::Infer()
{
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    CartesianPointVector vertexCandidatesU, vertexCandidatesV, vertexCandidatesW;
    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        HitType view{pCaloHitList->front()->GetHitType()};
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

        LArDLHelper::TorchInput input;
        PixelVector pixelVector;
        this->MakeNetworkInputFromHits(*pCaloHitList, input, pixelVector);

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

        double thresholds[]{0., 0.00275, 0.00825, 0.01925, 0.03575, 0.05775, 0.08525, 0.12375, 0.15125, 0.20625, 0.26125, 0.31625, 0.37125,
            0.42625, 0.50875, 0.59125, 0.67375, 0.75625, 0.85, 1.0};
        int colOffset{0}, rowOffset{0}, canvasWidth{m_width}, canvasHeight{m_height};
        this->GetCanvasParameters(output, pixelVector, thresholds, colOffset, rowOffset, canvasWidth, canvasHeight);

        float **canvas{new float*[canvasHeight]};
        for (int row = 0; row < canvasHeight; ++row)
            canvas[row] = new float[canvasWidth]{};

        // output is a 1 x num_classes x height x width tensor
        // we want the maximum value in the num_classes dimension (1) for every pixel
        auto classes{torch::argmax(output, 1)};
        // the argmax result is a 1 x height x width tensor where each element is a class id
        auto classesAccessor{classes.accessor<long, 3>()};
        const double scaleFactor{std::sqrt(m_height * m_height + m_width * m_width)};
        std::map<int, bool> haveSeenMap;
        for (const auto [ row, col ] : pixelVector)
        {
            const auto cls{classesAccessor[0][row][col]};
            if (cls > 0 && cls < 19)
            {
                const int inner{static_cast<int>(std::round(std::ceil(scaleFactor * thresholds[cls - 1])))};
                const int outer{static_cast<int>(std::round(std::ceil(scaleFactor * thresholds[cls])))};
                this->DrawRing(canvas, row + rowOffset, col + colOffset, inner, outer, 1.f / (outer * outer - inner * inner));
            }
        }

        CartesianPointVector positionVector;
        MakeWirePlaneCoordinatesFromCanvas(*pCaloHitList, canvas, canvasWidth, canvasHeight, colOffset, rowOffset, positionVector);
        if (isU)
            vertexCandidatesU.emplace_back(positionVector.front());
        else if (isV)
            vertexCandidatesV.emplace_back(positionVector.front());
        else
            vertexCandidatesW.emplace_back(positionVector.front());

        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
        for (const auto pos : positionVector)
        {
            Color color{isU ? RED : isV ? GREEN : BLUE};
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "hit", color, 1));

        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

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
        vertexTuples.emplace_back(VertexTuple(this->GetPandora(), vertexCandidatesU.front(), vertexCandidatesV.front(),
            vertexCandidatesW.front()));
    }
    else if (nEmptyLists == 1)
    {
        if (vertexCandidatesU.empty())
        {   // V and W available
            vertexTuples.emplace_back(VertexTuple(this->GetPandora(), vertexCandidatesV.front(), vertexCandidatesW.front(), TPC_VIEW_V,
                TPC_VIEW_W));
        }
        else if (vertexCandidatesV.empty())
        {   // U and W available
            vertexTuples.emplace_back(VertexTuple(this->GetPandora(), vertexCandidatesU.front(), vertexCandidatesW.front(), TPC_VIEW_U,
                TPC_VIEW_W));
        }
        else
        {   // U and V available
            vertexTuples.emplace_back(VertexTuple(this->GetPandora(), vertexCandidatesU.front(), vertexCandidatesV.front(), TPC_VIEW_U,
                TPC_VIEW_W));
        }
    }
    else
    {   // Not enough views to reconstruct a 3D vertex
        std::cout << "Insufficient 2D vertices to reconstruct a 3D vertex" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    if (m_visualise)
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexTuples.front().GetPosition(), "candidate", GREEN, 1));

    if (!vertexTuples.empty())
    {
        const CartesianVector &vertex{vertexTuples.front().GetPosition()};
        CartesianPointVector vertexCandidates;
        vertexCandidates.emplace_back(vertex);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->MakeCandidateVertexList(vertexCandidates));
    }

    if (m_visualise)
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeNetworkInputFromHits(const pandora::CaloHitList &caloHits, LArDLHelper::TorchInput &networkInput,
    PixelVector &pixelVector) const
{
    // Determine the range of coordinates for the view
    float xMin{0.f}, xMax{0.f}, zMin{0.f}, zMax{0.f};
    GetHitRegion(caloHits, xMin, xMax, zMin, zMax);

    // Determine the bin edges - need double precision here for consistency with Python binning
    std::vector<double> xBinEdges(m_width + 1);
    std::vector<double> zBinEdges(m_height + 1);
    xBinEdges[0] = xMin - PY_EPSILON;
    const double dx = ((xMax + PY_EPSILON) - (xMin - PY_EPSILON)) / m_width;
    for (int i = 1; i < m_width + 1; ++i)
        xBinEdges[i] = xBinEdges[i - 1] + dx;
    zBinEdges[0] = zMin - PY_EPSILON;
    const double dz = ((zMax + PY_EPSILON) - (zMin - PY_EPSILON)) / m_height;
    for (int i = 1; i < m_height + 1; ++i)
        zBinEdges[i] = zBinEdges[i - 1] + dz;

    LArDLHelper::InitialiseInput({1, 1, m_height, m_width}, networkInput);
    auto accessor = networkInput.accessor<float, 4>();

    float maxValue{0.f};
    for (const CaloHit *pCaloHit : caloHits)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        const float adc{pCaloHit->GetInputEnergy() < 500 ? pCaloHit->GetInputEnergy() : 500};
        const int pixelX{static_cast<int>(std::floor((x - xBinEdges[0]) / dx))};
        const int pixelZ{(m_height - 1) - static_cast<int>(std::floor((z - zBinEdges[0]) / dz))};
        accessor[0][0][pixelZ][pixelX] += adc;
        if (accessor[0][0][pixelZ][pixelX] > maxValue)
            maxValue = accessor[0][0][pixelZ][pixelX];
    }
    if (maxValue > 0)
    {
        for (int row = 0; row < m_height; ++row)
        {
            for (int col = 0; col < m_width; ++col)
            {
                const float value{accessor[0][0][row][col]};
                accessor[0][0][row][col] = 255.f * value / maxValue;
                if (value > 0)
                    pixelVector.emplace_back(std::make_pair(row, col));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeWirePlaneCoordinatesFromPixels(const CaloHitList &caloHits, const PixelVector &pixelVector,
    CartesianPointVector &positionVector) const
{
    // Determine the range of coordinates for the view
    float xMin{0.f}, xMax{0.f}, zMin{0.f}, zMax{0.f};
    GetHitRegion(caloHits, xMin, xMax, zMin, zMax);

    // Determine the bin size - need double precision here for consistency with Python binning
    const double dx{((xMax + PY_EPSILON) - (xMin - PY_EPSILON)) / m_width};
    xMin -= PY_EPSILON;
    const double dz{((zMax + PY_EPSILON) - (zMin - PY_EPSILON)) / m_height};
    zMin -= PY_EPSILON;
    // Original hit mapping applies a floor operation, so add half a pixel width to get pixel centre
    const float xShift{static_cast<float>(dx * 0.5f)};
    const float zShift{static_cast<float>(dz * 0.5f)};

    for (const auto [ row, col ] : pixelVector)
    {
        const float x{static_cast<float>(col * dx + xMin + xShift)};
        const float z{static_cast<float>(dz * ((m_height - 1) - row) + zMin + zShift)};
        CartesianVector pt(x, 0.f, z);
        positionVector.emplace_back(pt);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeWirePlaneCoordinatesFromCanvas(const CaloHitList &caloHits, float **canvas, const int canvasWidth,
    const int canvasHeight, const int columnOffset, const int rowOffset, CartesianPointVector &positionVector) const
{
    // Determine the range of coordinates for the view
    float xMin{0.f}, xMax{0.f}, zMin{0.f}, zMax{0.f};
    GetHitRegion(caloHits, xMin, xMax, zMin, zMax);

    // Determine the bin size - need double precision here for consistency with Python binning
    const double dx{((xMax + PY_EPSILON) - (xMin - PY_EPSILON)) / m_width};
    xMin -= PY_EPSILON;
    const double dz{((zMax + PY_EPSILON) - (zMin - PY_EPSILON)) / m_height};
    zMin -= PY_EPSILON;
    // Original hit mapping applies a floor operation, so add half a pixel width to get pixel centre
    const float xShift{static_cast<float>(dx * 0.5f)};
    const float zShift{static_cast<float>(dz * 0.5f)};

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
    (void)rowBest;
    (void)colBest;

/*    for (int row = 0; row < canvasHeight; ++row)
        for (int col = 0; col < canvasWidth; ++col)
            if (canvas[row][col] > (best * 0.9f))
            {
                const float x{static_cast<float>((col - columnOffset) * dx + xMin + xShift)};
                const float z{static_cast<float>(dz * ((m_height - 1) - (row - rowOffset)) + zMin + zShift)};
                CartesianVector pt(x, 0.f, z);
                positionVector.emplace_back(pt);
            }*/
    const float x{static_cast<float>((colBest - columnOffset) * dx + xMin + xShift)};
    const float z{static_cast<float>(dz * ((m_height - 1) - (rowBest - rowOffset)) + zMin + zShift)};
    CartesianVector pt(x, 0.f, z);
    positionVector.emplace_back(pt);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::GetCanvasParameters(const LArDLHelper::TorchOutput &networkOutput, const PixelVector &pixelVector,
    const double *const thresholds, int &colOffset, int &rowOffset, int &width, int &height) const
{
    const double scaleFactor{std::sqrt(m_height * m_height + m_width * m_width)};
    // output is a 1 x num_classes x height x width tensor
    // we want the maximum value in the num_classes dimension (1) for every pixel
    auto classes{torch::argmax(networkOutput, 1)};
    // the argmax result is a 1 x height x width tensor where each element is a class id
    auto classesAccessor{classes.accessor<long, 3>()};
    int colOffsetMin{0}, colOffsetMax{0}, rowOffsetMin{0}, rowOffsetMax{0};
    for (const auto [ row, col ] : pixelVector)
    {
        const auto cls{classesAccessor[0][row][col]};
        const double threshold{thresholds[cls]};
        if (threshold > 0. && threshold < 1.)
        {
            const int distance = static_cast<int>(std::round(std::ceil(scaleFactor * threshold)));
            if ((row - distance) < rowOffsetMin)
                rowOffsetMin = row - distance;
            if ((row + distance) > rowOffsetMax)
                rowOffsetMax = row + distance;
            if ((col - distance) < colOffsetMin)
                colOffsetMin = col - distance;
            if ((col + distance) > colOffsetMax)
                colOffsetMax = col + distance;
        }
    }
    colOffset = colOffsetMin < 0 ? -colOffsetMin : 0;
    rowOffset = rowOffsetMin < 0 ? -rowOffsetMin : 0;
    width = std::max(colOffsetMax + colOffset + 1, m_width);
    height = std::max(rowOffsetMax + rowOffset + 1, m_height);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::DrawRing(float **canvas, const int row, const int col, const int inner, const int outer, const float weight) const
{
    // Set the starting position for each circle bounding the ring
    int c1{inner}, r1{0}, c2{outer}, r2{0};
    int inner2{inner * inner}, outer2{outer * outer};
    while (c2 >= r2)
    {
        // Set the output pixel location
        int rp2{r2}, cp2{c2};
        // We're still within the octant for the inner ring, so use the inner pixel location (see Update comment below)
        // Note also that the inner row is always the same as the outer row, so no need to define rp1
        int cp1{c1};
        if (c1 <= r1)
        {   // We've completed the arc of the inner ring already, so just move radially out from here (see Update comment below)
            cp1 = r2;
        }
        // Fill the pixels from inner to outer in the current row and their mirror pixels in the other octants
        for (int c = cp1; c <= cp2; ++c)
        {
            canvas[row + rp2][col + c] += weight;
            if (rp2 != c)
                canvas[row + c][col + rp2] += weight;
            if (rp2 != 0 && cp2 != 0)
            {
                canvas[row - rp2][col - c] += weight;
                if (rp2 != c)
                    canvas[row - c][col - rp2] += weight;
            }
            if (rp2 != 0)
            {
                canvas[row - rp2][col + c] += weight;
                if (rp2 != c)
                    canvas[row + c][col - rp2] += weight;
            }
            if (cp2 != 0)
            {
                canvas[row + rp2][col - c] += weight;
                if (rp2 != c)
                    canvas[row - c][col + rp2] += weight;
            }
        }
        // Only update the inner location while it remains in the octant (outer ring also remains in the octant of course, but the logic of
        // the update means that the inner ring can leave its octant before the outer ring is complete, so we need to stop that)
        if (c1 > r1)
            this->Update(inner2, c1, r1);
        // Update the outer location - increase the row position with every step, decrease the column position if conditions are met
        this->Update(outer2, c2, r2);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::Update(const int radius2, int &col, int &row) const
{
    // Bresenham midpoint circle algorithm to determine if we should update the column position
    // This obscure looking block of code uses bit shifts and integer arithmetic to perform this check as efficiently as possible
    const int a{1 - (col << 2)};
    const int b{col * col + row * row - radius2 + (row << 2) + 1};
    const int c{(a << 2) * b + a * a};
    if (c < 0)
    {
        --col;
        ++row;
    }
    else
        ++row;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList2D, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCParticleList& mcHierarchy)
    const
{
    try
    {
        for (const auto [ mc, hits ] : mcToHitsMap)
        {   (void) hits;
            mcHierarchy.push_back(mc);
            LArMCParticleHelper::GetAllAncestorMCParticles(mc, mcHierarchy);
        }
    }
    catch (const StatusCodeException &e)
    {
        return e.GetStatusCode();
    }

    // Move the neutrino to the front of the list
    auto pivot = std::find_if(mcHierarchy.begin(), mcHierarchy.end(),
        [](const MCParticle* mc) -> bool { return LArMCParticleHelper::IsNeutrino(mc); });
    (void) pivot;
    if (pivot != mcHierarchy.end())
        std::rotate(mcHierarchy.begin(), pivot, std::next(pivot));
    else
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const
{
    xMin = std::numeric_limits<float>::max(); xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max(); zMax = -std::numeric_limits<float>::max();
    // Find the range of x and z values in the view
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        if (x < xMin)
            xMin = x;
        if (x > xMax)
            xMax = x;
        if (z < zMin)
            zMin = z;
        if (z > zMax)
            zMax = z;
    }
    // Avoid unreasonable rescaling of very small hit regions
    const float xRange{xMax - xMin}, zRange{zMax - zMin};
    if (2 * xRange < m_width)
    {
        const float padding{0.5f * (0.5f * m_width - xRange)};
        xMin -= padding;
        xMax += padding;
    }
    if (2 * zRange < m_height)
    {
        const float padding{0.5f * (0.5f * m_height - zRange)};
        zMin -= padding;
        zMax += padding;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeCandidateVertexList(const CartesianPointVector &positions)
{
    const VertexList *pVertexList{nullptr};
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (const CartesianVector position : positions)
    {
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pVertex(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode",
        m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise",
        m_visualise));

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameU", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelU);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameV", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelV);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameW", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelW);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PixelShift", m_pixelShift));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PixelScale", m_pixelScale));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingAlgorithm::VertexTuple::VertexTuple(const Pandora &pandora, const CartesianVector &vertexU, const CartesianVector &vertexV,
    const CartesianVector &vertexW) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    LArGeometryHelper::MergeThreePositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, vertexU, vertexV, vertexW, m_pos, m_chi2);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingAlgorithm::VertexTuple::VertexTuple(const Pandora &pandora, const CartesianVector &vertex1, const CartesianVector &vertex2,
    const HitType view1, const HitType view2) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    LArGeometryHelper::MergeTwoPositions3D(pandora, view1, view2, vertex1, vertex2, m_pos, m_chi2);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &DlVertexingAlgorithm::VertexTuple::GetPosition() const
{
    return m_pos;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float DlVertexingAlgorithm::VertexTuple::GetChi2() const
{
    return m_chi2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::string DlVertexingAlgorithm::VertexTuple::ToString() const
{
    const float x{m_pos.GetX()}, y{m_pos.GetY()}, z{m_pos.GetZ()};
    return "3D pos: (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")   X2 = " + std::to_string(m_chi2);
}

} // namespace lar_dl_content

