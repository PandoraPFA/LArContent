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
                if (LArMCParticleHelper::GetHierarchyTier(mc) == 1)
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
                }
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
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    CaloHitList vertexCandidatesU, vertexCandidatesV, vertexCandidatesW;
    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        HitType view{pCaloHitList->front()->GetHitType()};
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

        if (isW && m_visualise)
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitList, listName, BLACK));

        const int imageWidth = 256;
        const int imageHeight = 256;
        // Get bounds of hit region
        float xMin{}; float xMax{}; float zMin{}; float zMax{};
        GetHitRegion(*pCaloHitList, xMin, xMax, zMin, zMax);
        const float xRange = (xMax - xMin) < std::numeric_limits<float>::epsilon() ? 1.f : xMax - xMin;
        const float zRange = (zMax - zMin) < std::numeric_limits<float>::epsilon() ? 1.f : zMax - zMin;

        // ATTN - Pixel normalisation must match that used during training
        const float baseline{255.f};
        const float pixelValue = (baseline - m_pixelShift) / m_pixelScale;

        PixelToCaloHitsMap pixelToCaloHits;
        LArDLHelper::TorchInput input;
        LArDLHelper::InitialiseInput({1, 1, 256, 256}, input);
        auto accessor = input.accessor<float, 4>();
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const float x{pCaloHit->GetPositionVector().GetX()};
            const float z{pCaloHit->GetPositionVector().GetZ()};
            // Determine which pixel the hit will be assigned to, being careful to avoid an out of bounds index
            const float relativeX{(x - xMin) / xRange};
            const float relativeZ{1.f - (z - zMin) / zRange};
            const int pixelX{static_cast<int>(relativeX < 1.f ? std::floor(imageWidth * relativeX) : std::floor(imageWidth * relativeX - 0.5f))};
            const int pixelZ{static_cast<int>(relativeZ < 1.f ? std::floor(imageHeight * relativeZ) : std::floor(imageHeight * relativeZ - 0.5f))};
            accessor[0][0][pixelZ][pixelX] = pixelValue;
            // We want the get or create if not found behaviour
            pixelToCaloHits[std::make_pair(pixelZ, pixelX)].push_back(pCaloHit);
        }

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
        auto outputAccessor = output.accessor<float, 4>();

        // Find the maximum vertex score
        float maxScore{0.f};
        for (int pz = 0; pz < imageHeight; ++pz)
        {
            for (int px = 0; px < imageWidth; ++px)
            {
                const float actPrimary = outputAccessor[0][1][pz][px];
                if (actPrimary > maxScore)
                    maxScore = actPrimary;
            }
        }
        const float threshold{maxScore - std::numeric_limits<float>::epsilon()};

        // Tag the hits associated with each active pixel as candidate vertices
        CaloHitList vertexHits;
        for (int pz = 0; pz < imageHeight; ++pz)
        {
            for (int px = 0; px < imageWidth; ++px)
            {
                const float actPrimary = outputAccessor[0][1][pz][px];
                //const float actSecondary = exp(outputAccessor[0][2][pz][px]);
                const float actNull = outputAccessor[0][0][pz][px];
                if (actPrimary > actNull && actPrimary > threshold)
                {
                    // We want the get or create if not found behaviour because vertex pixels are not guaranteed to contain hits.
                    // Nothing will be displayed in this case, but the loop will continue to run without issue.
                    for (const CaloHit *pCaloHit : pixelToCaloHits[std::make_pair(pz, px)])
                    {
                        vertexHits.push_back(pCaloHit);
                        if (isW)
                            vertexCandidatesW.emplace_back(pCaloHit);
                        else if (isV)
                            vertexCandidatesV.emplace_back(pCaloHit);
                        else
                            vertexCandidatesU.emplace_back(pCaloHit);
                    }
                }
            }
        }
    }

    int nEmptyLists{0};
    if (vertexCandidatesU.empty())
        ++nEmptyLists;
    if (vertexCandidatesV.empty())
        ++nEmptyLists;
    if (vertexCandidatesW.empty())
        ++nEmptyLists;
    std::vector<CaloHitTuple> hits3D;
    if (nEmptyLists == 0)
    {
        for (const CaloHit *pCaloHitU : vertexCandidatesU)
        {
            const CartesianVector &posU(pCaloHitU->GetPositionVector());
            for (const CaloHit *pCaloHitV : vertexCandidatesV)
            {
                const CartesianVector &posV(pCaloHitV->GetPositionVector());
                for (const CaloHit *pCaloHitW : vertexCandidatesW)
                {
                    const CartesianVector &posW(pCaloHitW->GetPositionVector());
                    const float dxUV{std::fabs(posU.GetX() - posV.GetX())}, dxUW{std::fabs(posU.GetX() - posW.GetX())}, dxVW{std::fabs(posV.GetX() - posW.GetX())};
                    if (dxUV < 1.5 && dxUW < 1.5 && dxVW < 1.5)
                    {
                        CaloHitTuple hit3D(this->GetPandora(), pCaloHitU, pCaloHitV, pCaloHitW);
                        hits3D.push_back(hit3D);
                    }
                }
            }
        }
    }
    else if (nEmptyLists == 1)
    {
        if (vertexCandidatesU.empty())
        {   // V and W available
            for (const CaloHit *pCaloHitV : vertexCandidatesV)
            {
                const CartesianVector &posV(pCaloHitV->GetPositionVector());
                for (const CaloHit *pCaloHitW : vertexCandidatesW)
                {
                    const CartesianVector &posW(pCaloHitW->GetPositionVector());
                    const float dxVW{std::fabs(posV.GetX() - posW.GetX())};
                    if (dxVW < 1.5)
                    {
                        CaloHitTuple hit3D(this->GetPandora(), nullptr, pCaloHitV, pCaloHitW);
                        hits3D.push_back(hit3D);
                    }
                }
            }
        }
        else if (vertexCandidatesV.empty())
        {   // U and W available
            for (const CaloHit *pCaloHitU : vertexCandidatesU)
            {
                const CartesianVector &posU(pCaloHitU->GetPositionVector());
                for (const CaloHit *pCaloHitW : vertexCandidatesW)
                {
                    const CartesianVector &posW(pCaloHitW->GetPositionVector());
                    const float dxUW{std::fabs(posU.GetX() - posW.GetX())};
                    if (dxUW < 1.5)
                    {
                        CaloHitTuple hit3D(this->GetPandora(), pCaloHitU, nullptr, pCaloHitW);
                        hits3D.push_back(hit3D);
                    }
                }
            }
        }
        else
        {   // U and V available
            for (const CaloHit *pCaloHitU : vertexCandidatesU)
            {
                const CartesianVector &posU(pCaloHitU->GetPositionVector());
                for (const CaloHit *pCaloHitV : vertexCandidatesV)
                {
                    const CartesianVector &posV(pCaloHitV->GetPositionVector());
                    const float dxUV{std::fabs(posU.GetX() - posV.GetX())};
                    if (dxUV < 1.5)
                    {
                        CaloHitTuple hit3D(this->GetPandora(), pCaloHitU, pCaloHitV, nullptr);
                        hits3D.push_back(hit3D);
                    }
                }
            }
        }
    }
    else
    {   // Not enough hits in different views to return any candidates
        return STATUS_CODE_SUCCESS;
    }

    std::sort(hits3D.begin(), hits3D.end(), [] (const CaloHitTuple &t1, const CaloHitTuple &t2) -> bool { return t1.GetChi2() < t2.GetChi2(); });
    // Consider a random sample of a smaller number of hits to avoid combinatorial explosion
    std::map<const CaloHit*, bool> uHitMap, vHitMap, wHitMap;
    std::vector<CaloHitTuple> candidates3D;
    FloatVector drs;
    CartesianPointVector positions;
    for (const CaloHitTuple &hit3D : hits3D)
    {
        const CaloHit *pCaloHitU{hit3D.GetCaloHitU()}, *pCaloHitV{hit3D.GetCaloHitV()}, *pCaloHitW{hit3D.GetCaloHitW()};
        if (pCaloHitU && uHitMap.find(pCaloHitU) != uHitMap.end())
            continue;
        if (pCaloHitV && uHitMap.find(pCaloHitV) != uHitMap.end())
            continue;
        if (pCaloHitW && uHitMap.find(pCaloHitW) != uHitMap.end())
            continue;
        candidates3D.push_back(hit3D);
        if (pCaloHitU)
            uHitMap[pCaloHitU] = true;
        if (pCaloHitV)
            vHitMap[pCaloHitV] = true;
        if (pCaloHitW)
            wHitMap[pCaloHitW] = true;
        const CartesianVector &position{hit3D.GetPosition()};
        positions.emplace_back(position);

        if (m_visualise)
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &hit3D.GetPosition(), "candidate", GREEN, 1));
    }

    if (!positions.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->MakeCandidateVertexList(positions));

    if (m_visualise)
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
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

void DlVertexingAlgorithm::GetHitRegion(const CaloHitList& caloHitList, float& xMin, float& xMax, float& zMin, float& zMax) const
{
    xMin = std::numeric_limits<float>::max(); xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max(); zMax = -std::numeric_limits<float>::max();
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
        LArDLHelper::LoadModel(modelName, m_modelU);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameV", modelName));
        LArDLHelper::LoadModel(modelName, m_modelV);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameW", modelName));
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

DlVertexingAlgorithm::CaloHitTuple::CaloHitTuple(const pandora::Pandora &pandora, const CaloHit *pCaloHitU, const CaloHit *pCaloHitV,
    const CaloHit *pCaloHitW) :
    m_pCaloHitU{pCaloHitU},
    m_pCaloHitV{pCaloHitV},
    m_pCaloHitW{pCaloHitW},
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    const bool hasU{m_pCaloHitU}, hasV{m_pCaloHitV}, hasW{m_pCaloHitW};
    if (hasW && hasV && hasU)
        LArGeometryHelper::MergeThreePositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W,
            m_pCaloHitU->GetPositionVector(), m_pCaloHitV->GetPositionVector(), m_pCaloHitW->GetPositionVector(), m_pos, m_chi2);
    else if (hasW && hasV)
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_V, TPC_VIEW_W,
            m_pCaloHitV->GetPositionVector(), m_pCaloHitW->GetPositionVector(), m_pos, m_chi2);
    else if (hasW && hasU)
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_U, TPC_VIEW_W,
            m_pCaloHitU->GetPositionVector(), m_pCaloHitW->GetPositionVector(), m_pos, m_chi2);
    else if (hasV && hasU)
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V,
            m_pCaloHitU->GetPositionVector(), m_pCaloHitV->GetPositionVector(), m_pos, m_chi2);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &DlVertexingAlgorithm::CaloHitTuple::GetPosition() const
{
    return m_pos;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float DlVertexingAlgorithm::CaloHitTuple::GetChi2() const
{
    return m_chi2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::string DlVertexingAlgorithm::CaloHitTuple::ToString() const
{
    const CartesianVector &uPos{m_pCaloHitU ? m_pCaloHitU->GetPositionVector() : CartesianVector(0, 0, 0)};
    const CartesianVector &vPos{m_pCaloHitV ? m_pCaloHitV->GetPositionVector() : CartesianVector(0, 0, 0)};
    const CartesianVector &wPos{m_pCaloHitW ? m_pCaloHitW->GetPositionVector() : CartesianVector(0, 0, 0)};
    const float ux{uPos.GetX()}, uz{uPos.GetZ()};
    const float vx{vPos.GetX()}, vz{vPos.GetZ()};
    const float wx{wPos.GetX()}, wz{wPos.GetZ()};
    const float x{m_pos.GetX()}, y{m_pos.GetY()}, z{m_pos.GetZ()};

    return "U: (" + std::to_string(ux) + ", " + std::to_string(uz) +  ")   V: (" + std::to_string(vx) + ", " + std::to_string(vz) +
        ")   W: (" + std::to_string(wx) + ", " + std::to_string(wz) + ")   3D: (" +
        std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")   X2 = " + std::to_string(m_chi2);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *DlVertexingAlgorithm::CaloHitTuple::GetCaloHitU() const
{
    return m_pCaloHitU;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *DlVertexingAlgorithm::CaloHitTuple::GetCaloHitV() const
{
    return m_pCaloHitV;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *DlVertexingAlgorithm::CaloHitTuple::GetCaloHitW() const
{
    return m_pCaloHitW;
}

} // namespace lar_dl_content

