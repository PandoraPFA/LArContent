/**
 *  @file   larpandoradlcontent/LArVertex/DlVertexingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning vertexing algorithm.
 *
 *  $Log: $
 */

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoradlcontent/LArVertex/DlVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlVertexingAlgorithm::DlVertexingAlgorithm():
    m_trainingMode{false},
    m_trainingOutputFile{""},
    m_pixelShift{0.f},
    m_pixelScale{1.f}
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
            {
                const double w{transform->YZtoW(vertex.GetY(), vertex.GetZ())};
                featureVector.push_back(w);
            }
            else if (isV)
            {
                const double v{transform->YZtoV(vertex.GetY(), vertex.GetZ())};
                featureVector.push_back(v);
            }
            else
            {
                const double u{transform->YZtoU(vertex.GetY(), vertex.GetZ())};
                featureVector.push_back(u);
            }
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

                const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()};
                featureVector.push_back(static_cast<double>(x));
                featureVector.push_back(static_cast<double>(z));
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

    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);
        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;
        // Temporary
        if (!isW)
            continue;

        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitList, listName, BLACK));

        const int imageWidth = 256;
        const int imageHeight = 256;
        // Get bounds of hit region
        float xMin{}; float xMax{}; float zMin{}; float zMax{};
        GetHitRegion(*pCaloHitList, xMin, xMax, zMin, zMax);
        const float xRange = (xMax - xMin) < std::numeric_limits<float>::epsilon() ? 1.f : xMax - xMin;
        const float zRange = (zMax - zMin) < std::numeric_limits<float>::epsilon() ? 1.f : zMax - zMin;

        // Pixel normalisation must match that used during training
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

        // Tag the hits associated with each active pixel as candidate vertices
        CaloHitList vertexHits;
        for (int pz = 0; pz < imageHeight; ++pz)
        {
            for (int px = 0; px < imageWidth; ++px)
            {
                const float actPrimary = outputAccessor[0][1][pz][px];
                //const float actSecondary = exp(outputAccessor[0][2][pz][px]);
                const float actNull = outputAccessor[0][0][pz][px];
                if (actPrimary > actNull)
                {
                    // We want the get or create if not found behaviour because vertex pixels are not guaranteed to contain hits.
                    // Nothing will be displayed in this case, but the loop will continue to run without issue.
                    for (const CaloHit *pCaloHit : pixelToCaloHits[std::make_pair(pz, px)])
                        vertexHits.push_back(pCaloHit);
                }
            }
        }
        std::string vertexName{"vtx_" + listName};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vertexHits, vertexName, RED));
    }
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
    parameters.m_foldBackHierarchy = false;
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
        {
            (void) hits;
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

StatusCode DlVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode",
        m_trainingMode));

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
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content

