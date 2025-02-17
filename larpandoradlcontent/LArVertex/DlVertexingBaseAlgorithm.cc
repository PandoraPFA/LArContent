/**
 *  @file   larpandoradlcontent/LArVertex/DlVertexingBaseAlgorithm.cc
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

#include "larpandoradlcontent/LArVertex/DlVertexingBaseAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlVertexingBaseAlgorithm::DlVertexingBaseAlgorithm() :
    m_trainingMode{false},
    m_trainingOutputFile{""},
    m_pass{1},
    m_nClasses{0},
    m_height{256},
    m_width{256},
    m_driftStep{0.5f},
    m_volumeType{"dune_fd_hd"}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingBaseAlgorithm::~DlVertexingBaseAlgorithm()
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingBaseAlgorithm::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();
    // Find the range of x and z values in the view
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        xMin = std::min(x, xMin);
        xMax = std::max(x, xMax);
        zMin = std::min(z, zMin);
        zMax = std::max(z, zMax);
    }

    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const HitType view{caloHitList.front()->GetHitType()};
    const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
    if (!(isU || isV || isW))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());

    if (m_pass > 1)
    {
        const VertexList *pVertexList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pVertexList));
        if (pVertexList->empty())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const CartesianVector &vertex{pVertexList->front()->GetPosition()};

        // Get hit distribution left/right asymmetry
        int nHitsLeft{0}, nHitsRight{0};
        const double xVtx{vertex.GetX()};
        for (const std::string &listname : m_caloHitListNames)
        {
            const CaloHitList *pCaloHitList(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
            if (pCaloHitList->empty())
                continue;
            for (const CaloHit *const pCaloHit : *pCaloHitList)
            {
                const CartesianVector &pos{pCaloHit->GetPositionVector()};
                if (pos.GetX() <= xVtx)
                    ++nHitsLeft;
                else
                    ++nHitsRight;
            }
        }
        const int nHitsTotal{nHitsLeft + nHitsRight};
        if (nHitsTotal == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const float xAsymmetry{nHitsLeft / static_cast<float>(nHitsTotal)};

        // Vertices
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        double zVtx{0.};
        if (isW)
            zVtx += transform->YZtoW(vertex.GetY(), vertex.GetZ());
        else if (isV)
            zVtx += transform->YZtoV(vertex.GetY(), vertex.GetZ());
        else
            zVtx = transform->YZtoU(vertex.GetY(), vertex.GetZ());

        // Get hit distribution upstream/downstream asymmetry
        int nHitsUpstream{0}, nHitsDownstream{0};
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            if (pos.GetZ() <= zVtx)
                ++nHitsUpstream;
            else
                ++nHitsDownstream;
        }
        const int nHitsViewTotal{nHitsUpstream + nHitsDownstream};
        if (nHitsViewTotal == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const float zAsymmetry{nHitsUpstream / static_cast<float>(nHitsViewTotal)};

        const float xSpan{m_driftStep * (m_width - 1)};
        xMin = xVtx - xAsymmetry * xSpan;
        xMax = xMin + (m_driftStep * (m_width - 1));
        const float zSpan{pitch * (m_height - 1)};
        zMin = zVtx - zAsymmetry * zSpan;
        zMax = zMin + zSpan;
    }

    // Avoid unreasonable rescaling of very small hit regions, pixels are assumed to be 0.5cm in x and wire pitch in z
    // ATTN: Rescaling is to a size 1 pixel smaller than the intended image to ensure all hits fit within an imaged binned
    // to be one pixel wider than this
    const float xRange{xMax - xMin}, zRange{zMax - zMin};
    const float minXSpan{m_driftStep * (m_width - 1)};
    if (xRange < minXSpan)
    {
        const float padding{0.5f * (minXSpan - xRange)};
        xMin -= padding;
        xMax += padding;
    }
    const float minZSpan{pitch * (m_height - 1)};
    if (zRange < minZSpan)
    {
        const float padding{0.5f * (minZSpan - zRange)};
        zMin -= padding;
        zMax += padding;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingBaseAlgorithm::GetCanvasParameters(const LArDLHelper::TorchOutput &networkOutput, const PixelVector &pixelVector,
    int &colOffset, int &rowOffset, int &width, int &height) const
{
    const double scaleFactor{std::sqrt(m_height * m_height + m_width * m_width)};
    // output is a 1 x num_classes x height x width tensor
    // we want the maximum value in the num_classes dimension (1) for every pixel
    auto classes{torch::argmax(networkOutput, 1)};
    // the argmax result is a 1 x height x width tensor where each element is a class id
    auto classesAccessor{classes.accessor<int64_t, 3>()};
    int colOffsetMin{0}, colOffsetMax{0}, rowOffsetMin{0}, rowOffsetMax{0};
    for (const auto &[row, col] : pixelVector)
    {
        const auto cls{classesAccessor[0][row][col]};
        const double threshold{m_thresholds[cls]};
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

StatusCode DlVertexingBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Pass", m_pass));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_height));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_width));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DriftStep", m_driftStep));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DistanceThresholds", m_thresholds));
    m_nClasses = m_thresholds.size() - 1;
    if (m_pass > 1)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));
    }

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
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VolumeType", m_volumeType));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
