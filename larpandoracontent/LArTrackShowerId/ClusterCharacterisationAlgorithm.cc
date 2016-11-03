/**
 *  @file   larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm() :
    m_overwriteExistingId(false),
    m_useUnavailableClusters(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ClusterCharacterisationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!m_overwriteExistingId && (UNKNOWN_PARTICLE_TYPE != pCluster->GetParticleId()))
                continue;

            if (!m_useUnavailableClusters && !PandoraContentApi::IsAvailable(*this, pCluster))
                continue;

            PandoraContentApi::Cluster::Metadata metadata;

            if (this->IsClearTrack(pCluster))
            {
                metadata.m_particleId = MU_MINUS;
            }
            else
            {
                metadata.m_particleId = E_MINUS;
            }

            if (pCluster->GetParticleId() != metadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < 6)
        return false;

    // Quantities related to sliding linear fit
    float straightLineLength(-1.f), integratedPathLength(-1.f);
    float rTMin(+std::numeric_limits<float>::max()), rTMax(-std::numeric_limits<float>::max());

    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, 5, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector globalMaxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());
        straightLineLength = (globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude();

        integratedPathLength = 0.f;
        CartesianVector previousFitPosition(globalMinLayerPosition);
        const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

        for (const auto &mapEntry : layerFitResultMap)
        {
            rTMin = std::min(rTMin, static_cast<float>(mapEntry.second.GetFitT()));
            rTMax = std::max(rTMax, static_cast<float>(mapEntry.second.GetFitT()));

            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
            integratedPathLength += (thisFitPosition - previousFitPosition).GetMagnitude();
            previousFitPosition = thisFitPosition;
        }
    }
    catch (const StatusCodeException &)
    {
    }

    if (straightLineLength < std::numeric_limits<float>::epsilon())
        return false;

    if (straightLineLength > 80.f)
        return true;

    if (integratedPathLength / straightLineLength > 1.005f)
        return false;

    const float rTWidth(rTMax - rTMin);

    if (rTWidth / straightLineLength > 0.05f)
        return false;

    // Distance to interaction vertex
    const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*this, pVertexList);

    if (pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == pVertexList->front()->GetVertexType()))
    {
        const Vertex *const pVertex(pVertexList->front());
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
        const float vertexDistance(LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster));

        if (vertexDistance / straightLineLength > 0.5f)
            return false;
    }

    // Shower fit width
    float showerFitWidth(-1.f);

    try
    {
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        if (layerFitResultMapS.size() > 1)
        {
            showerFitWidth = 0.f;

            for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
            {
                LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
                LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

                if ((layerFitResultMapP.end() == iterP) || (layerFitResultMapN.end() == iterN))
                    continue;

                showerFitWidth += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
            }
        }
    }
    catch (const StatusCodeException &)
    {
    }

    if (showerFitWidth < std::numeric_limits<float>::epsilon())
        return false;

    if (showerFitWidth / straightLineLength > 0.35f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverwriteExistingId", m_overwriteExistingId));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnavailableClusters", m_useUnavailableClusters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
