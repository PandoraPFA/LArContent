/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm() :
    m_slidingFitWindow(10),
    m_minHitsInCluster(20),
    m_maxLayerGapFraction(0.2f),
    m_maxWidthPerUnitLength(0.15f),
    m_maxShowerLength(1000.f),
    m_useDetectorGaps(true)
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
            if (this->IsClearTrack(pCluster))
            {
                PandoraContentApi::Cluster::Metadata metadata;
                metadata.m_particleId = MU_MINUS;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AlterMetadata(*this, pCluster, metadata));
            }
            else
            {
                PandoraContentApi::Cluster::Metadata metadata;
                metadata.m_particleId = E_MINUS;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AlterMetadata(*this, pCluster, metadata));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minHitsInCluster)
        return false;

    try
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        if (layerFitResultMapS.size() < 2)
            return false;

        CartesianVector globalMinLayerPositionOnAxis(0.f, 0.f, 0.f), globalMaxLayerPositionOnAxis(0.f, 0.f, 0.f);
        showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, globalMinLayerPositionOnAxis);
        showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, globalMaxLayerPositionOnAxis);
        const float straightLinePathLength((globalMaxLayerPositionOnAxis - globalMinLayerPositionOnAxis).GetMagnitude());

        if (straightLinePathLength < std::numeric_limits<float>::epsilon())
            return false;

        // TODO: Improve on a straight selection cut above a certain distance
        if (straightLinePathLength > m_maxShowerLength)
            return true;

        float widthSum(0.f);
        CartesianVector previousLayerPosition(globalMinLayerPositionOnAxis);
        float biggestGapLength(std::numeric_limits<float>::epsilon());

        for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
        {
            CartesianVector thisLayerPosition(0.f, 0.f, 0.f);
            showerFitResult.GetShowerFitResult().GetGlobalPosition(iterS->second.GetL(), 0.f, thisLayerPosition);
            const float thisGapLength((thisLayerPosition - previousLayerPosition).GetMagnitude());

            if (thisGapLength > biggestGapLength)
            {
                const float minZ(std::min(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));
                const float maxZ(std::max(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));

                if ((maxZ - minZ) < std::numeric_limits<float>::epsilon())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const float gapZ(m_useDetectorGaps ? LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, LArClusterHelper::GetClusterHitType(pCluster)) : 0.f);
                const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));

                if (correctedGapLength > biggestGapLength)
                    biggestGapLength = correctedGapLength;
            }

            previousLayerPosition = thisLayerPosition;

            LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
            LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

            if ((layerFitResultMapP.end() == iterP) || (layerFitResultMapN.end() == iterN))
                continue;

            widthSum += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
        }

        const float layerGapFraction(biggestGapLength / straightLinePathLength);

        if (layerGapFraction > m_maxLayerGapFraction)
            return false;

        // ATTN: Could also try to recover long (sideways) tracks, with cuts on nHits and nHitsPerUnitLength
        const float widthPerUnitLength(widthSum / straightLinePathLength);

        if (widthPerUnitLength < m_maxWidthPerUnitLength)
            return true;
    }
    catch (StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsInCluster", m_minHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLayerGapFraction", m_maxLayerGapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxWidthPerUnitLength", m_maxWidthPerUnitLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerLength", m_maxShowerLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDetectorGaps", m_useDetectorGaps));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
