/**
 *  @file   larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.cc
 *
 *  @brief  Implementation of the track shower id feature fool class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerFitFeatureTool::ShowerFitFeatureTool() :
    m_slidingShowerFitWindow(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerFitFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    featureVector.push_back(CutClusterCharacterisationAlgorithm::GetShowerFitWidth(pAlgorithm, pCluster, m_slidingShowerFitWindow));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingShowerFitWindow", m_slidingShowerFitWindow));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void NHitsFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster * const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    featureVector.push_back(pCluster->GetNCaloHits());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NHitsFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

    LinearFitFeatureTool::LinearFitFeatureTool() :
    m_slidingLinearFitWindow(3),
    m_addDiffWithStraightLine(true),
    m_addDTDLWidth(true),
    m_addMaxFitGapLength(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LinearFitFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
const pandora::Cluster * const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float dTdLWidth(-1.f), straightLineLengthLarge(-1.f), diffWithStraightLine(0.f), maxFitGapLength(-1.f);
    this->CalculateVariablesSlidingLinearFit(pCluster, straightLineLengthLarge, diffWithStraightLine, dTdLWidth, maxFitGapLength);

    featureVector.push_back(straightLineLengthLarge);

    if (m_addDiffWithStraightLine)
        featureVector.push_back(diffWithStraightLine);

    if (m_addDTDLWidth)
        featureVector.push_back(dTdLWidth);

    if (m_addMaxFitGapLength)
        featureVector.push_back(maxFitGapLength);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster,
float &straightLineLengthLarge, float &diffWithStraightLine, float &dTdLWidth, float &maxFitGapLength) const
{
    try
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        straightLineLengthLarge = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();

        CartesianVector previousFitPosition(slidingFitResult.GetGlobalMinLayerPosition());
        float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());

        for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
        {
            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);

            LayerFitResultMap::const_iterator iterLarge = slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(mapEntry.second.GetL()));

            if (slidingFitResultLarge.GetLayerFitResultMap().end() == iterLarge)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            diffWithStraightLine += std::fabs(mapEntry.second.GetFitT() - iterLarge->second.GetFitT());

            const float thisGapLength((thisFitPosition - previousFitPosition).GetMagnitude());
            const float minZ(std::min(thisFitPosition.GetZ(), previousFitPosition.GetZ()));
            const float maxZ(std::max(thisFitPosition.GetZ(), previousFitPosition.GetZ()));

            if ((maxZ - minZ) > std::numeric_limits<float>::epsilon())
            {
                const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, hitType));
                const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));

                if (correctedGapLength > maxFitGapLength)
                    maxFitGapLength = correctedGapLength;
            }

            dTdLMin = std::min(dTdLMin, static_cast<float>(mapEntry.second.GetGradient()));
            dTdLMax = std::max(dTdLMax, static_cast<float>(mapEntry.second.GetGradient()));
            previousFitPosition = thisFitPosition;
        }

        dTdLWidth = dTdLMax - dTdLMin;
    }
    catch (const StatusCodeException &)
    {
        dTdLWidth = -1.f;
        maxFitGapLength = -1.f;
        diffWithStraightLine = -1.f;
        straightLineLengthLarge = -1.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LinearFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddDiffWithStraightLine", m_addDiffWithStraightLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddDTDLWidth", m_addDTDLWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddMaxFitGapLength", m_addMaxFitGapLength));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

 NNearbyClustersFeatureTool::NNearbyClustersFeatureTool() :
    m_minClusterCaloHits(6),
    m_nearbyClusterDistance(2.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NNearbyClustersFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterList candidateClusters;

    for (const std::string &clusterListName : m_clusterListNames)
    {
        const pandora::ClusterList *pClusterList = nullptr;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, clusterListName, pClusterList));

        for (const Cluster *const pCandidateCluster : *pClusterList)
        {
            if ((LArClusterHelper::GetClusterHitType(pCluster)) != (LArClusterHelper::GetClusterHitType(pCandidateCluster)))
                continue;

            if ((pCandidateCluster != pCluster) && (pCandidateCluster->GetNCaloHits() >= m_minClusterCaloHits))
                candidateClusters.push_back(pCandidateCluster);
        }
    }

    int nNearbyClusters(0);

    for (const Cluster *const pCandidateCluster : candidateClusters)
    {
        if ((pCluster != pCandidateCluster) && (LArClusterHelper::GetClosestDistance(pCluster, pCandidateCluster) <= m_nearbyClusterDistance))
            ++nNearbyClusters;
    }

    featureVector.push_back(nNearbyClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NNearbyClustersFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NearbyClusterDistance", m_nearbyClusterDistance));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MipEnergyFeatureTool::MipEnergyFeatureTool() :
    m_mipCorrectionPerHit(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MipEnergyFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float mipEnergy(0.f);

    for (const OrderedCaloHitList::value_type &layerIter : pCluster->GetOrderedCaloHitList())
    {
        for (const CaloHit *const pCaloHit : *layerIter.second)
            mipEnergy += pCaloHit->GetMipEquivalentEnergy();
    }

    featureVector.push_back(mipEnergy * m_mipCorrectionPerHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MipEnergyFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MipCorrectionPerHit", m_mipCorrectionPerHit));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void VertexDistanceFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    featureVector.push_back(CutClusterCharacterisationAlgorithm::GetVertexDistance(pAlgorithm, pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexDistanceFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
