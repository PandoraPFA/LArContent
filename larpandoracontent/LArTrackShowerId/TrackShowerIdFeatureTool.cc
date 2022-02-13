/**
 *  @file   larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.cc
 *
 *  @brief  Implementation of the track shower id feature fool class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h"

using namespace pandora;

namespace lar_content
{

TwoDShowerFitFeatureTool::TwoDShowerFitFeatureTool() : m_slidingShowerFitWindow(3), m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDShowerFitFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float ratio(-1.f);
    try
    {
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const float straightLineLength =
            (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        if (straightLineLength > std::numeric_limits<float>::epsilon())
            ratio = (CutClusterCharacterisationAlgorithm::GetShowerFitWidth(pAlgorithm, pCluster, m_slidingShowerFitWindow)) / straightLineLength;
    }
    catch (const StatusCodeException &)
    {
        ratio = -1.f;
    }
    featureVector.push_back(ratio);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDShowerFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingShowerFitWindow", m_slidingShowerFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoDLinearFitFeatureTool::TwoDLinearFitFeatureTool() : m_slidingLinearFitWindow(3), m_slidingLinearFitWindowLarge(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDLinearFitFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float dTdLWidth(-1.f), straightLineLengthLarge(-1.f), diffWithStraightLineMean(-1.f), diffWithStraightLineSigma(-1.f),
        maxFitGapLength(-1.f), rmsSlidingLinearFit(-1.f);
    this->CalculateVariablesSlidingLinearFit(pCluster, straightLineLengthLarge, diffWithStraightLineMean, diffWithStraightLineSigma,
        dTdLWidth, maxFitGapLength, rmsSlidingLinearFit);

    if (straightLineLengthLarge > std::numeric_limits<float>::epsilon())
    {
        diffWithStraightLineMean /= straightLineLengthLarge;
        diffWithStraightLineSigma /= straightLineLengthLarge;
        dTdLWidth /= straightLineLengthLarge;
        maxFitGapLength /= straightLineLengthLarge;
        rmsSlidingLinearFit /= straightLineLengthLarge;
    }

    featureVector.push_back(straightLineLengthLarge);
    featureVector.push_back(diffWithStraightLineMean);
    featureVector.push_back(diffWithStraightLineSigma);
    featureVector.push_back(dTdLWidth);
    featureVector.push_back(maxFitGapLength);
    featureVector.push_back(rmsSlidingLinearFit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDLinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge,
    float &diffWithStraightLineMean, float &diffWithStraightLineSigma, float &dTdLWidth, float &maxFitGapLength, float &rmsSlidingLinearFit) const
{
    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResultLarge(
            pCluster, m_slidingLinearFitWindowLarge, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        if (slidingFitResult.GetLayerFitResultMap().empty())
            throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

        straightLineLengthLarge =
            (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        rmsSlidingLinearFit = 0.f;

        FloatVector diffWithStraightLineVector;
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        CartesianVector previousFitPosition(slidingFitResult.GetGlobalMinLayerPosition());
        float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());

        for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
        {
            const LayerFitResult &layerFitResult(mapEntry.second);
            rmsSlidingLinearFit += layerFitResult.GetRms();

            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(layerFitResult.GetL(), layerFitResult.GetFitT(), thisFitPosition);

            LayerFitResultMap::const_iterator iterLarge =
                slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(layerFitResult.GetL()));

            if (slidingFitResultLarge.GetLayerFitResultMap().end() == iterLarge)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            diffWithStraightLineVector.push_back(static_cast<float>(std::fabs(layerFitResult.GetFitT() - iterLarge->second.GetFitT())));

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

            dTdLMin = std::min(dTdLMin, static_cast<float>(layerFitResult.GetGradient()));
            dTdLMax = std::max(dTdLMax, static_cast<float>(layerFitResult.GetGradient()));
            previousFitPosition = thisFitPosition;
        }

        if (diffWithStraightLineVector.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        diffWithStraightLineMean = 0.f;
        diffWithStraightLineSigma = 0.f;

        for (const float diffWithStraightLine : diffWithStraightLineVector)
            diffWithStraightLineMean += diffWithStraightLine;

        diffWithStraightLineMean /= static_cast<float>(diffWithStraightLineVector.size());

        for (const float diffWithStraightLine : diffWithStraightLineVector)
            diffWithStraightLineSigma += (diffWithStraightLine - diffWithStraightLineMean) * (diffWithStraightLine - diffWithStraightLineMean);

        if (diffWithStraightLineSigma < 0.f)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        diffWithStraightLineSigma = std::sqrt(diffWithStraightLineSigma / static_cast<float>(diffWithStraightLineVector.size()));
        dTdLWidth = dTdLMax - dTdLMin;
    }
    catch (const StatusCodeException &)
    {
        straightLineLengthLarge = -1.f;
        diffWithStraightLineMean = -1.f;
        diffWithStraightLineSigma = -1.f;
        dTdLWidth = -1.f;
        maxFitGapLength = -1.f;
        rmsSlidingLinearFit = -1.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDLinearFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoDVertexDistanceFeatureTool::TwoDVertexDistanceFeatureTool() : m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDVertexDistanceFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float straightLineLength(-1.f), ratio(-1.f);
    try
    {
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        straightLineLength = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        if (straightLineLength > std::numeric_limits<float>::epsilon())
            ratio = (CutClusterCharacterisationAlgorithm::GetVertexDistance(pAlgorithm, pCluster)) / straightLineLength;
    }
    catch (const StatusCodeException &)
    {
        ratio = -1.f;
    }
    featureVector.push_back(ratio);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDVertexDistanceFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PfoHierarchyFeatureTool::PfoHierarchyFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    CaloHitList parent3DHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, parent3DHitList);
    const unsigned int nParentHits3D(parent3DHitList.size());

    PfoList allDaughtersPfoList;
    LArPfoHelper::GetAllDownstreamPfos(pInputPfo, allDaughtersPfoList);
    const unsigned int nDaughterPfos(allDaughtersPfoList.empty() ? 0 : allDaughtersPfoList.size() - 1);

    unsigned int nDaughterHits3DTotal(0);

    if (nDaughterPfos > 0)
    {
        // ATTN This relies on knowing that the first pfo in allDaughtersPfoList is the input pfo
        allDaughtersPfoList.pop_front();

        for (const ParticleFlowObject *const pDaughterPfo : allDaughtersPfoList)
        {
            CaloHitList daughter3DHitList;
            LArPfoHelper::GetCaloHits(pDaughterPfo, TPC_3D, daughter3DHitList);
            nDaughterHits3DTotal += daughter3DHitList.size();
        }
    }

    const LArMvaHelper::MvaFeature nDaughters(static_cast<double>(nDaughterPfos));
    const LArMvaHelper::MvaFeature nDaughterHits3D(static_cast<double>(nDaughterHits3DTotal));
    const LArMvaHelper::MvaFeature daughterParentNHitsRatio(
        (nParentHits3D > 0) ? static_cast<double>(nDaughterHits3DTotal) / static_cast<double>(nParentHits3D) : 0.);

    featureVector.push_back(nDaughters);
    featureVector.push_back(nDaughterHits3D);
    featureVector.push_back(daughterParentNHitsRatio);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoHierarchyFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDLinearFitFeatureTool::ThreeDLinearFitFeatureTool() : m_slidingLinearFitWindow(3), m_slidingLinearFitWindowLarge(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, clusterList);
    float diffWithStraightLineMean(0.f), maxFitGapLength(0.f), rmsSlidingLinearFit(0.f);
    LArMvaHelper::MvaFeature length, diff, gap, rms;
    unsigned int nClustersUsed(0);

    for (const Cluster *const pCluster : clusterList)
    {
        float straightLineLengthLargeCluster(-1.f), diffWithStraightLineMeanCluster(-1.f), maxFitGapLengthCluster(-1.f),
            rmsSlidingLinearFitCluster(-1.f);

        this->CalculateVariablesSlidingLinearFit(
            pCluster, straightLineLengthLargeCluster, diffWithStraightLineMeanCluster, maxFitGapLengthCluster, rmsSlidingLinearFitCluster);

        if (straightLineLengthLargeCluster > std::numeric_limits<float>::epsilon())
        {
            diffWithStraightLineMeanCluster /= straightLineLengthLargeCluster;
            maxFitGapLengthCluster /= straightLineLengthLargeCluster;
            rmsSlidingLinearFitCluster /= straightLineLengthLargeCluster;

            diffWithStraightLineMean += diffWithStraightLineMeanCluster;
            maxFitGapLength += maxFitGapLengthCluster;
            rmsSlidingLinearFit += rmsSlidingLinearFitCluster;

            ++nClustersUsed;
        }
    }

    if (nClustersUsed > 0)
    {
        const float nClusters(static_cast<float>(nClustersUsed));
        length = std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pInputPfo));
        diff = diffWithStraightLineMean / nClusters;
        gap = maxFitGapLength / nClusters;
        rms = rmsSlidingLinearFit / nClusters;
    }

    featureVector.push_back(length);
    featureVector.push_back(diff);
    featureVector.push_back(gap);
    featureVector.push_back(rms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::RunWithMap(LArMvaHelper::MvaFeatureMap &featureMap, LArMvaHelper::MvaFeatureVector &featureVector, std::string featureToolName,
					    const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run( toolFeatureVec, pAlgorithm, pInputPfo );

    if ( featureMap.find(featureToolName+"_Length") != featureMap.end() ||
	 featureMap.find(featureToolName+"_DiffStraightLineMean") != featureMap.end() ||
	 featureMap.find(featureToolName+"_MaxFitGapLength") != featureMap.end() ||
         featureMap.find(featureToolName+"_SlidingLinearFitRMS") != featureMap.end() ){
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
	return;
    }

    featureMap[ featureToolName+"_Length" ] = toolFeatureVec[0].Get();
    featureMap[ featureToolName+"_DiffStraightLineMean" ] = toolFeatureVec[1].Get();
    featureMap[ featureToolName+"_MaxFitGapLength" ] = toolFeatureVec[2].Get();
    featureMap[ featureToolName+"_SlidingLinearFitRMS" ] = toolFeatureVec[3].Get();

    for ( auto const& feature : toolFeatureVec )
      featureVector.push_back(feature);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge,
    float &diffWithStraightLineMean, float &maxFitGapLength, float &rmsSlidingLinearFit) const
{
    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResultLarge(
            pCluster, m_slidingLinearFitWindowLarge, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        if (slidingFitResult.GetLayerFitResultMap().empty())
            throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

        straightLineLengthLarge =
            (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        rmsSlidingLinearFit = 0.f;

        FloatVector diffWithStraightLineVector;
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        CartesianVector previousFitPosition(slidingFitResult.GetGlobalMinLayerPosition());
        float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());

        for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
        {
            const LayerFitResult &layerFitResult(mapEntry.second);
            rmsSlidingLinearFit += layerFitResult.GetRms();

            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(layerFitResult.GetL(), layerFitResult.GetFitT(), thisFitPosition);

            LayerFitResultMap::const_iterator iterLarge =
                slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(layerFitResult.GetL()));

            if (slidingFitResultLarge.GetLayerFitResultMap().end() == iterLarge)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            diffWithStraightLineVector.push_back(static_cast<float>(std::fabs(layerFitResult.GetFitT() - iterLarge->second.GetFitT())));

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
            else
            {
                maxFitGapLength = 0.f;
            }

            dTdLMin = std::min(dTdLMin, static_cast<float>(layerFitResult.GetGradient()));
            dTdLMax = std::max(dTdLMax, static_cast<float>(layerFitResult.GetGradient()));
            previousFitPosition = thisFitPosition;
        }

        if (diffWithStraightLineVector.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        diffWithStraightLineMean = 0.f;

        for (const float diffWithStraightLine : diffWithStraightLineVector)
            diffWithStraightLineMean += diffWithStraightLine;

        diffWithStraightLineMean /= static_cast<float>(diffWithStraightLineVector.size());
    }
    catch (const StatusCodeException &)
    {
        straightLineLengthLarge = -1.f;
        diffWithStraightLineMean = -1.f;
        maxFitGapLength = -1.f;
        rmsSlidingLinearFit = -1.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDLinearFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDVertexDistanceFeatureTool::ThreeDVertexDistanceFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexDistanceFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    LArMvaHelper::MvaFeature vertexDistance;

    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || pVertexList->empty())
    {
        featureVector.push_back(vertexDistance);
        return;
    }

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }

    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        try
        {
            vertexDistance = (pInteractionVertex->GetPosition() - LArPfoHelper::GetVertex(pInputPfo)->GetPosition()).GetMagnitude();
        }
        catch (const StatusCodeException &)
        {
            CaloHitList threeDCaloHitList;
            LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

            if (!threeDCaloHitList.empty())
                vertexDistance = (pInteractionVertex->GetPosition() - (threeDCaloHitList.front())->GetPositionVector()).GetMagnitude();
        }
    }

    featureVector.push_back(vertexDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexDistanceFeatureTool::RunWithMap(LArMvaHelper::MvaFeatureMap &featureMap, LArMvaHelper::MvaFeatureVector &featureVector, std::string featureToolName,
					  const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run( toolFeatureVec, pAlgorithm, pInputPfo );

    if ( featureMap.find(featureToolName+"_VertexDistance")!=featureMap.end() ) {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
	 return;
    }

    featureMap[ featureToolName+"_VertexDistance" ] = toolFeatureVec[0].Get();

    for ( auto const& feature : toolFeatureVec )
      featureVector.push_back(feature);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDVertexDistanceFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDOpeningAngleFeatureTool::ThreeDOpeningAngleFeatureTool() : m_hitFraction(0.5f), m_defaultValue(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Need the 3D hits to calculate PCA components
    CaloHitList threeDCaloHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

    LArMvaHelper::MvaFeature diffAngle;
    if (!threeDCaloHitList.empty())
    {
        CartesianPointVector pointVectorStart, pointVectorEnd;
        this->Divide3DCaloHitList(pAlgorithm, threeDCaloHitList, pointVectorStart, pointVectorEnd);

        // Able to calculate angles only if > 1 point provided
        if ((pointVectorStart.size() > 1) && (pointVectorEnd.size() > 1))
        {
            try
            {
                // Run the PCA analysis twice
                CartesianVector centroidStart(0.f, 0.f, 0.f), centroidEnd(0.f, 0.f, 0.f);
                LArPcaHelper::EigenVectors eigenVecsStart, eigenVecsEnd;
                LArPcaHelper::EigenValues eigenValuesStart(0.f, 0.f, 0.f), eigenValuesEnd(0.f, 0.f, 0.f);

                LArPcaHelper::RunPca(pointVectorStart, centroidStart, eigenValuesStart, eigenVecsStart);
                LArPcaHelper::RunPca(pointVectorEnd, centroidEnd, eigenValuesEnd, eigenVecsEnd);

                const float openingAngle(this->OpeningAngle(eigenVecsStart.at(0), eigenVecsStart.at(1), eigenValuesStart));
                const float closingAngle(this->OpeningAngle(eigenVecsEnd.at(0), eigenVecsEnd.at(1), eigenValuesEnd));
                diffAngle = std::fabs(openingAngle - closingAngle);
            }
            catch (const StatusCodeException &)
            {
            }
        }
        else
        {
            diffAngle = m_defaultValue;
        }
    }

    featureVector.push_back(diffAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::RunWithMap(LArMvaHelper::MvaFeatureMap &featureMap, LArMvaHelper::MvaFeatureVector &featureVector, std::string featureToolName,
					const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run( toolFeatureVec, pAlgorithm, pInputPfo );

    if ( featureMap.find(featureToolName+"_AngleDiff")!=featureMap.end() ) {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
	return;
    }

    featureMap[ featureToolName+"_AngleDiff" ] = toolFeatureVec[0].Get();

    for ( auto const& feature : toolFeatureVec )
      featureVector.push_back(feature);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::Divide3DCaloHitList(const Algorithm *const pAlgorithm, const CaloHitList &threeDCaloHitList,
    CartesianPointVector &pointVectorStart, CartesianPointVector &pointVectorEnd)
{
    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (threeDCaloHitList.empty() || !pVertexList || pVertexList->empty())
        return;

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }

    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        // Order by distance to vertex, so first ones are closer to nuvertex
        CaloHitVector threeDCaloHitVector(threeDCaloHitList.begin(), threeDCaloHitList.end());
        std::sort(threeDCaloHitVector.begin(), threeDCaloHitVector.end(),
            ThreeDChargeFeatureTool::VertexComparator(pInteractionVertex->GetPosition()));

        unsigned int iHit(1);
        const unsigned int nHits(threeDCaloHitVector.size());

        for (const CaloHit *const pCaloHit : threeDCaloHitVector)
        {
            if (static_cast<float>(iHit) / static_cast<float>(nHits) <= m_hitFraction)
                pointVectorStart.push_back(pCaloHit->GetPositionVector());

            if (static_cast<float>(iHit) / static_cast<float>(nHits) >= 1.f - m_hitFraction)
                pointVectorEnd.push_back(pCaloHit->GetPositionVector());

            ++iHit;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDOpeningAngleFeatureTool::OpeningAngle(const CartesianVector &principal, const CartesianVector &secondary, const CartesianVector &eigenValues) const
{
    const float principalMagnitude(principal.GetMagnitude());
    const float secondaryMagnitude(secondary.GetMagnitude());

    if (std::fabs(principalMagnitude) < std::numeric_limits<float>::epsilon())
    {
        std::cout << "ThreeDOpeningAngleFeatureTool::OpeningAngle - The principal eigenvector is 0." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
    else if (std::fabs(secondaryMagnitude) < std::numeric_limits<float>::epsilon())
    {
        return 0.f;
    }

    const float cosTheta(principal.GetDotProduct(secondary) / (principalMagnitude * secondaryMagnitude));

    if (cosTheta > 1.f)
    {
        std::cout << "PcaShowerParticleBuildingAlgorithm::OpeningAngle - cos(theta) reportedly greater than 1." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const float sinTheta(std::sqrt(1.f - cosTheta * cosTheta));

    if (eigenValues.GetX() < std::numeric_limits<float>::epsilon())
    {
        std::cout << "PcaShowerParticleBuildingAlgorithm::OpeningAngle - principal eigenvalue less than or equal to 0." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
    else if (eigenValues.GetY() < std::numeric_limits<float>::epsilon())
    {
        return 0.f;
    }

    return std::atan(std::sqrt(eigenValues.GetY()) * sinTheta / std::sqrt(eigenValues.GetX()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDOpeningAngleFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitFraction", m_hitFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DefaultValue", m_defaultValue));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDPCAFeatureTool::ThreeDPCAFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDPCAFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    LArMvaHelper::MvaFeature pca1, pca2;

    // Need the 3D hits to calculate PCA components
    CaloHitList threeDCaloHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

    if (!threeDCaloHitList.empty())
    {
        try
        {
            CartesianVector centroid(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVecs;
            LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);

            LArPcaHelper::RunPca(threeDCaloHitList, centroid, eigenValues, eigenVecs);
            const float principalEigenvalue(eigenValues.GetX()), secondaryEigenvalue(eigenValues.GetY()), tertiaryEigenvalue(eigenValues.GetZ());

            if (principalEigenvalue > std::numeric_limits<float>::epsilon())
            {
                pca1 = secondaryEigenvalue / principalEigenvalue;
                pca2 = tertiaryEigenvalue / principalEigenvalue;
            }
            else
            {
                // ATTN if n3dHits == 1 then principal, secondary, and tertiary eigenvalues are zero hence default to zero
                pca1 = 0.;
                pca2 = 0.;
            }
        }
        catch (const StatusCodeException &)
        {
        }
    }

    featureVector.push_back(pca1);
    featureVector.push_back(pca2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDPCAFeatureTool::RunWithMap(LArMvaHelper::MvaFeatureMap &featureMap, LArMvaHelper::MvaFeatureVector &featureVector, std::string featureToolName,
				  const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run( toolFeatureVec, pAlgorithm, pInputPfo );

    if ( featureMap.find(featureToolName+"_SecondaryPCARatio")!=featureMap.end() ||
	 featureMap.find(featureToolName+"_TertiaryPCARatio")!=featureMap.end() ) {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
	return;
    }

    featureMap[ featureToolName+"_SecondaryPCARatio" ] = toolFeatureVec[0].Get();
    featureMap[ featureToolName+"_TertiaryPCARatio" ] = toolFeatureVec[1].Get();

    for ( auto const& feature : toolFeatureVec )
      featureVector.push_back(feature);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDPCAFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDChargeFeatureTool::ThreeDChargeFeatureTool() : m_endChargeFraction(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float totalCharge(-1.f), chargeSigma(-1.f), chargeMean(-1.f), endCharge(-1.f);
    LArMvaHelper::MvaFeature charge1, charge2;

    ClusterList clusterListW;
    LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, clusterListW);

    if (!clusterListW.empty())
        this->CalculateChargeVariables(pAlgorithm, clusterListW.front(), totalCharge, chargeSigma, chargeMean, endCharge);

    if (chargeMean > std::numeric_limits<float>::epsilon())
        charge1 = chargeSigma / chargeMean;

    if (totalCharge > std::numeric_limits<float>::epsilon())
        charge2 = endCharge / totalCharge;

    featureVector.push_back(charge1);
    featureVector.push_back(charge2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::RunWithMap(LArMvaHelper::MvaFeatureMap &featureMap, LArMvaHelper::MvaFeatureVector &featureVector, std::string featureToolName,
				  const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
  LArMvaHelper::MvaFeatureVector toolFeatureVec;
  this->Run( toolFeatureVec, pAlgorithm, pInputPfo );

  if ( featureMap.find(featureToolName+"_FractionalSpread")!=featureMap.end() ||
       featureMap.find(featureToolName+"_EndFraction")!=featureMap.end() ) {
      std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
      return;
  }

  featureMap[ featureToolName+"_FractionalSpread" ] = toolFeatureVec[0].Get();
  featureMap[ featureToolName+"_EndFraction" ] = toolFeatureVec[1].Get();

  for ( auto const& feature : toolFeatureVec )
    featureVector.push_back(feature);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::CalculateChargeVariables(const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster,
    float &totalCharge, float &chargeSigma, float &chargeMean, float &endCharge)
{
    totalCharge = 0.f;
    chargeSigma = 0.f;
    chargeMean = 0.f;
    endCharge = 0.f;

    CaloHitList orderedCaloHitList;
    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, pCluster, orderedCaloHitList);

    FloatVector chargeVector;
    unsigned int hitCounter(0);
    const unsigned int nTotalHits(orderedCaloHitList.size());

    for (const CaloHit *const pCaloHit : orderedCaloHitList)
    {
        ++hitCounter;
        const float pCaloHitCharge(pCaloHit->GetInputEnergy());

        if (pCaloHitCharge >= 0.f)
        {
            totalCharge += pCaloHitCharge;
            chargeVector.push_back(pCaloHitCharge);

            if (hitCounter >= std::floor(static_cast<float>(nTotalHits) * (1.f - m_endChargeFraction)))
                endCharge += pCaloHitCharge;
        }
    }

    if (!chargeVector.empty())
    {
        chargeMean = totalCharge / static_cast<float>(chargeVector.size());

        for (const float charge : chargeVector)
            chargeSigma += (charge - chargeMean) * (charge - chargeMean);

        chargeSigma = std::sqrt(chargeSigma / static_cast<float>(chargeVector.size()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::OrderCaloHitsByDistanceToVertex(
    const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, CaloHitList &caloHitList)
{
    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || pVertexList->empty())
        return;

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }

    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pInteractionVertex->GetPosition(), hitType));

        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        clusterCaloHitList.sort(ThreeDChargeFeatureTool::VertexComparator(vertexPosition2D));
        caloHitList.insert(caloHitList.end(), clusterCaloHitList.begin(), clusterCaloHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDChargeFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EndChargeFraction", m_endChargeFraction));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDChargeFeatureTool::VertexComparator::VertexComparator(const CartesianVector vertexPosition2D) : m_neutrinoVertex(vertexPosition2D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDChargeFeatureTool::VertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
    const float distanceL((left->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    const float distanceR((right->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    return distanceL < distanceR;
}

} // namespace lar_content
