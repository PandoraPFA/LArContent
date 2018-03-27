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
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoDShowerFitFeatureTool::TwoDShowerFitFeatureTool() :
    m_slidingShowerFitWindow(3),
    m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDShowerFitFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float ratio(-1.f);
    try
    {
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const float straightLineLength = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        if (straightLineLength > std::numeric_limits<float>::epsilon())
            ratio = (CutClusterCharacterisationAlgorithm::GetShowerFitWidth(pAlgorithm, pCluster, m_slidingShowerFitWindow))/straightLineLength;
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingShowerFitWindow", m_slidingShowerFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoDLinearFitFeatureTool::TwoDLinearFitFeatureTool() :
    m_slidingLinearFitWindow(3),
    m_slidingLinearFitWindowLarge(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDLinearFitFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
const pandora::Cluster * const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float dTdLWidth(-1.f), straightLineLengthLarge(-1.f), diffWithStraightLineMean(-1.f), diffWithStraightLineSigma(-1.f), maxFitGapLength(-1.f), rmsSlidingLinearFit(-1.f);
    this->CalculateVariablesSlidingLinearFit(pCluster, straightLineLengthLarge, diffWithStraightLineMean, diffWithStraightLineSigma, dTdLWidth, maxFitGapLength, rmsSlidingLinearFit);

    if (straightLineLengthLarge > std::numeric_limits<float>::epsilon())
    {
        diffWithStraightLineMean  /= straightLineLengthLarge;
        diffWithStraightLineSigma /= straightLineLengthLarge;
        dTdLWidth                 /= straightLineLengthLarge;
        maxFitGapLength           /= straightLineLengthLarge;
        rmsSlidingLinearFit       /= straightLineLengthLarge;
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
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindowLarge, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        if (slidingFitResult.GetLayerFitResultMap().empty())
            throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

        straightLineLengthLarge = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
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

            LayerFitResultMap::const_iterator iterLarge = slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(layerFitResult.GetL()));

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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoDVertexDistanceFeatureTool::TwoDVertexDistanceFeatureTool() :
    m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDVertexDistanceFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float straightLineLength(-1.f), ratio(-1.f);
    try
    {
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        straightLineLength = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        if (straightLineLength > std::numeric_limits<float>::epsilon())
            ratio = (CutClusterCharacterisationAlgorithm::GetVertexDistance(pAlgorithm, pCluster))/straightLineLength;
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDLinearFitFeatureTool::ThreeDLinearFitFeatureTool() :
    m_slidingLinearFitWindow(3),
    m_slidingLinearFitWindowLarge(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, clusterList);
    float diffWithStraightLineMean(0.f), maxFitGapLength(0.f), rmsSlidingLinearFit(0.f);

    for (const Cluster *const pCluster : clusterList)
    {
        float straightLineLengthLargeCluster(-1.f), diffWithStraightLineMeanCluster(-1.f), maxFitGapLengthCluster(-1.f), rmsSlidingLinearFitCluster(-1.f);

        this->CalculateVariablesSlidingLinearFit(pCluster, straightLineLengthLargeCluster, diffWithStraightLineMeanCluster, maxFitGapLengthCluster, rmsSlidingLinearFitCluster);

        if (straightLineLengthLargeCluster > std::numeric_limits<float>::epsilon())
        {
            diffWithStraightLineMeanCluster  /= straightLineLengthLargeCluster;
            maxFitGapLengthCluster           /= straightLineLengthLargeCluster;
            rmsSlidingLinearFitCluster       /= straightLineLengthLargeCluster;
        }

        diffWithStraightLineMean   += diffWithStraightLineMeanCluster;
        maxFitGapLength            += maxFitGapLengthCluster;
        rmsSlidingLinearFit        += rmsSlidingLinearFitCluster;
    }

    if (!clusterList.empty())
    {
        const float nClusters(static_cast<float>(clusterList.size()));  
        diffWithStraightLineMean   /= nClusters;
        maxFitGapLength            /= nClusters;
        rmsSlidingLinearFit        /= nClusters;
    }

    const float lenghtSquare(LArPfoHelper::GetThreeDLengthSquared(pInputPfo));
    featureVector.push_back(std::sqrt(lenghtSquare));
    featureVector.push_back(diffWithStraightLineMean);
    featureVector.push_back(maxFitGapLength);
    featureVector.push_back(rmsSlidingLinearFit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge,
    float &diffWithStraightLineMean, float &maxFitGapLength, float &rmsSlidingLinearFit) const
{
    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindowLarge, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        if (slidingFitResult.GetLayerFitResultMap().empty())
            throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

        straightLineLengthLarge = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
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

            LayerFitResultMap::const_iterator iterLarge = slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(layerFitResult.GetL()));

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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDVertexDistanceFeatureTool::ThreeDVertexDistanceFeatureTool() :
    m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexDistanceFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    //find the neutrino vertex
    const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (pVertexList->empty() || (pVertexList->size() != 1) || (VERTEX_3D != pVertexList->front()->GetVertexType()))
    {
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    const Vertex *const nuVertex(pVertexList->front());
    //find the particle vertex
   float vertexDistance(-1.f);
    try
    {
        const Vertex *const pVertex = LArPfoHelper::GetVertex(pInputPfo);
        const CartesianVector nuPosition(nuVertex->GetPosition()), pfoPosition(pVertex->GetPosition());
        vertexDistance = (nuPosition - pfoPosition).GetMagnitude();
    }
    catch (const StatusCodeException &)
    {
        vertexDistance = -1.f;
    }

    featureVector.push_back(vertexDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDVertexDistanceFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDOpeningAngleFeatureTool::ThreeDOpeningAngleFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Need the 3D clusters and hits to calculate PCA components
    ClusterList threeDClusterList;
    LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);

    CaloHitList threeDCaloHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

    if (threeDCaloHitList.empty())
        return;

    CartesianPointVector pointVectorStart, pointVectorEnd;

    this->Divide3DCaloHitList(pAlgorithm, threeDCaloHitList, pointVectorStart, pointVectorEnd);
    // Run the PCA analysis twice
    CartesianVector centroidStart(0.f, 0.f, 0.f), centroidEnd(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecsStart, eigenVecsEnd;
    LArPcaHelper::EigenValues eigenValuesStart(0.f, 0.f, 0.f), eigenValuesEnd(0.f, 0.f, 0.f);

    LArPcaHelper::RunPca(pointVectorStart, centroidStart, eigenValuesStart, eigenVecsStart);
    LArPcaHelper::RunPca(pointVectorEnd, centroidEnd, eigenValuesEnd, eigenVecsEnd);

    const float openingAngle(this->OpeningAngle(eigenVecsStart.at(0), eigenVecsStart.at(1), eigenValuesStart));
    const float closingAngle(this->OpeningAngle(eigenVecsEnd.at(0), eigenVecsEnd.at(1), eigenValuesEnd));

    featureVector.push_back(std::fabs(openingAngle-closingAngle));

}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::Divide3DCaloHitList(const Algorithm *const pAlgorithm, CaloHitList &threeDCaloHitList, CartesianPointVector &pointVectorStart, CartesianPointVector &pointVectorEnd)
{
    const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || (pVertexList->size() != 1) || (VERTEX_3D != pVertexList->front()->GetVertexType()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const CartesianVector nuVertex(pVertexList->front()->GetPosition());
    CaloHitVector threeDCaloHitVector(threeDCaloHitList.begin(), threeDCaloHitList.end());

    //order by distance to vertex, so first ones are closer to nuvertex
    std::sort(threeDCaloHitVector.begin(), threeDCaloHitVector.end(), ThreeDChargeFeatureTool::VertexComparator(nuVertex));
    CaloHitList orderedCaloHitList(threeDCaloHitVector.begin(),threeDCaloHitVector.end());
    const unsigned int nhits(orderedCaloHitList.size());

    for (const CaloHit *const pCaloHit : orderedCaloHitList)
    {
        CartesianPointVector &targetVector((pointVectorStart.size() < (nhits / 2)) ? pointVectorStart : pointVectorEnd);
        targetVector.push_back(pCaloHit->GetPositionVector());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDOpeningAngleFeatureTool::OpeningAngle(const CartesianVector &principal, const CartesianVector &secondary,
    const CartesianVector &eigenValues) const
{
    const float principalMagnitude(principal.GetMagnitude());
    const float secondaryMagnitude(secondary.GetMagnitude());

    if (std::fabs(principalMagnitude) < std::numeric_limits<float>::epsilon())
    {
        std::cout << "PcaShowerParticleBuildingAlgorithm::OpeningAngle - The principal eigenvector is 0." << std::endl;
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

    if (std::fabs(eigenValues.GetX()) < std::numeric_limits<float>::epsilon())
    {
        std::cout << "PcaShowerParticleBuildingAlgorithm::OpeningAngle - principal eigenvalue less than or equal to 0." << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    }
    else if (std::fabs(eigenValues.GetY()) < std::numeric_limits<float>::epsilon())
    {
        return 0.f;
    }

    return std::atan(std::sqrt(eigenValues.GetY()) * sinTheta / std::sqrt(eigenValues.GetX()));
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDOpeningAngleFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDPCAFeatureTool::ThreeDPCAFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDPCAFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pInputPfo)
{
     if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Need the 3D cluster and hits to calculate PCA components
    ClusterList threeDClusterList;
    LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);

    if (threeDClusterList.empty())
        return;

    CaloHitList threeDCaloHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

    if (threeDCaloHitList.empty())
        return;

    // Run the PCA analysis
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(threeDCaloHitList, centroid, eigenValues, eigenVecs);
    const float principalEigenvalue(eigenValues.GetX()), secondaryEigenvalue(eigenValues.GetY()), tertiaryEigenvalue(eigenValues.GetZ());

    float pca1(-1.f), pca2(-1.f);
    if (principalEigenvalue > std::numeric_limits<float>::epsilon())
    {
        pca1 = secondaryEigenvalue/principalEigenvalue;
        pca2 = tertiaryEigenvalue/principalEigenvalue;
    }

    featureVector.push_back(pca1);
    featureVector.push_back(pca2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDPCAFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
ThreeDChargeFeatureTool::ThreeDChargeFeatureTool() :
    m_endChargeFraction(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float totalCharge(-1.f), chargeSigma(-1.f), chargeMean(-1.f), endCharge(-1.f);
    ClusterList pClusterList;
    LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, pClusterList);
    if ((!pClusterList.empty()) && (pClusterList.size() == 1))
    {
        const Cluster *const pCluster(pClusterList.front());
        this->CalculateChargeVariables(pAlgorithm, pCluster, totalCharge, chargeSigma, chargeMean, endCharge, m_endChargeFraction);
    }

    if (chargeMean > std::numeric_limits<float>::epsilon())
        chargeSigma /= chargeMean;

    if (totalCharge > std::numeric_limits<float>::epsilon())
        endCharge /= totalCharge;

    featureVector.push_back(chargeSigma);
    featureVector.push_back(endCharge);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::CalculateChargeVariables(const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, float &totalCharge,
    float &chargeSigma, float &chargeMean, float &endCharge, float endChargeFraction)
{

    CaloHitList orderedCaloHitList;
    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, pCluster, orderedCaloHitList);

    const int totalHits(pCluster->GetNCaloHits());
    FloatVector chargeVector;
    int hitCounter(0);
    totalCharge = 0.f;

    for (const CaloHit *const pCaloHit : orderedCaloHitList)
    {
        hitCounter++;
        const float pCaloHitCharge(pCaloHit->GetInputEnergy());

        if (pCaloHitCharge < 0)
        {
            std::cout << "Found a hit with negative charge! " << std::endl;
        }
        else
        {
            totalCharge    += pCaloHitCharge;
            chargeVector.push_back(pCaloHitCharge);

            if (hitCounter >= std::floor(totalHits*(1.f-endChargeFraction)))
            {
                endCharge += pCaloHitCharge;
            }
        }
    }

    chargeMean = 0.f;
    chargeSigma = 0.f;

    for (const float charge : chargeVector)
        chargeMean += charge;

    chargeMean /= static_cast<float>(chargeVector.size());

    for (const float charge : chargeVector)
        chargeSigma += (charge - chargeMean) * (charge - chargeMean);

    if (chargeSigma < 0.f)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    chargeSigma = std::sqrt(chargeSigma / static_cast<float>(chargeVector.size()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::OrderCaloHitsByDistanceToVertex(const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, CaloHitList &caloHitList)
{
    //find the neutrino vertex and sort hits by distance to vertex
    const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || (pVertexList->size() != 1) || (VERTEX_3D != pVertexList->front()->GetVertexType()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const Vertex *const pVertex(pVertexList->front());
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pVertex->GetPosition(), hitType));

    CaloHitList clusterCaloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);
    CaloHitVector clusterCaloHitVector(clusterCaloHitList.begin(), clusterCaloHitList.end());

    //TODO: might give problems if vertex in the middle of the cluster ?
    std::sort(clusterCaloHitVector.begin(), clusterCaloHitVector.end(), VertexComparator(vertexPosition2D));
    caloHitList.insert(caloHitList.end(), clusterCaloHitVector.begin(), clusterCaloHitVector.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDChargeFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EndChargeFraction", m_endChargeFraction));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDChargeFeatureTool::VertexComparator::VertexComparator(const CartesianVector vertexPosition2D) :
    m_neutrinoVertex(vertexPosition2D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDChargeFeatureTool::VertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
    float distanceL((left->GetPositionVector()-m_neutrinoVertex).GetMagnitudeSquared());
    float distanceR((right->GetPositionVector()-m_neutrinoVertex).GetMagnitudeSquared());
    return distanceL < distanceR;
}

} // namespace lar_content
