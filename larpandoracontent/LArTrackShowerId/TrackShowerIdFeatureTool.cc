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
#include <vector>
#include <list>
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

//-------------------------------------------------------------------------------------------------------------------------------
PfoHierarchyFeatureTool::PfoHierarchyFeatureTool() 
{
}

//-------------------------------------------------------------------------------------------------------------------------------
void PfoHierarchyFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
  	if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
    	std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
  
    LArMvaHelper::MvaFeature nAllDaughter, nHits3DDaughterTotal, daughterParentNhitsRatio;

    CaloHitList nHits3DParentList;    
    size_t nHits3DDaughter(0);
    PfoList newPfoList(1, pInputPfo);
    size_t nHits3DParent(1.f);
    PfoList allDaughtersPfoList;
    float nHits3DDaughterTotalNumber(0);

    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, nHits3DParentList);
    nHits3DParent = nHits3DParentList.size();
    LArPfoHelper::GetAllDownstreamPfos(newPfoList, allDaughtersPfoList);
    nAllDaughter = allDaughtersPfoList.size() - 1;
    if (nAllDaughter.Get() > 0.0)
    {
    	allDaughtersPfoList.pop_front();
		for(const ParticleFlowObject *const pDaughterPfo : allDaughtersPfoList)
			{
        		CaloHitList nHits3DDaughterList;
               	LArPfoHelper::GetCaloHits(pDaughterPfo, TPC_3D, nHits3DDaughterList);
                nHits3DDaughter = nHits3DDaughterList.size();
                nHits3DDaughterTotalNumber += nHits3DDaughter;
			}
    }
    else if (nAllDaughter.Get() == 0.0)
	{
     	nHits3DDaughter = 0.0;
    }

    nHits3DDaughterTotal = nHits3DDaughterTotalNumber;
    daughterParentNhitsRatio = ((nHits3DDaughterTotal.Get()))/(static_cast<double>(nHits3DParent));
  
  //---------------push_back into feature vector-----------------------------------------------------------------------------
	featureVector.push_back(nAllDaughter);
	featureVector.push_back(nHits3DDaughterTotal);
	featureVector.push_back(daughterParentNhitsRatio);
}

StatusCode PfoHierarchyFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{


    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
ThreeDLinearFitFeatureTool::ThreeDLinearFitFeatureTool() :
    m_slidingLinearFitWindow(3),
    m_slidingLinearFitWindowLarge(10000)
{
}
//------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
const pandora::ParticleFlowObject *const pInputPfo)
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
        float straightLineLengthLargeCluster(-1.f), diffWithStraightLineMeanCluster(-1.f), maxFitGapLengthCluster(-1.f), rmsSlidingLinearFitCluster(-1.f);

        this->CalculateVariablesSlidingLinearFit(pCluster, straightLineLengthLargeCluster, diffWithStraightLineMeanCluster, maxFitGapLengthCluster, rmsSlidingLinearFitCluster);

        if (straightLineLengthLargeCluster > std::numeric_limits<float>::epsilon())
        {
            diffWithStraightLineMeanCluster  /= straightLineLengthLargeCluster;
            maxFitGapLengthCluster           /= straightLineLengthLargeCluster;
            rmsSlidingLinearFitCluster       /= straightLineLengthLargeCluster;

            diffWithStraightLineMean   += diffWithStraightLineMeanCluster;
            maxFitGapLength            += maxFitGapLengthCluster;
            rmsSlidingLinearFit        += rmsSlidingLinearFitCluster;

            ++nClustersUsed;
        }
    }

    if (nClustersUsed > 0)
    {
        const float nClusters(static_cast<float>(nClustersUsed));
        length = std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pInputPfo));
        diff   = diffWithStraightLineMean / nClusters;
        gap    = maxFitGapLength          / nClusters;
        rms    = rmsSlidingLinearFit      / nClusters;
    }

    featureVector.push_back(length);
    featureVector.push_back(diff);
    featureVector.push_back(gap);
    featureVector.push_back(rms);
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDVertexDistanceFeatureTool::ThreeDVertexDistanceFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexDistanceFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    LArMvaHelper::MvaFeature vertexDistance;
    const VertexList *pVertexList(nullptr);
    (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if ((!pVertexList->empty()) && (pVertexList->size() == 1) && (VERTEX_3D == pVertexList->front()->GetVertexType()))
    {
        try
        {
            vertexDistance = (pVertexList->front()->GetPosition() - LArPfoHelper::GetVertex(pInputPfo)->GetPosition()).GetMagnitude();
        }
        catch (const StatusCodeException &) {}
    }

    featureVector.push_back(vertexDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

  StatusCode ThreeDVertexDistanceFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDOpeningAngleFeatureTool::ThreeDOpeningAngleFeatureTool() :
  m_hitFraction(0.5)
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

    LArMvaHelper::MvaFeature diffAngle;
    if (!threeDCaloHitList.empty())
    {
        CartesianPointVector pointVectorStart, pointVectorEnd;
        this->Divide3DCaloHitList(pAlgorithm, threeDCaloHitList, pointVectorStart, pointVectorEnd);

        //able to calculate angles only if > 1 point provided
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
                diffAngle = std::fabs(openingAngle-closingAngle);
            }
            catch (const StatusCodeException &){}
        }
		else
		{
		throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
		}    

    }

    featureVector.push_back(diffAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::Divide3DCaloHitList(const Algorithm *const pAlgorithm, CaloHitList &threeDCaloHitList, CartesianPointVector &pointVectorStart, CartesianPointVector &pointVectorEnd)
{
    const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if ((!pVertexList->empty()) && (pVertexList->size() == 1) && (VERTEX_3D == pVertexList->front()->GetVertexType()))
    {
        const CartesianVector nuVertex(pVertexList->front()->GetPosition());
        CaloHitVector threeDCaloHitVector(threeDCaloHitList.begin(), threeDCaloHitList.end());

        //order by distance to vertex, so first ones are closer to nuvertex
        std::sort(threeDCaloHitVector.begin(), threeDCaloHitVector.end(), ThreeDChargeFeatureTool::VertexComparator(nuVertex));
        CaloHitList orderedCaloHitList(threeDCaloHitVector.begin(),threeDCaloHitVector.end());
        const unsigned int nhits(orderedCaloHitList.size());
		int iHit = 1;
		for (const CaloHit *const pCaloHit : orderedCaloHitList)
	  		{
	    		if((float)iHit / nhits <= m_hitFraction)
              		pointVectorStart.push_back(pCaloHit->GetPositionVector());
	    		if((float)iHit / nhits >= 1.0 - m_hitFraction)
              		pointVectorEnd.push_back(pCaloHit->GetPositionVector());
	    		iHit++;
	  		}
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

StatusCode ThreeDOpeningAngleFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                                                                                                       "HitFraction", m_hitFraction));

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

    LArMvaHelper::MvaFeature pca1, pca2;
    // Need the 3D cluster and hits to calculate PCA components
    ClusterList threeDClusterList;
    LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);

    CaloHitList threeDCaloHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

    if ((!threeDClusterList.empty()) && (!threeDCaloHitList.empty()))
    {
        // Run the PCA analysis
        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecs;
        LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        try
        {
            LArPcaHelper::RunPca(threeDCaloHitList, centroid, eigenValues, eigenVecs);
            const float principalEigenvalue(eigenValues.GetX()), secondaryEigenvalue(eigenValues.GetY()), tertiaryEigenvalue(eigenValues.GetZ());
            if (principalEigenvalue > std::numeric_limits<float>::epsilon())
            {
                pca1 = secondaryEigenvalue/principalEigenvalue;
                pca2 = tertiaryEigenvalue/principalEigenvalue;
            }
        }
        catch (const StatusCodeException &){}
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
ThreeDPCAVariablesFeatureTool::ThreeDPCAVariablesFeatureTool() :
  m_pfoFraction(0.2),
  m_minHits(3),
  m_segmentWidth(1.0),
  m_MoliereRadius(10.1),
  m_MoliereFraction(0.05),
  m_slidingLinearFitWindow(3) //added by mousam
  {
  }

//------------------------------------------------------------------------------------------------------------------------------------------
void ThreeDPCAVariablesFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo) 

{
  if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
    std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

  LArMvaHelper::MvaFeature concentration, concentration2, conicalness, invHitDensityY, invHitDensityZ, 
    invHitDensityYRatio, invHitDensityZRatio, nSegmentsDoubleHits, nSegmentsDoubleHitsRatio, nHits;
  float conc = 1.f;
  float conc2 = 1.f;
  float conic = 1.f;
  float densY = 0.f;
  float densZ = 0.f;
  float densYRatio = 0.f;
  float densZRatio = 0.f;
  int nSegDoubleHits = 0;
  float nSegDoubleHitsRatio = 0.f;
  nHits = 0;
  //Need the 3D cluster and hits to calculate PCA components                                                                                                                                           
  ClusterList threeDClusterList;
  LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);

  CaloHitList threeDCaloHitList;
  LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

  if ((!threeDClusterList.empty()) && (!threeDCaloHitList.empty()))
    {
      //Run the PCA analysis                                                                                                                                                                           
      CartesianVector centroid(0.f, 0.f, 0.f);
      LArPcaHelper::EigenVectors eigenVecs;
      LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        try
          {
	    LArPcaHelper::RunPca(threeDCaloHitList, centroid, eigenValues, eigenVecs);
          }
        catch (const StatusCodeException &){}
	
	std::map<float,Eigen::Vector3f> pcaPositions;

	//get matrix of eigenvectors and its transpose
	Eigen::Matrix3f eVecs;
	eVecs << eigenVecs.at(0).GetX(), eigenVecs.at(1).GetX(), eigenVecs.at(2).GetX(),
	  eigenVecs.at(0).GetY(), eigenVecs.at(1).GetY(), eigenVecs.at(2).GetY(),
	  eigenVecs.at(0).GetZ(), eigenVecs.at(1).GetZ(), eigenVecs.at(2).GetZ();

	Eigen::Matrix3f eVecsTranspose = eVecs.transpose();

	Eigen::Vector3f PcaCentroid(centroid.GetX(), centroid.GetY(), centroid.GetZ());

	for (const CaloHit *const pCaloHit : threeDCaloHitList)
	  {
	    //get spacepoint in Cartesian coordinates
	    const CartesianVector spacePt = pCaloHit->GetPositionVector();
	    Eigen::Vector3f spacePoint(spacePt.GetX(), spacePt.GetY(), spacePt.GetZ());

	    //calculate position of hit in PCA coordinates
	    Eigen::Vector3f PcaPoint = eVecsTranspose * (spacePoint - PcaCentroid);

	    //it might seem odd to put the PCA x positions in both key and value
	    //this is done to order hit PCA positions by x coordinate which is the principal axis
	    pcaPositions.insert(std::pair<float,Eigen::Vector3f>(PcaPoint.x(),PcaPoint));

	  }

	const float paStart = pcaPositions.begin()->first;
	const float paEnd = pcaPositions.rbegin()->first;
	const float pfoLengthOnPA = fabs(paEnd - paStart);

	this->CalculateConcentrationConicalness(threeDCaloHitList, eVecsTranspose, PcaCentroid, paStart, pfoLengthOnPA, conc, conc2, conic);
	
	this->CalculateDensity(paStart, pfoLengthOnPA, pcaPositions,
                               densY, densZ, densYRatio, densZRatio);

	this->SegmentsWithDoubleHits(paStart, pcaPositions, nSegDoubleHits, nSegDoubleHitsRatio);
    }

  //-------------------------------------------------------------------
  concentration = conc;
  concentration2 = conc2;
  conicalness = conic;
  invHitDensityY = densY;
  invHitDensityZ = densZ;
  invHitDensityYRatio = densYRatio;
  invHitDensityZRatio = densZRatio;
  nSegmentsDoubleHits = nSegDoubleHits;
  nSegmentsDoubleHitsRatio = nSegDoubleHitsRatio;
  nHits = threeDCaloHitList.size();

  featureVector.push_back(concentration);
  featureVector.push_back(concentration2);
  featureVector.push_back(conicalness);
  featureVector.push_back(invHitDensityY);
  featureVector.push_back(invHitDensityZ);
  featureVector.push_back(invHitDensityYRatio);
  featureVector.push_back(invHitDensityZRatio);
  featureVector.push_back(nSegmentsDoubleHits);
  featureVector.push_back(nSegmentsDoubleHitsRatio);
  featureVector.push_back(nHits); 
}

//------------------------------------------------------------------------------------------------------------------------------------------
void ThreeDPCAVariablesFeatureTool::CalculateConcentrationConicalness(CaloHitList &threeDCaloHitList, Eigen::Matrix3f &eVecsTranspose, Eigen::Vector3f &PcaCentroid,
								      const float paStart, const float pfoLengthOnPA, float &concentration, float &concentration2, float &conicalness)
{

  concentration = 1.f;
  concentration2 = 1.f;
  conicalness = 1.f;

  float distanceFromPA = 0.f;
  float totalCharge = 0.f;
  float chargeOverDist = 0.f;
  float chargeTimesDist = 0.f;
  float chargeStart = 0.f;
  float chargeEnd = 0.f;

  float conStart = 0.f;
  float conEnd = 0.f;
  int nHits = 0;
  int nHitsStart = 0;
  int nHitsEnd = 0;

  for (const CaloHit *const pCaloHit : threeDCaloHitList)
    {
      nHits++;
      //get spacepoint in Cartesian coordinates
      const CartesianVector spacePt = pCaloHit->GetPositionVector();
      Eigen::Vector3f spacePoint(spacePt.GetX(), spacePt.GetY(), spacePt.GetZ());

      //calculate position of hit in PCA coordinates
      Eigen::Vector3f PcaPoint = eVecsTranspose * (spacePoint - PcaCentroid);

      distanceFromPA = sqrt(PcaPoint.y() * PcaPoint.y() + PcaPoint.z() * PcaPoint.z());

      const float pCaloHitCharge(pCaloHit->GetInputEnergy());

      if(pCaloHitCharge >= 0)
	{
	  totalCharge += pCaloHitCharge;

	  if(distanceFromPA > std::numeric_limits<float>::epsilon())
	    {
	    chargeOverDist += pCaloHitCharge / distanceFromPA;
	    chargeTimesDist += pCaloHitCharge * distanceFromPA;
	    }

	  if(fabs((PcaPoint.x() - paStart) / pfoLengthOnPA) <= m_pfoFraction)
	    {
	      conStart += distanceFromPA * distanceFromPA * pCaloHitCharge;
	      chargeStart += pCaloHitCharge;
	      nHitsStart++;
	    }
	  if(fabs((PcaPoint.x() - paStart) / pfoLengthOnPA) >= 1.0 - m_pfoFraction)
	    {
	      conEnd += distanceFromPA * distanceFromPA * pCaloHitCharge;
	      chargeEnd += pCaloHitCharge;
	      nHitsEnd++;
	    }
	}
    }

  if(threeDCaloHitList.size() >= m_minHits && totalCharge > std::numeric_limits<float>::epsilon())
    {
    concentration = chargeOverDist / totalCharge;
    concentration2 = chargeTimesDist / totalCharge;
    }

  if(nHitsStart >= m_minHits && nHitsEnd >= m_minHits && chargeEnd > std::numeric_limits<float>::epsilon() && sqrt(conStart) / chargeStart > std::numeric_limits<float>::epsilon())
    conicalness = (sqrt(conEnd) / chargeEnd) / (sqrt(conStart) / chargeStart);
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

void ThreeDPCAVariablesFeatureTool::CalculateDensity(const float paStart, const float pfoLengthOnPA, std::map<float,Eigen::Vector3f> &pcaPositions, 
                                                       float &invHitDensityY, float &invHitDensityZ,
						       float &invHitDensityYRatio, float &invHitDensityZRatio)
{

  invHitDensityY = 0.f;
  invHitDensityZ = 0.f;
  invHitDensityYRatio = 0.f;
  invHitDensityZRatio = 0.f;

  std::map<float,Eigen::Vector3f>::const_iterator pcaIt;

  float maxY = 0.f;
  float maxZ = 0.f;
  float minY = 0.f;
  float minZ = 0.f;

  int iSegment = 1;

  float areaY = 0.f;
  float areaZ = 0.f;

  float areaYStart = 0.f;
  float areaYEnd = 0.f;
  float areaZStart = 0.f;
  float areaZEnd = 0.f;

  std::vector<float> pcaY;
  pcaY.clear();
  std::vector<float> pcaZ;
  pcaZ.clear();
  std::vector<float>::const_iterator pcaSegIt;

  const int nHits = pcaPositions.size();
  int nHitsStart = 0;
  int nHitsEnd = 0;

  for(pcaIt = pcaPositions.begin(); pcaIt != pcaPositions.end(); ++pcaIt)
    {
      if((pcaIt->first - paStart) / pfoLengthOnPA <= m_pfoFraction)
	nHitsStart++;
      if((pcaIt->first - paStart) / pfoLengthOnPA >= 1.0 - m_pfoFraction)
        nHitsEnd++; 

      if(pcaIt->first - paStart >= iSegment * m_segmentWidth)
	{
	  if(pcaY.size() >= 1)
	    {
	      sort(pcaY.begin(), pcaY.end());
	      pcaSegIt = pcaY.begin();
	      if((*pcaSegIt) < 0)
		minY = (*pcaSegIt);
	      sort(pcaY.begin(), pcaY.end(),std::greater<float>());
	      pcaSegIt = pcaY.begin();
	      if((*pcaSegIt) > 0)
		maxY = (*pcaSegIt);
	    }
	  if(pcaZ.size() >= 1)
	    {
	      sort(pcaZ.begin(), pcaZ.end());
	      pcaSegIt = pcaZ.begin();
	      if((*pcaSegIt) < 0)
		minZ = (*pcaSegIt);
	      sort(pcaZ.begin(), pcaZ.end(),std::greater<float>());
	      pcaSegIt = pcaZ.begin();
	      if((*pcaSegIt) > 0)
		maxZ = (*pcaSegIt);
	    }

	  areaY += fabs((maxY - minY) * m_segmentWidth);
	  areaZ += fabs((maxZ - minZ) * m_segmentWidth);
	  if((pcaIt->first - paStart) / pfoLengthOnPA <= m_pfoFraction)
	    {
	      areaYStart += fabs((maxY - minY) * m_segmentWidth);
	      areaZStart += fabs((maxZ - minZ) * m_segmentWidth);
	    }
	  if((pcaIt->first - paStart) / pfoLengthOnPA >= 1.0 - m_pfoFraction)
            {
              areaYEnd += fabs((maxY - minY) * m_segmentWidth);
              areaZEnd += fabs((maxZ - minZ) * m_segmentWidth);
            }

	  maxY = 0.f;
	  maxZ = 0.f;
	  minY = 0.f;
	  minZ = 0.f;
	  pcaY.clear();
	  pcaZ.clear();
	  iSegment = 1 + (int)((pcaIt->first - paStart) / m_segmentWidth);
	}

      pcaY.push_back(pcaIt->second.y());
      pcaZ.push_back(pcaIt->second.z()); 
    }

  if(nHits >= m_minHits)
    {
      invHitDensityY = areaY / nHits;
      invHitDensityZ = areaZ / nHits;
    }
  if(nHitsStart >= m_minHits && nHitsEnd >= m_minHits)
    {
      if(areaYStart / nHitsStart > std::numeric_limits<float>::epsilon())
        invHitDensityYRatio = (areaYEnd / nHitsEnd) / (areaYStart / nHitsStart);
      if(areaZStart / nHitsStart > std::numeric_limits<float>::epsilon())
        invHitDensityZRatio = (areaZEnd / nHitsEnd) / (areaZStart / nHitsStart);
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------
  void ThreeDPCAVariablesFeatureTool::SegmentsWithDoubleHits(const float paStart, std::map<float,Eigen::Vector3f> &pcaPositions, 
                                                             int &nSegmentsDoubleHits, float &nSegmentsDoubleHitsRatio)
{

  std::map<float,Eigen::Vector3f>::const_iterator pcaIt;

  std::vector<float> pcaYPlus;
  pcaYPlus.clear();
  std::vector<float> pcaZPlus;
  pcaZPlus.clear();
  std::vector<float> pcaYMinus;
  pcaYMinus.clear();
  std::vector<float> pcaZMinus;
  pcaZMinus.clear();
  std::vector<float>::const_iterator pcaSegIt;

  nSegmentsDoubleHits = 0;
  bool segmentDoubleHitsFound = false;

  int iSegment = 1;
  float prevDist = 0.f;

  for(pcaIt = pcaPositions.begin(); pcaIt != pcaPositions.end(); ++pcaIt)
    {

      if(pcaIt->first - paStart >= iSegment * m_segmentWidth)
	{
	  if(pcaYPlus.size() >= 2)
	    {
	      sort(pcaYPlus.begin(), pcaYPlus.end());
	      prevDist = 0.f;

	      for(pcaSegIt = pcaYPlus.begin(); pcaSegIt != pcaYPlus.end(); ++pcaSegIt)
		{
		  if(pcaSegIt != pcaYPlus.begin() && (*pcaSegIt) - prevDist > m_MoliereFraction * m_MoliereRadius && !segmentDoubleHitsFound)
		    {
		      nSegmentsDoubleHits++;
		      segmentDoubleHitsFound = true;
		    }
		  prevDist = (*pcaSegIt);
		}
	    }
	  if(pcaZPlus.size() >= 2)
	    {
	      sort(pcaZPlus.begin(), pcaZPlus.end());
	      prevDist = 0.f;

	      for(pcaSegIt = pcaZPlus.begin(); pcaSegIt != pcaZPlus.end(); ++pcaSegIt)
		{
		  if(pcaSegIt != pcaZPlus.begin() && (*pcaSegIt) - prevDist > m_MoliereFraction * m_MoliereRadius && !segmentDoubleHitsFound)
		    {
		      nSegmentsDoubleHits++;
		      segmentDoubleHitsFound = true;
		    }
		  prevDist = (*pcaSegIt);
		}
	    }
	  if(pcaYMinus.size() >= 2)
	    {
	      sort(pcaYMinus.begin(), pcaYMinus.end(), std::greater<float>());
	      prevDist = 0.f;

	      for(pcaSegIt = pcaYMinus.begin(); pcaSegIt != pcaYMinus.end(); ++pcaSegIt)
		{
		  if(pcaSegIt != pcaYMinus.begin() && fabs((*pcaSegIt) - prevDist) > m_MoliereFraction * m_MoliereRadius && !segmentDoubleHitsFound)
		    {
		      nSegmentsDoubleHits++;
		      segmentDoubleHitsFound = true;
		    }
		  prevDist = (*pcaSegIt);
		}
	    }
	  if(pcaZMinus.size() >= 2)
	    {
	      sort(pcaZMinus.begin(), pcaZMinus.end(), std::greater<float>());
	      prevDist = 0.f;
	      for(pcaSegIt = pcaZMinus.begin(); pcaSegIt != pcaZMinus.end(); ++pcaSegIt)
		{
		  if(pcaSegIt != pcaZMinus.begin() && fabs((*pcaSegIt) - prevDist) > m_MoliereFraction * m_MoliereRadius && !segmentDoubleHitsFound)
		    {
		      nSegmentsDoubleHits++;
		      segmentDoubleHitsFound = true;
		    }
		  prevDist = (*pcaSegIt);
		}
	    }

	  pcaYPlus.clear();
	  pcaZPlus.clear();
	  pcaYMinus.clear();
	  pcaZMinus.clear();
	  segmentDoubleHitsFound = false;
	  iSegment = 1 + (int)((pcaIt->first - paStart) / m_segmentWidth);
	}

      if(pcaIt->second.y() > 0)
	pcaYPlus.push_back(pcaIt->second.y());
      if(pcaIt->second.z() > 0)
	pcaZPlus.push_back(pcaIt->second.z());
      if(pcaIt->second.y() < 0)
	pcaYMinus.push_back(pcaIt->second.y());
      if(pcaIt->second.z() < 0)
	pcaZMinus.push_back(pcaIt->second.z());

    }

  if(iSegment >= m_minHits)
    {
      nSegmentsDoubleHitsRatio = (float)nSegmentsDoubleHits / iSegment;
    } 	
}

//------------------------------------------------------------------------------------------------------------------------------------------ 
StatusCode ThreeDPCAVariablesFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                                                                                                       "PfoFraction", m_pfoFraction));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                                                                                                       "MinHits", m_minHits));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
												       "SegmentWidth", m_segmentWidth));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
												       "MoliereFraction", m_MoliereFraction));

  return STATUS_CODE_SUCCESS;
}

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
    LArMvaHelper::MvaFeature charge1, charge2;

    ClusterList pClusterList;
    LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, pClusterList);
    if ((!pClusterList.empty()) && (pClusterList.size() == 1))
    {
        const Cluster *const pCluster(pClusterList.front());
        this->CalculateChargeVariables(pAlgorithm, pCluster, totalCharge, chargeSigma, chargeMean, endCharge);
    }

    if (chargeMean > std::numeric_limits<float>::epsilon())
        charge1 = chargeSigma / chargeMean;

    if (totalCharge > std::numeric_limits<float>::epsilon())
        charge2 = endCharge / totalCharge;

    featureVector.push_back(charge1);
    featureVector.push_back(charge2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::CalculateChargeVariables(const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, float &totalCharge,
    float &chargeSigma, float &chargeMean, float &endCharge)
{

    CaloHitList orderedCaloHitList;
    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, pCluster, orderedCaloHitList);

    const int totalHits(pCluster->GetNCaloHits());
    FloatVector chargeVector;
    int hitCounter(0);
    totalCharge = 0.f;
    endCharge = 0.f;

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

            if (hitCounter >= std::floor(totalHits*(1.f-m_endChargeFraction)))
            {
                endCharge += pCaloHitCharge;
            }
        }
    }

    if (!chargeVector.empty())
    {
        chargeMean = 0.f;
        chargeSigma = 0.f;

        for (const float charge : chargeVector)
            chargeMean += charge;

        chargeMean /= static_cast<float>(chargeVector.size());

        for (const float charge : chargeVector)
            chargeSigma += (charge - chargeMean) * (charge - chargeMean);

        chargeSigma = std::sqrt(chargeSigma / static_cast<float>(chargeVector.size()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::OrderCaloHitsByDistanceToVertex(const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, CaloHitList &caloHitList)
{
    //find the neutrino vertex and sort hits by distance to vertex
    const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if ((!pVertexList->empty()) && (pVertexList->size() == 1) && (VERTEX_3D == pVertexList->front()->GetVertexType()))
    {
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
