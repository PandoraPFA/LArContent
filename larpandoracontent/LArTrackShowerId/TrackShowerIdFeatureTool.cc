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

ShowerFitFeatureTool::ShowerFitFeatureTool() :
    m_ratioVariables(true),
	m_slidingShowerFitWindow(3),
	m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerFitFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
	
	if (m_ratioVariables)
	{
		float ratio(-1.f);
		try
		{
			const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
			float straightLineLength = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
			if (straightLineLength > std::numeric_limits<double>::epsilon())
				ratio = (CutClusterCharacterisationAlgorithm::GetShowerFitWidth(pAlgorithm, pCluster, m_slidingShowerFitWindow))/straightLineLength;
		}
		catch (const StatusCodeException &)
		{
			ratio = -1.f;
		}
		featureVector.push_back(ratio);
	}
	else
	{
		featureVector.push_back(CutClusterCharacterisationAlgorithm::GetShowerFitWidth(pAlgorithm, pCluster, m_slidingShowerFitWindow));
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingShowerFitWindow", m_slidingShowerFitWindow));
	
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RatioVariables", m_ratioVariables));
	
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));
		
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
    m_slidingLinearFitWindowLarge(10000),
	m_ratioVariables(true),
	m_addStraightLineLength(true),
    m_addDiffWithStraightLineMean(true),
    m_addDiffWithStraightLineSigma(false),
    m_addDTDLWidth(true),
    m_addMaxFitGapLength(true),
    m_addRMSLinearFit(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LinearFitFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
const pandora::Cluster * const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float dTdLWidth(-1.f), straightLineLengthLarge(-1.f), diffWithStraightLineMean(-1.f), diffWithStraightLineSigma(-1.f), maxFitGapLength(-1.f), rmsSlidingLinearFit(-1.f);
    this->CalculateVariablesSlidingLinearFit(pCluster, straightLineLengthLarge, diffWithStraightLineMean, diffWithStraightLineSigma, dTdLWidth, maxFitGapLength, rmsSlidingLinearFit);

    if (m_ratioVariables && (straightLineLengthLarge > std::numeric_limits<double>::epsilon()))
	{
		diffWithStraightLineMean  /= straightLineLengthLarge;
		diffWithStraightLineSigma /= straightLineLengthLarge;
		dTdLWidth                 /= straightLineLengthLarge;
		maxFitGapLength           /= straightLineLengthLarge;
		rmsSlidingLinearFit       /= straightLineLengthLarge;
	}

	if(m_addStraightLineLength)
		featureVector.push_back(straightLineLengthLarge);

    if (m_addDiffWithStraightLineMean)
        featureVector.push_back(diffWithStraightLineMean);

    if (m_addDiffWithStraightLineSigma)
        featureVector.push_back(diffWithStraightLineSigma);

    if (m_addDTDLWidth)
        featureVector.push_back(dTdLWidth);

    if (m_addMaxFitGapLength)
        featureVector.push_back(maxFitGapLength);

    if (m_addRMSLinearFit)
        featureVector.push_back(rmsSlidingLinearFit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge,
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
            rmsSlidingLinearFit += mapEntry.second.GetRms();

            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);

            LayerFitResultMap::const_iterator iterLarge = slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(mapEntry.second.GetL()));

            if (slidingFitResultLarge.GetLayerFitResultMap().end() == iterLarge)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            diffWithStraightLineVector.push_back(static_cast<float>(std::fabs(mapEntry.second.GetFitT() - iterLarge->second.GetFitT())));

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

StatusCode LinearFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));
	
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddStraightLineLength", m_addStraightLineLength));	

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddDiffWithStraightLineMean", m_addDiffWithStraightLineMean));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddDiffWithStraightLineSigma", m_addDiffWithStraightLineSigma));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddDTDLWidth", m_addDTDLWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddMaxFitGapLength", m_addMaxFitGapLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddRMSLinearFit", m_addRMSLinearFit));
	
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RatioVariables", m_ratioVariables));	

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

VertexDistanceFeatureTool::VertexDistanceFeatureTool() :
    m_ratioVariables(true),
	m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexDistanceFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
	
	if (m_ratioVariables)
	{
		float straightLineLength(-1.f), ratio(-1.f);
		try
		{
			const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
			straightLineLength = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
			if (straightLineLength > std::numeric_limits<double>::epsilon())
				ratio = (CutClusterCharacterisationAlgorithm::GetVertexDistance(pAlgorithm, pCluster))/straightLineLength;
		}
		catch (const StatusCodeException &)
		{
			ratio = -1.f;
		}
		featureVector.push_back(ratio);
	}
	else
	{
		featureVector.push_back((CutClusterCharacterisationAlgorithm::GetVertexDistance(pAlgorithm, pCluster)));
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexDistanceFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RatioVariables", m_ratioVariables));
		
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));
		
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

    ThreeDLinearFitFeatureTool::ThreeDLinearFitFeatureTool() :
    m_slidingLinearFitWindow(3),
    m_slidingLinearFitWindowLarge(10000),
	m_ratioVariables(true),
	m_addStraightLineLength(true),
    m_addDiffWithStraightLineMean(true),
    m_addDiffWithStraightLineSigma(false),
    m_addDTDLWidth(true),
    m_addMaxFitGapLength(true),
    m_addRMSLinearFit(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

	ClusterList clusterList;                                                                                                                                    
	LArPfoHelper::GetTwoDClusterList(pInputPfo, clusterList); 
	float dTdLWidth(0.f),  diffWithStraightLineMean(0.f), diffWithStraightLineSigma(0.f), maxFitGapLength(0.f), rmsSlidingLinearFit(0.f);  
	
	for (const Cluster *const pCluster : clusterList)
	{
		float dTdLWidthCluster(-1.f), straightLineLengthLargeCluster(-1.f), diffWithStraightLineMeanCluster(-1.f), diffWithStraightLineSigmaCluster(-1.f), 
			maxFitGapLengthCluster(-1.f), rmsSlidingLinearFitCluster(-1.f);
		this->CalculateVariablesSlidingLinearFit(pCluster, straightLineLengthLargeCluster, diffWithStraightLineMeanCluster, diffWithStraightLineSigmaCluster, 
			dTdLWidthCluster, maxFitGapLengthCluster, rmsSlidingLinearFitCluster);
		
		if (m_ratioVariables && (straightLineLengthLargeCluster > std::numeric_limits<double>::epsilon()))
		{
			diffWithStraightLineMeanCluster  /= straightLineLengthLargeCluster;
			diffWithStraightLineSigmaCluster /= straightLineLengthLargeCluster;
			dTdLWidthCluster                 /= straightLineLengthLargeCluster;
			maxFitGapLengthCluster           /= straightLineLengthLargeCluster;
			rmsSlidingLinearFitCluster       /= straightLineLengthLargeCluster;
		}
	
			diffWithStraightLineMean   += diffWithStraightLineMeanCluster;
			diffWithStraightLineSigma  += diffWithStraightLineSigmaCluster;
			dTdLWidth                  += dTdLWidthCluster;
			maxFitGapLength            += maxFitGapLengthCluster;
			rmsSlidingLinearFit        += rmsSlidingLinearFitCluster;		
	}
	
	const int nClusters(clusterList.size());
	diffWithStraightLineMean   /= nClusters;
	diffWithStraightLineSigma  /= nClusters;
	dTdLWidth                  /= nClusters;
	maxFitGapLength            /= nClusters;
	rmsSlidingLinearFit        /= nClusters;	

	if(m_addStraightLineLength)
	{
		ClusterList clusterList3D;                                                                                                                                    
		LArPfoHelper::GetThreeDClusterList(pInputPfo, clusterList3D);                                                                                                      
                                                                                                                                                                        
		if (1 != clusterList3D.size())                                                                                                                                
			throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);     
                                                                                                                                                                        
		const Cluster *const pCluster(clusterList3D.front());       
		const ThreeDSlidingFitResult sliding3DFitResultLarge(pCluster, m_slidingLinearFitWindowLarge, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
		const TwoDSlidingFitResult slidingFitResult3DLarge(sliding3DFitResultLarge.GetFirstFitResult());
		featureVector.push_back((slidingFitResult3DLarge.GetGlobalMaxLayerPosition() - slidingFitResult3DLarge.GetGlobalMinLayerPosition()).GetMagnitude());
	}
    if (m_addDiffWithStraightLineMean)
		featureVector.push_back(diffWithStraightLineMean);
	
    if (m_addDiffWithStraightLineSigma)
		featureVector.push_back(diffWithStraightLineSigma);
	
    if (m_addDTDLWidth)
		featureVector.push_back(dTdLWidth);

    if (m_addMaxFitGapLength)
		featureVector.push_back(maxFitGapLength);

    if (m_addRMSLinearFit)
		featureVector.push_back(rmsSlidingLinearFit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge,
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
            rmsSlidingLinearFit += mapEntry.second.GetRms();

            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);

            LayerFitResultMap::const_iterator iterLarge = slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(mapEntry.second.GetL()));

            if (slidingFitResultLarge.GetLayerFitResultMap().end() == iterLarge)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            diffWithStraightLineVector.push_back(static_cast<float>(std::fabs(mapEntry.second.GetFitT() - iterLarge->second.GetFitT())));

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

StatusCode ThreeDLinearFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));

	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddStraightLineLength", m_addStraightLineLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddDiffWithStraightLineMean", m_addDiffWithStraightLineMean));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddDiffWithStraightLineSigma", m_addDiffWithStraightLineSigma));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddDTDLWidth", m_addDTDLWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddMaxFitGapLength", m_addMaxFitGapLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddRMSLinearFit", m_addRMSLinearFit));
		
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RatioVariables", m_ratioVariables));	

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDVertexDistanceFeatureTool::ThreeDVertexDistanceFeatureTool() :
    m_ratioVariables(false),
	m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexDistanceFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
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
	const Vertex *const pVertex = LArPfoHelper::GetVertex(pInputPfo);
	const CartesianVector nuPosition(nuVertex->GetPosition()), pPosition(pVertex->GetPosition());

	float vertexDistance((nuPosition - pPosition).GetMagnitudeSquared());	 
	
	if (m_ratioVariables)
	{
		float ratio(-1.f);
		try
		{
			ClusterList clusterList;                                                                                                                                    
			LArPfoHelper::GetThreeDClusterList(pInputPfo, clusterList);                                                                                                      
                                                                                                                                                                        
			if (1 == clusterList.size())   
			{																																																							
				const Cluster *const pCluster(clusterList.front()); 
	
				const ThreeDSlidingFitResult sliding3DFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
				const TwoDSlidingFitResult slidingFitResult(sliding3DFitResult.GetFirstFitResult());
				float straightLineLength = (slidingFitResult.GetGlobalMaxLayerPosition() - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitude();
				if (straightLineLength > std::numeric_limits<double>::epsilon())
					ratio = vertexDistance / straightLineLength;
			}
		}
		catch (const StatusCodeException &)
		{
			ratio = -1.f;
		}
		featureVector.push_back(ratio);
	}
	else
	{
		featureVector.push_back(vertexDistance);
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDVertexDistanceFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RatioVariables", m_ratioVariables));
		
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingLinearFitWindow", m_slidingLinearFitWindow));
		
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

OpeningAngleFeatureTool::OpeningAngleFeatureTool() 
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OpeningAngleFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pInputPfo)
{
	if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
		
	// Need the 3D clusters and hits to calculate PCA components
	ClusterList threeDClusterList;
	LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);
	
	CaloHitList threeDCaloHitList; // Might be able to get through LArPfoHelper instead
	const Cluster *const pThreeDCluster(threeDClusterList.front());
	pThreeDCluster->GetOrderedCaloHitList().FillCaloHitList(threeDCaloHitList);
        
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
	
	float openingAngle(this->OpeningAngle(eigenVecsStart.at(0), eigenVecsStart.at(1), eigenValuesStart));
	float closingAngle(this->OpeningAngle(eigenVecsEnd.at(0), eigenVecsEnd.at(1), eigenValuesEnd));
	
	featureVector.push_back(std::abs(openingAngle-closingAngle));

}

//------------------------------------------------------------------------------------------------------------------------------------------

void OpeningAngleFeatureTool::Divide3DCaloHitList(const Algorithm *const pAlgorithm, CaloHitList threeDCaloHitList, CartesianPointVector &pointVectorStart, CartesianPointVector &pointVectorEnd)
{
	
	const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || (pVertexList->size() != 1) || (VERTEX_3D != pVertexList->front()->GetVertexType()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const Vertex *const pVertex(pVertexList->front());
	const CartesianVector nuVertex(pVertex->GetPosition().GetX(),pVertex->GetPosition().GetY(),pVertex->GetPosition().GetZ());
	CaloHitVector threeDCaloHitVector(threeDCaloHitList.begin(), threeDCaloHitList.end());
	
	//order by distance to vertex, so first ones are closer to nuvertex
	std::sort(threeDCaloHitVector.begin(), threeDCaloHitVector.end(), std::bind(ChargeFeatureTool::SortByDistanceToVertex, std::placeholders::_1, std::placeholders::_2, nuVertex));
	CaloHitList orderedCaloHitList(threeDCaloHitVector.begin(),threeDCaloHitVector.end());
	const int nhits(orderedCaloHitList.size());

	for (const CaloHit *const pCaloHit : orderedCaloHitList)
	{	
		if (pointVectorStart.size()<(nhits/2))//half in start, half in end
		{
			pointVectorStart.push_back(pCaloHit->GetPositionVector());	
		}
		else
		{
			pointVectorEnd.push_back(pCaloHit->GetPositionVector());
		}
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------
//TODO - use from PCAShowerBuilding, maybe move to PCAHelper?

float OpeningAngleFeatureTool::OpeningAngle(const CartesianVector &principal, const CartesianVector &secondary,
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

StatusCode OpeningAngleFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
		
    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PCAFeatureTool::PCAFeatureTool() :
    m_ratioVariables(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PCAFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pInputPfo)
{
	
	 if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
		
        // Need the 3D clusters and hits to calculate PCA components
        ClusterList threeDClusterList;
        LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);

        if (threeDClusterList.empty())
            return;

        CaloHitList threeDCaloHitList; // Might be able to get through LArPfoHelper instead
        const Cluster *const pThreeDCluster(threeDClusterList.front());
        pThreeDCluster->GetOrderedCaloHitList().FillCaloHitList(threeDCaloHitList);

        if (threeDCaloHitList.empty())
            return;

        // Run the PCA analysis
        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecs;
        LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        LArPcaHelper::RunPca(threeDCaloHitList, centroid, eigenValues, eigenVecs);
		const float principalEigenvalue(eigenValues.GetX()), secondaryEigenvalue(eigenValues.GetY()), tertiaryEigenvalue(eigenValues.GetZ());
		
		if (m_ratioVariables)
		{
			float pca1(-1.f), pca2(-1.f);	
			if(principalEigenvalue > std::numeric_limits<float>::epsilon())	
			{
				pca1 = secondaryEigenvalue/principalEigenvalue;
				pca2 = tertiaryEigenvalue/principalEigenvalue;
			}
			featureVector.push_back(pca1);	
			featureVector.push_back(pca2);
		}
		else
		{
			featureVector.push_back(principalEigenvalue);
			featureVector.push_back(secondaryEigenvalue);
			featureVector.push_back(tertiaryEigenvalue);
		}
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PCAFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RatioVariables", m_ratioVariables));
		
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ChargeFeatureTool::ChargeFeatureTool() :
	m_ratioVariables(true),
    m_addTotalCharge(false),
	m_addChargeSigma(true),
	m_addChargeIncrements(false),
	m_addStartCharge(false),
	m_addEndCharge(true),
	m_endChargeFraction(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChargeFeatureTool::CalculateChargeVariables(const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, float &totalCharge, float &chargeSigma, 
	float &chargeMean, float &startCharge, float &endCharge, float &incrementCharge, float endChargeFraction) 
{

	CaloHitList orderedCaloHitList(ChargeFeatureTool::OrderCaloHitsByDistanceToVertex(pAlgorithm, pCluster));
	
	const int totalHits(pCluster->GetNCaloHits());	
	FloatVector chargeVector;
	int hitCounter(0);
	totalCharge = 0.f;
	startCharge = 0.f;
	endCharge = 0.f; 
	incrementCharge = 0.f;
	float previousCharge(0.f);
	
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
			totalCharge	+= pCaloHitCharge;
			chargeVector.push_back(pCaloHitCharge);
			
			if (hitCounter <= std::ceil(totalHits*endChargeFraction)) //ensure at least 1 hit is used
			{
				startCharge += pCaloHitCharge;
			}
			if (hitCounter >= std::floor(totalHits*(1.f-endChargeFraction)))
			{	
				endCharge += pCaloHitCharge;
			}		
			
			if(hitCounter>1)
				incrementCharge = (pCaloHitCharge-previousCharge)*(pCaloHitCharge-previousCharge);
				
			previousCharge = pCaloHitCharge;
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
	incrementCharge = std::sqrt(incrementCharge / static_cast<float>(chargeVector.size()));
		
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList ChargeFeatureTool::OrderCaloHitsByDistanceToVertex(const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster) 
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

	//TODO: might give problems if vertex in the middle of the cluster 
	std::sort(clusterCaloHitVector.begin(), clusterCaloHitVector.end(), std::bind(SortByDistanceToVertex, std::placeholders::_1, std::placeholders::_2, vertexPosition2D));
	CaloHitList orderedCaloHitList(clusterCaloHitVector.begin(),clusterCaloHitVector.end());
	
	return orderedCaloHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------               

bool ChargeFeatureTool::SortByDistanceToVertex(const CaloHit *const left, const CaloHit *const right, const CartesianVector &nuVertex2D)
{
    float distanceL((left->GetPositionVector()-nuVertex2D).GetMagnitude());
    float distanceR((right->GetPositionVector()-nuVertex2D).GetMagnitude());
    return distanceL < distanceR;                                                                                                                                  
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChargeFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::Cluster *const pCluster)
{
	if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
		
	float totalCharge(-1.f), chargeSigma(-1.f), chargeMean(-1.f), startCharge(-1.f), endCharge(-1.f), incrementCharge(-1.f);
    ChargeFeatureTool::CalculateChargeVariables(pAlgorithm, pCluster, totalCharge, chargeSigma, chargeMean, startCharge, endCharge, incrementCharge, m_endChargeFraction);
	
	if (m_ratioVariables) 
	{
		if (chargeMean > std::numeric_limits<double>::epsilon())
			chargeSigma /= chargeMean;	
	
		if (endCharge > std::numeric_limits<double>::epsilon())
			startCharge /= endCharge;
		
		if (totalCharge > std::numeric_limits<double>::epsilon())
			endCharge /= totalCharge;
			
		if (totalCharge > std::numeric_limits<double>::epsilon())
			incrementCharge /= totalCharge;
	}
	
    if (m_addTotalCharge)
        featureVector.push_back(totalCharge);

    if (m_addChargeSigma)
        featureVector.push_back(chargeSigma);
		
	if 	(m_addChargeIncrements)
		featureVector.push_back(incrementCharge);

    if (m_addStartCharge)
        featureVector.push_back(startCharge);

    if (m_addEndCharge)
        featureVector.push_back(endCharge);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ChargeFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
	
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddTotalCharge", m_addTotalCharge));
		
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddChargeSigma", m_addChargeSigma));
	
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddChargeIncrements", m_addChargeIncrements));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddStartCharge", m_addStartCharge));	
	
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddEndCharge", m_addEndCharge));	

	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EndChargeFraction", m_endChargeFraction));	
   
	return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
// TODO - make it inherit from ChargeFeatureTool and convert into static ordered by distance to vertex etc ?
ThreeDChargeFeatureTool::ThreeDChargeFeatureTool() :
    m_ratioVariables(true),
	m_addTotalCharge(false),
	m_addChargeSigma(true),
	m_addChargeIncrements(false),
	m_addStartCharge(false),
	m_addEndCharge(true),
	m_endChargeFraction(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const Algorithm *const pAlgorithm,
    const pandora::ParticleFlowObject *const pInputPfo)
{
	if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
		
	float totalCharge(-1.f), chargeSigma(-1.f), chargeMean(-1.f), startCharge(-1.f), endCharge(-1.f), incrementCharge(-1.f);
	ClusterList pClusterList;
    LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, pClusterList);
	if ((!pClusterList.empty()) && (pClusterList.size() == 1)) // TODO, think what to do if there is no W cluster (for the charge info) 
	{
		const Cluster *const pCluster(pClusterList.front());
		ChargeFeatureTool::CalculateChargeVariables(pAlgorithm, pCluster, totalCharge, chargeSigma, chargeMean, startCharge, endCharge, incrementCharge, m_endChargeFraction);
	}

	if (m_ratioVariables) 
	{
		if (chargeMean > std::numeric_limits<double>::epsilon())
			chargeSigma /= chargeMean;	
	
		if (endCharge > std::numeric_limits<double>::epsilon())
			startCharge /= endCharge;
		
		if (totalCharge > std::numeric_limits<double>::epsilon())
			endCharge /= totalCharge;
			
		if (totalCharge > std::numeric_limits<double>::epsilon())
			incrementCharge /= totalCharge;
	}

    if (m_addTotalCharge)
        featureVector.push_back(totalCharge);

    if (m_addChargeSigma)
        featureVector.push_back(chargeSigma);
		
	if 	(m_addChargeIncrements)
		featureVector.push_back(incrementCharge);

    if (m_addStartCharge)
        featureVector.push_back(startCharge);

    if (m_addEndCharge)
        featureVector.push_back(endCharge);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDChargeFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
	
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddTotalCharge", m_addTotalCharge));
		
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddChargeSigma", m_addChargeSigma));
	
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddChargeIncrements", m_addChargeIncrements));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddStartCharge", m_addStartCharge));	
	
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AddEndCharge", m_addEndCharge));		
   
	return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
