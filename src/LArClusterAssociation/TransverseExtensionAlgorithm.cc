/**
 *  @file   LArContent/src/LArClusterAssociation/TransverseExtensionAlgorithm.cc
 * 
 *  @brief  Implementation of the transverse extension algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/TransverseExtensionAlgorithm.h"

#include "LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar
{

StatusCode TransverseExtensionAlgorithm::Run()
{


    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));



    typedef std::map<pandora::Cluster*,LArClusterHelper::TwoDSlidingFitResult> TwoDSlidingFitResultMap;
    TwoDSlidingFitResultMap slidingFitResultMap;


    bool carryOn(true);

    while(carryOn)
    {
        carryOn = false;

        ClusterVector longVector, shortVector;
        this->SortInputClusters(pClusterList, longVector, shortVector);

        std::sort(longVector.begin(), longVector.end(), LArClusterHelper::SortByNHits);
        std::sort(shortVector.begin(), shortVector.end(), LArClusterHelper::SortByNHits);

        for (ClusterVector::const_iterator iter = longVector.begin(), iterEnd = longVector.end(); iter != iterEnd; ++iter)
        {
            if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
	    {             
                LArClusterHelper::TwoDSlidingFitResult slidingFitResult;
                LArClusterHelper::LArTwoDSlidingFit(*iter, 20, slidingFitResult);

                if (!slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
	    }
	}

        for (ClusterVector::iterator iterI = longVector.begin(), iterEndI = longVector.end(); iterI != iterEndI; ++iterI)
        {
            TwoDSlidingFitResultMap::const_iterator iterFit = slidingFitResultMap.find(*iterI);
            if (slidingFitResultMap.end() == iterFit)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult(iterFit->second);

            bool isExpanded(false);

            for (ClusterVector::iterator iterJ = shortVector.begin(), iterEndJ = shortVector.end(); iterJ != iterEndJ; ++iterJ)
            {
                if (NULL == *iterJ)
	            continue;

                if (this->IsAssociated(slidingFitResult, *iterJ))
	        {
// ClusterList tempListI, tempListJ;
// tempListI.insert(*iterI);
// tempListJ.insert(*iterJ);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempListI, "ClusterToEnlarge", RED));
// PANDORA_MONITORING_API(VisualizeClusters(&tempListJ, "ClusterToDelete", BLUE));
// PANDORA_MONITORING_API(ViewEvent());

                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, *iterI, *iterJ));
                    isExpanded = true;
                    *iterJ = NULL;
	        }

	    }

            if (isExpanded)
	    {
	        slidingFitResultMap.erase(*iterI);
                *iterI = NULL;
                carryOn = true;
	    }
	}
    }



    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::SortInputClusters(const ClusterList *const pClusterList, ClusterVector &longVector, ClusterVector &shortVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) > m_minClusterLength * m_minClusterLength)
	{
	    longVector.push_back(pCluster);
	}
        else
	{
	    shortVector.push_back(pCluster);
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseExtensionAlgorithm::IsAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const Cluster *const pCluster) const
{
    return (this->IsEndAssociated(slidingFitResult,pCluster) || this->IsMidAssociated(slidingFitResult,pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseExtensionAlgorithm::IsEndAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const Cluster *const pCluster) const
{

    CartesianVector firstCoordinate(0.f, 0.f, 0.f);
    CartesianVector secondCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster, firstCoordinate, secondCoordinate);

    for (unsigned int iForward = 0; iForward < 2; ++iForward)
    {
        const CartesianVector vertex(1==iForward ? slidingFitResult.GetGlobalMinLayerPosition() : slidingFitResult.GetGlobalMaxLayerPosition());
        const CartesianVector direction(1==iForward ? slidingFitResult.GetGlobalMinLayerDirection() : slidingFitResult.GetGlobalMaxLayerDirection() * -1.f);

        float firstL(0.f), firstT(0.f), secondT(0.f), secondL(0.f);
        LArVertexHelper::GetImpactParameters(vertex, direction, firstCoordinate, firstL, firstT);
        LArVertexHelper::GetImpactParameters(vertex, direction, secondCoordinate, secondL, secondT);

        const float innerL(firstL < secondL ? firstL : secondL);
        const float innerT(firstL < secondL ? firstT : secondT);
        const float outerT(firstL > secondL ? firstT : secondT);

        if ( innerL > 0.f && innerL < m_maxLongitudinalDisplacement && 
             innerT < m_maxTransverseDisplacement && outerT < 1.5f * m_maxTransverseDisplacement )
	    return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseExtensionAlgorithm::IsMidAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const Cluster *const pCluster) const
{
    bool isFirstAssociated(false), isSecondAssociated(false);
    float firstL(0.f), firstT(0.f), secondT(0.f), secondL(0.f); 

    CartesianVector firstCoordinate(0.f, 0.f, 0.f);
    CartesianVector secondCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pCluster, firstCoordinate, secondCoordinate);

    try{
        CartesianVector projectedFirstCoordinate(0.f, 0.f, 0.f);
        slidingFitResult.GetGlobalFitProjection(firstCoordinate, projectedFirstCoordinate);
        slidingFitResult.GetLocalPosition(projectedFirstCoordinate, firstL, firstT);
        isFirstAssociated = ((projectedFirstCoordinate - firstCoordinate).GetMagnitudeSquared() < m_maxTransverseDisplacement * m_maxTransverseDisplacement);
    }
    catch (StatusCodeException &)
    {
    }

    try{
        CartesianVector projectedSecondCoordinate(0.f, 0.f, 0.f);
        slidingFitResult.GetGlobalFitProjection(secondCoordinate, projectedSecondCoordinate);
        slidingFitResult.GetLocalPosition(projectedSecondCoordinate, secondL, secondT);
        isSecondAssociated = ((projectedSecondCoordinate - secondCoordinate).GetMagnitudeSquared() < m_maxTransverseDisplacement * m_maxTransverseDisplacement);
    }
    catch (StatusCodeException &)
    {
    }

    if (!isFirstAssociated || !isSecondAssociated)
        return false;


    const float minL(std::min(firstL,secondL));
    const float maxL(std::max(firstL,secondL));
    const int minLayer(slidingFitResult.GetLayer(minL));
    const int maxLayer(slidingFitResult.GetLayer(maxL));

    const TwoDSlidingFitResult::LayerFitContributionMap &layerFitContributionMap = slidingFitResult.GetLayerFitContributionMap();

    for (int iLayer = minLayer; iLayer<=maxLayer; ++iLayer)
    {
        if (layerFitContributionMap.end() != layerFitContributionMap.find(iLayer))
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterLength = 5.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    m_maxTransverseDisplacement = 0.75f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_maxLongitudinalDisplacement = 3.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
