/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray track consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRayTrackConsolidationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));


    ClusterVector shortClusters, longClusters;

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minTrackLength * m_minTrackLength)
        {
            shortClusters.push_back(pCluster);
        }
        else
        {
            longClusters.push_back(pCluster);
        }
    }

    std::sort(longClusters.begin(), longClusters.end(), LArClusterHelper::SortByNHits);

 

    const unsigned int m_halfWindowLayers(20);

   
    
    TwoDSlidingFitResultMap slidingFitResultMap;

    for (ClusterVector::const_iterator iter = longClusters.begin(), iterEnd = longClusters.end(); iter != iterEnd; ++iter)
    {
         TwoDSlidingFitResult slidingFitResult;
         LArClusterHelper::LArTwoDSlidingFit(*iter, m_halfWindowLayers, slidingFitResult);

         if (!slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
              throw StatusCodeException(STATUS_CODE_FAILURE);
    }     
    


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minTrackLength = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLength", m_minTrackLength));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
