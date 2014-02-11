/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding//LengthSeedFindingAlgorithm.cc
 * 
 *  @brief  Implementation of the length seed finding algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArSeedFinding/LengthSeedFindingAlgorithm.h"

using namespace pandora;

namespace lar
{

void LengthSeedFindingAlgorithm::GetSeedClusterList(const ClusterVector &candidateClusters, ClusterList &seedClusterList) const
{
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    for (ClusterVector::const_iterator iter = candidateClusters.begin(), iterEnd = candidateClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;
        const float clusterLength(LArClusterHelper::GetLength(pCluster));

        if (clusterLength > m_lengthCut)
            seedClusterList.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LengthSeedFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_lengthCut = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "LengthCut", m_lengthCut));

    return SeedFindingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
