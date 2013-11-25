/**
 *  @file   LArContent/src/LArTwoDSeed/LengthSeedFindingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex seed finding algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDSeed/LengthSeedFindingAlgorithm.h"

using namespace pandora;

namespace lar
{

void LengthSeedFindingAlgorithm::GetSeedClusterList(const ClusterVector &candidateClusters, ClusterList &seedClusterList) const
{
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    unsigned int nIter((NULL != pSeedClusterList) ? pSeedClusterList->size() : 0);

    for (ClusterVector::const_iterator iter = candidateClusters.begin(), iterEnd = candidateClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;
        unsigned int lengthCut((nIter < m_finalChangeIter) ? m_initialLengthCut : m_finalLengthCut);

        if ((nIter > m_initialChangeIter) && (nIter < m_finalChangeIter))
        {
            lengthCut = m_initialLengthCut + (nIter - m_initialChangeIter) * (m_finalLengthCut - m_initialLengthCut) / (m_finalChangeIter - m_initialChangeIter);
        }

        if (pCluster->GetOrderedCaloHitList().size() >= lengthCut)
        {
            seedClusterList.insert(pCluster);
        }

        ++nIter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LengthSeedFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_initialLengthCut = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "InitialLengthCut", m_initialLengthCut));

    m_finalLengthCut = 50;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "FinalLengthCut", m_finalLengthCut));

    m_initialChangeIter = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "InitialChangeIter", m_initialChangeIter));

    m_finalChangeIter = 6;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "FinalChangeIter", m_finalChangeIter));

    return SeedFindingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
