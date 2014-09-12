/**
 *  @file   LArContent/src/LArCheating/CosmicRayShowerMatchingAlg.cc
 * 
 *  @brief  Implementation of the cheater for the cosmic ray shower matching algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArCheating/CheatingCosmicRayShowerMatchingAlg.h"

using namespace pandora;

namespace lar_content
{

StatusCode CheatingCosmicRayShowerMatchingAlg::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    PfoVector pfoVector(pPfoList->begin(), pPfoList->end());
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);

    for (PfoVector::const_iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;
        const ClusterList &pfoClusterList(pPfo->GetClusterList());

        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pPfoCluster = *cIter;
            const HitType hitType(LArClusterHelper::GetClusterHitType(pPfoCluster));

            if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            {
                if (TPC_3D == hitType)
                    continue;

                std::cout << "CheatingCosmicRayShowerMatchingAlg: Encountered unexpected hit type " << std::endl;
                return STATUS_CODE_INVALID_PARAMETER;
            }

            const StringVector &clusterListNames((TPC_VIEW_U == hitType) ? m_inputClusterListNamesU :
                (TPC_VIEW_V == hitType) ? m_inputClusterListNamesV : m_inputClusterListNamesW);

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CosmicRayShowerMatching(clusterListNames, pPfoCluster, pPfo));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayShowerMatchingAlg::CosmicRayShowerMatching(const StringVector &clusterListNames, const Cluster *const pPfoCluster,
    ParticleFlowObject *pPfo) const
{
    try
    {
        const MCParticle *pPfoMCParticle(MCParticleHelper::GetMainMCParticle(pPfoCluster));
        const MCParticle *pPfoParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pPfoMCParticle));

        for (StringVector::const_iterator sIter = clusterListNames.begin(), sIterEnd = clusterListNames.end(); sIter != sIterEnd; ++sIter)
        {
            const ClusterList *pClusterList = NULL;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *sIter, pClusterList));

            for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
            {
                Cluster *pCluster = *cIter;

                if (!pCluster->IsAvailable() || (pPfoCluster == pCluster))
                    continue;

                try
                {
                    const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
                    const MCParticle *pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

                    if (!LArMCParticleHelper::IsNeutrino(pParentMCParticle) && (pPfoParentMCParticle == pParentMCParticle))
                    {
                        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster));
                    }
                }
                catch (StatusCodeException &)
                {
                }
            }
        }
    }
    catch (StatusCodeException &)
    {
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayShowerMatchingAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputPfoListName", m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesU", m_inputClusterListNamesU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesV", m_inputClusterListNamesV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesW", m_inputClusterListNamesW));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
