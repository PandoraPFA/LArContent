/**
 *  @file   LArContent/src/LArCheating/CheatingCosmicRayIdentificationAlg.cc
 * 
 *  @brief  Implementation of the cheater for the cosmic ray identification algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArCheating/CheatingCosmicRayIdentificationAlg.h"

using namespace pandora;

namespace lar
{

StatusCode CheatingCosmicRayIdentificationAlg::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    PfoList outputPfoList;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;

        bool isCosmicRay(false);
        const ClusterList &clusterList(pPfo->GetClusterList());

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;

            if (TPC_3D == LArThreeDHelper::GetClusterHitType(pCluster))
                continue;

            try
            {
                const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

                if (!LArMCParticleHelper::GetPrimaryNeutrino(pMCParticle))
                    isCosmicRay = true;
            }
            catch (StatusCodeException &)
            {
                isCosmicRay = false;
                break;
            }
        }

        if (isCosmicRay)
        {
            outputPfoList.insert(pPfo);
        }
    }

    if (!outputPfoList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_inputPfoListName, m_outputPfoListName, outputPfoList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayIdentificationAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
