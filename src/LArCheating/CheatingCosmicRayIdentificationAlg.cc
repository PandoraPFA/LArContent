/**
 *  @file   LArContent/src/LArCheating/CheatingCosmicRayIdentificationAlg.cc
 * 
 *  @brief  Implementation of the cheater for the cosmic ray identification algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArMCParticleHelper.h"

#include "LArCheating/CheatingCosmicRayIdentificationAlg.h"

using namespace pandora;

namespace lar
{

StatusCode CheatingCosmicRayIdentificationAlg::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    PfoList outputPfoList, outputDaughterPfoList;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;

        bool isCosmicRay(false);
        const ClusterList &clusterList(pPfo->GetClusterList());

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;

            if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
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
            const PfoList &daughterPfoList(pPfo->GetDaughterPfoList());

            for (PfoList::const_iterator dIter = daughterPfoList.begin(), dIterEnd = daughterPfoList.end(); dIter != dIterEnd; ++dIter)
            {
                ParticleFlowObject *pDaughterPfo = *dIter;

                if (outputPfoList.count(pDaughterPfo) || !pDaughterPfo->GetDaughterPfoList().empty())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                outputDaughterPfoList.insert(pDaughterPfo);
            }
        }
    }

    if (!outputPfoList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_inputPfoListName, m_outputPfoListName, outputPfoList));

        if (!outputDaughterPfoList.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_inputDaughterPfoListName, m_outputDaughterPfoListName, outputDaughterPfoList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayIdentificationAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputDaughterPfoListName", m_inputDaughterPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputDaughterPfoListName", m_outputDaughterPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
