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

namespace lar_content
{

StatusCode CheatingCosmicRayIdentificationAlg::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    if (NULL == pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingCosmicRayIdentificationAlg: pfo list " << m_inputPfoListName << " unavailable." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

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

                if (!LArMCParticleHelper::GetParentNeutrinoId(pMCParticle))
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

    m_inputDaughterPfoListName = m_inputPfoListName;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputDaughterPfoListName", m_inputDaughterPfoListName));

    m_outputDaughterPfoListName = m_outputPfoListName;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputDaughterPfoListName", m_outputDaughterPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
