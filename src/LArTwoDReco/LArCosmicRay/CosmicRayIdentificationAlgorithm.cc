/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray identification algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRayIdentificationAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    PfoList outputPfoList;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;

        bool isCosmicRay(true);
        const ClusterList &clusterList(pPfo->GetClusterList());

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;

            try
            {
                const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

                if (LArMCParticleHelper::GetPrimaryNeutrino(pMCParticle))
                    isCosmicRay = false;
            }
            catch (StatusCodeException &)
            {
                isCosmicRay = false;
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

StatusCode CosmicRayIdentificationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
