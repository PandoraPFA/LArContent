/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/UnattachedDeltaRaysAlgorithm.cc
 *
 *  @brief  Implementation of the unattached delta rays algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/UnattachedDeltaRaysAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode UnattachedDeltaRaysAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "UnattachedDeltaRaysAlgorithm: pfo list " << m_pfoListName << " unavailable." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    PfoList unattachedDeltaRays;

    for (const Pfo *const pPfo : *pPfoList)
    {
        if ((0 == pPfo->GetNParentPfos()) && LArPfoHelper::IsShower(pPfo))
            unattachedDeltaRays.push_back(pPfo);
    }

    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(unattachedDeltaRays, allPfosToDelete);

    for (const Pfo *const pPfoToDelete : allPfosToDelete)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, m_pfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode UnattachedDeltaRaysAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
