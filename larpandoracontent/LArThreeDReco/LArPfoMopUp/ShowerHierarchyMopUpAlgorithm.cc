/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/ShowerHierarchyMopUpAlgorithm.cc
 *
 *  @brief  Implementation of the shower hierarchy mop up algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/ShowerHierarchyMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ShowerHierarchyMopUpAlgorithm::Run()
{
    const PfoList *pLeadingPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_leadingPfoListName, pLeadingPfoList));

    if (!pLeadingPfoList || pLeadingPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerHierarchyMopUpAlgorithm: unable to find pfos in provided list, " << m_leadingPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    PfoList parentShowerPfos;
    this->FindParentShowerPfos(pLeadingPfoList, parentShowerPfos);
if (!parentShowerPfos.empty())
{
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1., -1., 1.);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &parentShowerPfos, "BeforeMerge", RED, true, false);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), pLeadingPfoList, "Leading", RED, true, true);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
}
    this->PerformPfoMerges(parentShowerPfos);
if (!parentShowerPfos.empty())
{
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1., -1., 1.);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &parentShowerPfos, "AfterMerge", GREEN, true, false);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), pLeadingPfoList, "Leading", GREEN, true, true);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
}
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHierarchyMopUpAlgorithm::FindParentShowerPfos(const PfoList *const pLeadingPfoList, PfoList &parentShowerPfos) const
{
    for (const Pfo *const pLeadingPfo : *pLeadingPfoList)
    {
        this->FindParentShowerPfos(pLeadingPfo, parentShowerPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHierarchyMopUpAlgorithm::FindParentShowerPfos(const Pfo *const pPfo, PfoList &parentShowerPfos) const
{
    if (LArPfoHelper::IsShower(pPfo))
    {
        if (pPfo->GetDaughterPfoList().empty())
            return;

        if (parentShowerPfos.end() != std::find(parentShowerPfos.begin(), parentShowerPfos.end(), pPfo))
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        
        parentShowerPfos.push_back(pPfo);
    }
    else
    {
        for (const Pfo *const pDaughterPfo : pPfo->GetDaughterPfoList())
            this->FindParentShowerPfos(pDaughterPfo, parentShowerPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerHierarchyMopUpAlgorithm::PerformPfoMerges(const PfoList &parentShowerPfos) const
{
    for (const Pfo *const pParentShowerPfo : parentShowerPfos)
    {
        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pParentShowerPfo, downstreamPfos);

        for (const Pfo *const pDownstreamPfo : downstreamPfos)
        {
            if (pDownstreamPfo != pParentShowerPfo)
                this->MergeAndDeletePfos(pParentShowerPfo, pDownstreamPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerHierarchyMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "LeadingPfoListName", m_leadingPfoListName));

    return PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
