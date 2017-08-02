/**
 *  @file   larpandoracontent/LArUtility/ParentAlgorithm.cc
 *
 *  @brief  Implementation of the parent algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArUtility/ParentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ParentAlgorithm::Run()
{
    PfoList parentCosmicRayPfos;
    this->RunAllHitsCosmicRayReconstruction(parentCosmicRayPfos);

    if (parentCosmicRayPfos.empty())
        return STATUS_CODE_SUCCESS;

    PfoList ambiguousPfos;
    this->FindAmbiguousPfos(parentCosmicRayPfos, ambiguousPfos);
    this->RemoveAmbiguousCosmicRayPfos(ambiguousPfos);

    SliceList sliceList;
    this->RunFastReconstructionAndSlicing(sliceList);

    if (sliceList.empty())
        return STATUS_CODE_SUCCESS;

    SliceIndexToPfoListMap sliceToCosmicRayPfosMap;
    SliceIndexToPropertiesMap sliceIndexToPropertiesMap;
    this->ReconstructSlices(sliceList, sliceToCosmicRayPfosMap, sliceIndexToPropertiesMap);

    if (sliceToCosmicRayPfosMap.empty())
        return STATUS_CODE_SUCCESS;

    unsigned int neutrinoSliceIndex(0);
    if (this->GetNeutrinoSliceIndex(sliceIndexToPropertiesMap, neutrinoSliceIndex))
    {
        this->RemoveSliceCosmicRayReconstruction(sliceToCosmicRayPfosMap, neutrinoSliceIndex);
        this->AddSliceNeutrinoReconstruction(sliceList, neutrinoSliceIndex);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::RunAllHitsCosmicRayReconstruction(PfoList &parentCosmicRayPfos) const
{
    SliceList sliceList;
    this->CopyAllHitsToSingleSlice(sliceList);
    this->CosmicRayReconstruction(sliceList.front(), std::string());

    const PfoList *pParentCRPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_crParentListName, pParentCRPfoList));

    if (pParentCRPfoList)
        parentCosmicRayPfos = *pParentCRPfoList;

    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "ParentAlgorithm: cosmic-ray reconstuction done" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos) const
{
    int counter(0); // TODO - use tool here
    PfoList ambiguousParentPfos;

    for (const Pfo *const pParentCosmicRayPfo : parentCosmicRayPfos)
    {
        if (++counter % 4 == 0)
            ambiguousParentPfos.push_back(pParentCosmicRayPfo);
    }

    LArPfoHelper::GetAllConnectedPfos(ambiguousParentPfos, ambiguousPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::RemoveAmbiguousCosmicRayPfos(const PfoList &ambiguousPfos) const
{
    const PfoList *pParentCosmicRayPfos(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_crParentListName, pParentCosmicRayPfos));

    const PfoList *pCosmicRayDaughterPfos(nullptr);
    (void) PandoraContentApi::GetList(*this, m_crDaughterListName, pCosmicRayDaughterPfos);

    for (const Pfo *const pPfo : ambiguousPfos)
    {
        const std::string listName(std::find(pParentCosmicRayPfos->begin(), pParentCosmicRayPfos->end(), pPfo) != pParentCosmicRayPfos->end() ? m_crParentListName :
            (pCosmicRayDaughterPfos && std::find(pCosmicRayDaughterPfos->begin(), pCosmicRayDaughterPfos->end(), pPfo) != pCosmicRayDaughterPfos->end()) ? m_crDaughterListName : "");

        if (listName.empty())
        {
            std::cout << "ParentAlgorithm::RemoveAmbiguousCosmicRayPfos: parent/daughter particles list mismatch" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, listName));
    }

    this->RunAlgorithm(m_crListPruningAlgorithm);
    this->RunAlgorithm(m_crListMovingAlgorithm);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::RunFastReconstructionAndSlicing(SliceList &sliceList) const
{
    if (m_shouldPerformSlicing)
    {
        this->PerformSlicing(sliceList);
    }
    else
    {
        this->CopyAllHitsToSingleSlice(sliceList);
    }

    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "ParentAlgorithm: slicing done, nSlices " << sliceList.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::ReconstructSlices(const SliceList &sliceList, SliceIndexToPfoListMap &sliceToCosmicRayPfosMap,
    SliceIndexToPropertiesMap &sliceIndexToPropertiesMap) const
{
    unsigned int sliceIndex(0);

    for (const Slice &slice : sliceList)
    {
        const unsigned int thisSliceIndex(sliceIndex++);
        const std::string sliceIndexString(TypeToString(thisSliceIndex));
        SliceProperties sliceProperties; // TODO - fill with nu and cr properties

        this->NeutrinoReconstruction(slice, sliceIndexString);
        this->RunAlgorithm(m_nuListDeletionAlgorithm);

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ParentAlgorithm: neutrino reconstruction done for slice " << thisSliceIndex << std::endl;

        this->CosmicRayReconstruction(slice, sliceIndexString);

        const PfoList *pParentCRPfosInSlice(nullptr);
        (void) PandoraContentApi::GetList(*this, m_crParentListName, pParentCRPfosInSlice);

        if (pParentCRPfosInSlice)
            sliceToCosmicRayPfosMap[thisSliceIndex] = *pParentCRPfosInSlice;

        this->RunAlgorithm(m_crListMovingAlgorithm);

        sliceIndexToPropertiesMap[thisSliceIndex] = sliceProperties;

        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ParentAlgorithm: cosmic-ray reconstruction done for slice " << thisSliceIndex << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParentAlgorithm::GetNeutrinoSliceIndex(const SliceIndexToPropertiesMap &/*sliceIndexToPropertiesMap*/, unsigned int &neutrinoSliceIndex) const
{
    neutrinoSliceIndex = 0;
    return true; // TODO - examine properties here, using tool maybe
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::RemoveSliceCosmicRayReconstruction(const SliceIndexToPfoListMap &sliceToCosmicRayPfosMap, const unsigned int neutrinoSliceIndex) const
{
    const PfoList *pParentCosmicRayPfos(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_outputListPrefix + m_crParentListName, pParentCosmicRayPfos));

    const PfoList *pCosmicRayDaughterPfos(nullptr);
    (void) PandoraContentApi::GetList(*this, m_outputListPrefix + m_crDaughterListName, pCosmicRayDaughterPfos);

    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(sliceToCosmicRayPfosMap.at(neutrinoSliceIndex), allPfosToDelete);
PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &allPfosToDelete , "allPfosToDelete", RED);
PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const Pfo *const pPfo : allPfosToDelete)
    {
        const std::string listName(std::find(pParentCosmicRayPfos->begin(), pParentCosmicRayPfos->end(), pPfo) != pParentCosmicRayPfos->end() ? m_outputListPrefix + m_crParentListName :
            (pCosmicRayDaughterPfos && std::find(pCosmicRayDaughterPfos->begin(), pCosmicRayDaughterPfos->end(), pPfo) != pCosmicRayDaughterPfos->end()) ? m_outputListPrefix + m_crDaughterListName : "");

        if (listName.empty())
        {
            std::cout << "ParentAlgorithm::RemoveSliceCosmicRayReconstruction: parent/daughter particles list mismatch" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, listName));
    }

    this->RunAlgorithm(m_outputListPruningAlgorithm);

    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "ParentAlgorithm: cosmic-ray reconstruction removed for slice " << neutrinoSliceIndex << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::AddSliceNeutrinoReconstruction(const SliceList &sliceList, const unsigned int neutrinoSliceIndex) const
{
    this->NeutrinoReconstruction(sliceList.at(neutrinoSliceIndex), TypeToString(neutrinoSliceIndex));
    this->RunAlgorithm(m_nuListMovingAlgorithm);

    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "ParentAlgorithm: neutrino reconstruction applied to slice " << neutrinoSliceIndex << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::FastReconstruction() const
{
    this->RunTwoDClustering(std::string(), m_nuClusteringAlgorithm, false, m_nuTwoDAlgorithms);
    this->RunAlgorithms(m_nuThreeDAlgorithms);
    this->RunAlgorithms(m_nuThreeDHitAlgorithms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::CosmicRayReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const
{
    this->SaveTwoDCaloHitLists(slice, sliceIndexString);
    this->RunTwoDClustering(sliceIndexString, m_crTrackClusteringAlgorithm, false, m_crTwoDAlgorithms);
    this->RunAlgorithms(m_crThreeDAlgorithms);
    this->RunAlgorithm(m_crListPruningAlgorithm);
    this->RunTwoDClustering(sliceIndexString, m_crDeltaRayClusteringAlgorithm, true);
    this->RunAlgorithms(m_crDeltaRayAlgorithms);

    for (const HitType hitType : m_hitTypeList)
    {
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, (PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType))));
        this->RunAlgorithms(m_crTwoDRemnantAlgorithms);
    }

    this->RunAlgorithms(m_crThreeDRemnantAlgorithms);
    this->RunAlgorithms(m_crThreeDHitAlgorithms);
    this->RunAlgorithms(m_crVertexAlgorithms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::NeutrinoReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const
{
    this->SaveTwoDCaloHitLists(slice, sliceIndexString);
    this->RunTwoDClustering(sliceIndexString, m_nuClusteringAlgorithm, false, m_nuTwoDAlgorithms);
    this->RunAlgorithms(m_nuVertexAlgorithms);
    this->RunAlgorithms(m_nuThreeDAlgorithms);
    this->RunAlgorithms(m_nuTwoDMopUpAlgorithms);
    this->RunAlgorithms(m_nuThreeDHitAlgorithms);
    this->RunAlgorithms(m_nuThreeDMopUpAlgorithms);
    this->RunAlgorithms(m_nuNeutrinoAlgorithms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CRParentListName", m_crParentListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CRDaughterListName", m_crDaughterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputListPrefix", m_outputListPrefix));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "OutputListPruning", m_outputListPruningAlgorithm));

    TiXmlElement *const pCosmicRayXmlElement(xmlHandle.FirstChild("CosmicRayReconstruction").Element());

    if (nullptr != pCosmicRayXmlElement)
    {
        const TiXmlHandle crXmlHandle(pCosmicRayXmlElement);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "TwoDTrackClustering", m_crTrackClusteringAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "TwoDDeltaRayClustering", m_crDeltaRayClusteringAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "TwoDAlgorithms", m_crTwoDAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "ThreeDAlgorithms", m_crThreeDAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "DeltaRayAlgorithms", m_crDeltaRayAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "TwoDRemnantAlgorithms", m_crTwoDRemnantAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "ThreeDRemnantAlgorithms", m_crThreeDRemnantAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "ThreeDHitAlgorithms", m_crThreeDHitAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "VertexAlgorithms", m_crVertexAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "ListPruning", m_crListPruningAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "ListDeletion", m_crListDeletionAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "ListMoving", m_crListMovingAlgorithm));
    }

    TiXmlElement *const pNeutrinoXmlElement(xmlHandle.FirstChild("NeutrinoReconstruction").Element());

    if (nullptr != pNeutrinoXmlElement)
    {
        const TiXmlHandle nuXmlHandle(pNeutrinoXmlElement);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, nuXmlHandle,
            "TwoDClustering", m_nuClusteringAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "TwoDAlgorithms", m_nuTwoDAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "ThreeDAlgorithms", m_nuThreeDAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "ThreeDHitAlgorithms", m_nuThreeDHitAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "VertexAlgorithms", m_nuVertexAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "TwoDMopUpAlgorithms", m_nuTwoDMopUpAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "ThreeDMopUpAlgorithms", m_nuThreeDMopUpAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "NeutrinoAlgorithms", m_nuNeutrinoAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, nuXmlHandle,
            "ListDeletion", m_nuListDeletionAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, nuXmlHandle,
            "ListMoving", m_nuListMovingAlgorithm));
    }

    return ParentSlicingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
