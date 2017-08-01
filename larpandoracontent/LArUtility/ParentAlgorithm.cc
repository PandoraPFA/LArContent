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
    // Run cosmic-ray reconstruction for all hits
    SliceList crSliceList;
    this->CopyAllHitsToSingleSlice(crSliceList);
    this->CosmicRayReconstruction(crSliceList.front(), std::string());
std::cout << "CRReco -- done" << std::endl;

    const PfoList *pParentCRPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_crParentListName, pParentCRPfoList));

    if (!pParentCRPfoList || pParentCRPfoList->empty())
    {
        std::cout << "ParentAlgorithm: No cosmic-ray muons reconstructed " << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    // Tag and delete all products associated with ambiguous cosmic rays, moving only clear cosmic rays to the output lists
    int counter(0); // TODO
    PfoList pfosToDelete;

    for (const Pfo *const pCosmicRayPfo : *pParentCRPfoList)
    {
        if (++counter % 2 == 0)
            pfosToDelete.push_back(pCosmicRayPfo);
    }

    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(pfosToDelete, allPfosToDelete);

    const PfoList *pCRDaughterPfoList(nullptr);
    (void) PandoraContentApi::GetList(*this, m_crDaughterListName, pCRDaughterPfoList);

    for (const Pfo *const pPfo : allPfosToDelete)
    {
        const std::string listName(std::find(pParentCRPfoList->begin(), pParentCRPfoList->end(), pPfo) != pParentCRPfoList->end() ? m_crParentListName :
            (pCRDaughterPfoList && std::find(pCRDaughterPfoList->begin(), pCRDaughterPfoList->end(), pPfo) != pCRDaughterPfoList->end()) ? m_crDaughterListName : "");

        if (listName.empty())
        {
            std::cout << "ParentAlgorithm: cosmic-ray parent/daughter particles list mismatch " << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, listName));
    }

    this->RunAlgorithm(m_crListPruningAlgorithm);
    this->RunAlgorithm(m_crListMovingAlgorithm);
std::cout << "CRTagging -- done" << std::endl;

    // Perform slicing, using simple fast reconstruction of the remaining available hits
    SliceList nuSliceList;

    if (m_shouldPerformSlicing)
    {
        this->PerformSlicing(nuSliceList);
    }
    else
    {
        this->CopyAllHitsToSingleSlice(nuSliceList);
    }

    if (nuSliceList.empty())
    {
        std::cout << "ParentAlgorithm: No slices found, so neutrino reconstructed " << std::endl;
        return STATUS_CODE_SUCCESS;
    }
std::cout << "Slicing -- done, nSlices " << nuSliceList.size() << std::endl;

    // Run cosmic-ray and neutrino reconstruction for each slice, keeping just the cosmic-ray output for now
    unsigned int sliceIndex(0);

    typedef std::map<unsigned int, PfoList> SliceToParentPfosMap;
    SliceToParentPfosMap sliceToParentCRPfosMap;

    for (const Slice &slice : nuSliceList)
    {
        const unsigned int thisSliceIndex(sliceIndex++);
        const std::string sliceIndexString(TypeToString(thisSliceIndex));

        this->NeutrinoReconstruction(slice, sliceIndexString);
        this->RunAlgorithm(m_nuListDeletionAlgorithm);
std::cout << "Neutrino Reco " << sliceIndexString << " -- done" << std::endl;

        this->CosmicRayReconstruction(slice, sliceIndexString);

        const PfoList *pParentCRPfosInSlice(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_crParentListName, pParentCRPfosInSlice));

        if (pParentCRPfosInSlice)
            sliceToParentCRPfosMap[thisSliceIndex] = *pParentCRPfosInSlice;

        this->RunAlgorithm(m_crListMovingAlgorithm);
std::cout << "CR Reco " << sliceIndexString << " -- done" << std::endl;
    }

    // For chosen nu slice, replace cosmic-ray reconstruction output with the neutrino reconstruction output
    const int neutrinoSliceIndex(0); // TODO

    const PfoList *pAllParentCRPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_outputListPrefix + m_crParentListName, pAllParentCRPfoList));

    const PfoList *pAllCRDaughterPfoList(nullptr);
    (void) PandoraContentApi::GetList(*this, m_outputListPrefix + m_crDaughterListName, pAllCRDaughterPfoList);

    PfoList allSlicePfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(sliceToParentCRPfosMap.at(neutrinoSliceIndex), allSlicePfosToDelete);
PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &allSlicePfosToDelete , "allSlicePfosToDelete", RED);
PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const Pfo *const pPfo : allSlicePfosToDelete)
    {
        const std::string listName(std::find(pAllParentCRPfoList->begin(), pAllParentCRPfoList->end(), pPfo) != pAllParentCRPfoList->end() ? m_outputListPrefix + m_crParentListName :
            (pAllCRDaughterPfoList && std::find(pAllCRDaughterPfoList->begin(), pAllCRDaughterPfoList->end(), pPfo) != pAllCRDaughterPfoList->end()) ? m_outputListPrefix + m_crDaughterListName : "");

        if (listName.empty())
        {
            std::cout << "ParentAlgorithm: cosmic-ray parent/daughter particles list mismatch " << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, listName));
    }

    this->RunAlgorithm(m_outputListPruningAlgorithm);
std::cout << "NeutrinoTagging " << neutrinoSliceIndex << " -- done" << std::endl;

    this->NeutrinoReconstruction(nuSliceList.at(neutrinoSliceIndex), TypeToString(neutrinoSliceIndex));
    this->RunAlgorithm(m_nuListMovingAlgorithm);
std::cout << "NeutrinoReco -- done" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::FastReconstruction() const
{
    // OR, leave blank and slice remainder of initial CR reco
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
