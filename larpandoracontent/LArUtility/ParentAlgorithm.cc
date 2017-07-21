/**
 *  @file   larpandoracontent/LArUtility/ParentAlgorithm.cc
 *
 *  @brief  Implementation of the parent algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ParentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ParentAlgorithm::Run()
{
    SliceList crSliceList;
    this->CopyAllHitsToSingleSlice(crSliceList);
std::cout << "Full CRReco" << std::endl;
    this->CosmicRayReconstruction(crSliceList.front(), std::string());
    this->RunAlgorithm(m_crListDeletionAlgorithm);
std::cout << "Full CRReco -- done" << std::endl;
    // Move unambiguous CRs into final list, delete/prune everything else, make new hit list from available hits

    SliceList nuSliceList;

    if (m_shouldPerformSlicing)
    {
std::cout << "Slicing" << std::endl;
        this->PerformSlicing(nuSliceList);
std::cout << "Slicing -- done, nSlices " << nuSliceList.size() << std::endl;
    }
    else
    {
        this->CopyAllHitsToSingleSlice(nuSliceList);
    }

    unsigned int sliceIndex(0);

    for (const Slice &slice : nuSliceList)
    {
        const std::string sliceIndexString(TypeToString(sliceIndex++));
std::cout << "Neutrino Reco " << sliceIndexString << std::endl;
        this->NeutrinoReconstruction(slice, sliceIndexString);
        this->RunAlgorithm(m_nuListDeletionAlgorithm);
std::cout << "Neutrino Reco " << sliceIndexString << " -- done" << std::endl;
        //this->RunAlgorithm(m_nuListMovingAlgorithm);

std::cout << "CR Reco " << sliceIndexString << std::endl;
        this->CosmicRayReconstruction(slice, sliceIndexString);
        //this->RunAlgorithm(m_crListDeletionAlgorithm);
        this->RunAlgorithm(m_crListMovingAlgorithm);
std::cout << "CR Reco " << sliceIndexString << " -- done" << std::endl;

        // Store map of top-level PFOs in each slice, for later access
    }

    // For chosen nu slice, delete all PFOs stored in relevant entry of map, above. Get slice index and re-run nu reco and list moving

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
