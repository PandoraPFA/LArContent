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
    this->CosmicRayReconstruction(crSliceList.front(), "_hack");

    // Delete everything for now
    this->RunAlgorithm(m_listDeletionAlgorithm);

    SliceList nuSliceList;

    if (m_shouldPerformSlicing)
    {
        this->PerformSlicing(nuSliceList);
    }
    else
    {
        this->CopyAllHitsToSingleSlice(nuSliceList);
    }

    unsigned int sliceIndex(0);

    for (const Slice &slice : nuSliceList)
    {
        const std::string sliceIndexString(TypeToString(sliceIndex++));
        this->NeutrinoReconstruction(slice, sliceIndexString);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::FastReconstruction() const
{
    this->RunTwoDClustering(std::string(), m_crTrackClusteringAlgorithm, false, m_crTwoDAlgorithms);
    this->RunAlgorithms(m_crThreeDAlgorithms);
    this->RunAlgorithms(m_crThreeDHitAlgorithms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::CosmicRayReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const
{
    for (const HitType hitType : m_hitTypeList)
    {
        const CaloHitList &caloHitList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU : (TPC_VIEW_V == hitType) ? slice.m_caloHitListV : slice.m_caloHitListW);
        const std::string workingCaloHitListName(m_caloHitListNames.at(hitType) + sliceIndexString);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList, workingCaloHitListName));
    }

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
    //this->RunAlgorithm(m_listMovingAlgorithm);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::NeutrinoReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const
{
    for (const HitType hitType : m_hitTypeList)
    {
        const CaloHitList &caloHitList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU : (TPC_VIEW_V == hitType) ? slice.m_caloHitListV : slice.m_caloHitListW);
        const std::string workingCaloHitListName(m_caloHitListNames.at(hitType) + sliceIndexString);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList, workingCaloHitListName));
    }

    this->RunTwoDClustering(sliceIndexString, m_nuClusteringAlgorithm, false, m_nuTwoDAlgorithms);
    this->RunAlgorithms(m_nuVertexAlgorithms);
    this->RunAlgorithms(m_nuThreeDAlgorithms);
    this->RunAlgorithms(m_nuTwoDMopUpAlgorithms);
    this->RunAlgorithms(m_nuThreeDHitAlgorithms);
    this->RunAlgorithms(m_nuThreeDMopUpAlgorithms);
    this->RunAlgorithms(m_nuNeutrinoAlgorithms);
    this->RunAlgorithm(m_listMovingAlgorithm);
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

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "ListPruning", m_crListPruningAlgorithm));

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
    }

    return ParentSlicingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
