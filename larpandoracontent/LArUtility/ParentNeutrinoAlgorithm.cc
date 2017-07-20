/**
 *  @file   larpandoracontent/LArUtility/ParentNeutrinoAlgorithm.cc
 *
 *  @brief  Implementation of the parent neutrino algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ParentNeutrinoAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ParentNeutrinoAlgorithm::Run()
{
    SliceList sliceList;

    if (m_shouldPerformSlicing)
    {
        this->PerformSlicing(sliceList);
    }
    else
    {
        this->CopyAllHitsToSingleSlice(sliceList);
    }

    unsigned int sliceIndex(0);

    for (const Slice &slice : sliceList)
    {
        const std::string sliceIndexString(TypeToString(sliceIndex++));
        this->NeutrinoReconstruction(slice, sliceIndexString);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentNeutrinoAlgorithm::FastReconstruction() const
{
    this->RunTwoDClustering(std::string(), m_clusteringAlgorithm, false, m_twoDAlgorithms);
    this->RunAlgorithms(m_threeDAlgorithms);
    this->RunAlgorithms(m_threeDHitAlgorithms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentNeutrinoAlgorithm::NeutrinoReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const
{
    for (const HitType hitType : m_hitTypeList)
    {
        const CaloHitList &caloHitList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU : (TPC_VIEW_V == hitType) ? slice.m_caloHitListV : slice.m_caloHitListW);
        const std::string workingCaloHitListName(m_caloHitListNames.at(hitType) + sliceIndexString);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList, workingCaloHitListName));
    }

    this->RunTwoDClustering(sliceIndexString, m_clusteringAlgorithm, false, m_twoDAlgorithms);
    this->RunAlgorithms(m_vertexAlgorithms);
    this->RunAlgorithms(m_threeDAlgorithms);
    this->RunAlgorithms(m_twoDMopUpAlgorithms);
    this->RunAlgorithms(m_threeDHitAlgorithms);
    this->RunAlgorithms(m_threeDMopUpAlgorithms);
    this->RunAlgorithms(m_neutrinoAlgorithms);
    this->RunAlgorithm(m_listMovingAlgorithm);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentNeutrinoAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "TwoDClustering", m_clusteringAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDAlgorithms", m_twoDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDAlgorithms", m_threeDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDHitAlgorithms", m_threeDHitAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "VertexAlgorithms", m_vertexAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDMopUpAlgorithms", m_twoDMopUpAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDMopUpAlgorithms", m_threeDMopUpAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "NeutrinoAlgorithms", m_neutrinoAlgorithms));

    return ParentSlicingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
