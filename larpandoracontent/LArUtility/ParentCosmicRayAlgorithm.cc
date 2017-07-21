/**
 *  @file   larpandoracontent/LArUtility/ParentCosmicRayAlgorithm.cc
 *
 *  @brief  Implementation of the parent cosmic ray algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ParentCosmicRayAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ParentCosmicRayAlgorithm::Run()
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
        this->CosmicRayReconstruction(slice, sliceIndexString);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentCosmicRayAlgorithm::FastReconstruction() const
{
    std::cout << "Fast reconstruction not implemented for dedicated cosmic-ray reconstruction." << std::endl;
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentCosmicRayAlgorithm::CosmicRayReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const
{
    this->SaveTwoDCaloHitLists(slice, sliceIndexString);
    this->RunTwoDClustering(sliceIndexString, m_trackClusteringAlgorithm, false, m_twoDAlgorithms);
    this->RunAlgorithms(m_threeDAlgorithms);
    this->RunAlgorithm(m_listPruningAlgorithm);
    this->RunTwoDClustering(sliceIndexString, m_deltaRayClusteringAlgorithm, true);
    this->RunAlgorithms(m_deltaRayAlgorithms);

    for (const HitType hitType : m_hitTypeList)
    {
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, (PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType))));
        this->RunAlgorithms(m_twoDRemnantAlgorithms);
    }

    this->RunAlgorithms(m_threeDRemnantAlgorithms);
    this->RunAlgorithms(m_threeDHitAlgorithms);
    this->RunAlgorithms(m_vertexAlgorithms);
    this->RunAlgorithm(m_listMovingAlgorithm);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentCosmicRayAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "TwoDTrackClustering", m_trackClusteringAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "TwoDDeltaRayClustering", m_deltaRayClusteringAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDAlgorithms", m_twoDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDAlgorithms", m_threeDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "DeltaRayAlgorithms", m_deltaRayAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDRemnantAlgorithms", m_twoDRemnantAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDRemnantAlgorithms", m_threeDRemnantAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDHitAlgorithms", m_threeDHitAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "VertexAlgorithms", m_vertexAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ListPruning", m_listPruningAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ListMoving", m_listMovingAlgorithm));

    return ParentSlicingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
