/**
 *  @file   larpandoracontent/LArUtility/ParentBaseAlgorithm.cc
 *
 *  @brief  Implementation of the parent base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ParentBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ParentBaseAlgorithm::ParentBaseAlgorithm()
{
    m_hitTypeList.push_back(TPC_VIEW_U);
    m_hitTypeList.push_back(TPC_VIEW_V);
    m_hitTypeList.push_back(TPC_VIEW_W);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParentBaseAlgorithm::~ParentBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentBaseAlgorithm::Initialize()
{
    m_caloHitListNames[TPC_VIEW_U] = m_caloHitListNameU;
    m_caloHitListNames[TPC_VIEW_V] = m_caloHitListNameV;
    m_caloHitListNames[TPC_VIEW_W] = m_caloHitListNameW;

    m_clusterListNames[TPC_VIEW_U] = m_clusterListNameU;
    m_clusterListNames[TPC_VIEW_V] = m_clusterListNameV;
    m_clusterListNames[TPC_VIEW_W] = m_clusterListNameW;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentBaseAlgorithm::RunAlgorithm(const std::string &algorithmName) const
{
    this->RunAlgorithms(StringVector(1, algorithmName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentBaseAlgorithm::RunAlgorithms(const StringVector &algorithmNames) const
{
    for (const std::string &algorithmName : algorithmNames)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListNameU", m_caloHitListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListNameV", m_caloHitListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListNameW", m_caloHitListNameW));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListNameU", m_clusterListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListNameV", m_clusterListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListNameW", m_clusterListNameW));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
