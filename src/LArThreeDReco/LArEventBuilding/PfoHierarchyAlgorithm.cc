/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/PfoHierarchyAlgorithm.cc
 * 
 *  @brief  Implementation of the pfo hierarchy algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPfoHelper.h"

#include "LArThreeDReco/LArEventBuilding/PfoHierarchyAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoHierarchyAlgorithm::PfoHierarchyAlgorithm() :
    m_halfWindowLayers(20)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoHierarchyAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoListName,
        pPfoList));

    if (NULL == pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoVertexBuildingAlgorithm: pfo list unavailable." << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // TODO
    const Vertex *pNeutrinoVertex(NULL);
    PfoInfoMap pfoInfoMap;

    for (PfoRelationToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        (*iter)->Run(this, pNeutrinoVertex, pfoInfoMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PfoHierarchyAlgorithm::PfoInfo::PfoInfo(const pandora::ParticleFlowObject *const pPfo, const unsigned int halfWindowLayers,
        const float layerPitch) :
    m_pThisPfo(pPfo),
    m_pCluster3D(NULL),
    m_pVertex3D(LArPfoHelper::GetVertex(pPfo)),
    m_pSlidingFitResult3D(NULL),
    m_isNeutrinoVertexAssociated(false),
    m_pParentPfo(NULL)
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    if (1 != clusterList3D.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_pCluster3D = *(clusterList3D.begin());
    m_pSlidingFitResult3D = new ThreeDSlidingFitResult(m_pCluster3D, halfWindowLayers, layerPitch);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoHierarchyAlgorithm::PfoInfo::~PfoInfo()
{
    delete m_pSlidingFitResult3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::SetNeutrinoVertexAssociation(const bool isNeutrinoVertexAssociated)
{
    m_isNeutrinoVertexAssociated = isNeutrinoVertexAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::SetParentPfo(const pandora::ParticleFlowObject *const pParentPfo)
{
    if (m_pParentPfo)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    m_pParentPfo = pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::RemoveParentPfo()
{
    m_pParentPfo = NULL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::AddDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo)
{
    if (!m_daughterPfoList.insert(pDaughterPfo).second)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::RemoveDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo)
{
    if (!m_daughterPfoList.erase(pDaughterPfo))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoHierarchyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "PfoRelationTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        PfoRelationTool *const pPfoRelationTool(dynamic_cast<PfoRelationTool*>(*iter));

        if (NULL == pPfoRelationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pPfoRelationTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
