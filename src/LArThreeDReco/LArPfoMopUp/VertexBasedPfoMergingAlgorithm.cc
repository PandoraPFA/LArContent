/**
 *  @file   LArContent/src/LArThreeDReco/LArPfoMopUp/VertexBasedPfoMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex based pfo merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArPfoMopUp/VertexBasedPfoMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexBasedPfoMergingAlgorithm::VertexBasedPfoMergingAlgorithm() :
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexTransverseDistance(1.5f),
    m_minVertexAssociatedHitTypes(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexBasedPfoMergingAlgorithm::Run()
{
    while (true)
    {
        PfoList vertexPfos, nonVertexPfos;
        this->GetInputPfos(vertexPfos, nonVertexPfos);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &vertexPfos, "vertexPfos", RED);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &nonVertexPfos, "nonVertexPfos", GREEN);
PandoraMonitoringApi::ViewEvent(this->GetPandora());
        PfoAssociationList pfoAssociationList;
        this->GetPfoAssociations(vertexPfos, nonVertexPfos, pfoAssociationList);

        const bool pfoMergeMade(this->ProcessPfoAssociations(pfoAssociationList));

        if (!pfoMergeMade)
            break;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoMergingAlgorithm::GetInputPfos(PfoList &vertexPfos, PfoList &nonVertexPfos) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex((1 == pVertexList->size()) ? *(pVertexList->begin()) : NULL);

    if (NULL == pVertex)
        return;

    StringVector listNames;
    listNames.push_back(m_trackPfoListName);
    listNames.push_back(m_showerPfoListName);

    for (StringVector::const_iterator iter = listNames.begin(), iterEnd = listNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList(NULL);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *iter, pPfoList))
            continue;

        for (PfoList::const_iterator pfoIter = pPfoList->begin(), pfoIterEnd = pPfoList->end(); pfoIter != pfoIterEnd; ++pfoIter)
        {
            Pfo *const pPfo(*pfoIter);
            PfoList &pfoTargetList(this->IsVertexAssociated(pPfo, pVertex) ? vertexPfos : nonVertexPfos);
            pfoTargetList.insert(pPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoMergingAlgorithm::GetPfoAssociations(const PfoList &/*vertexPfos*/, const PfoList &/*nonVertexPfos*/, PfoAssociationList &/*pfoAssociationList*/) const
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexBasedPfoMergingAlgorithm::ProcessPfoAssociations(const PfoAssociationList &/*pfoAssociationList*/) const
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexBasedPfoMergingAlgorithm::IsVertexAssociated(const Pfo *const pPfo, const Vertex *const pVertex) const
{
    if (VERTEX_3D != pVertex->GetVertexType())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    HitTypeSet hitTypeSet;
    const ClusterList &clusterList(pPfo->GetClusterList());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster(*iter);
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            continue;

        const CartesianVector vertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

        try
        {
            const LArPointingCluster pointingCluster(pCluster);

            if (LArPointingClusterHelper::IsNode(vertex2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
                LArPointingClusterHelper::IsNode(vertex2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance))
            {
                hitTypeSet.insert(hitType);
            }
        }
        catch (StatusCodeException &) {}
    }

    const unsigned int nVertexAssociatedHitTypes(hitTypeSet.size());
    return (nVertexAssociatedHitTypes >= m_minVertexAssociatedHitTypes);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexBasedPfoMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackPfoListName", m_trackPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerPfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexAssociatedHitTypes", m_minVertexAssociatedHitTypes));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
