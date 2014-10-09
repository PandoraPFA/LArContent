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
    const VertexList *pVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex((1 == pVertexList->size()) ? *(pVertexList->begin()) : NULL);

    if ((NULL == pVertex) || (VERTEX_3D != pVertex->GetVertexType()))
        return STATUS_CODE_SUCCESS;

    while (true)
    {
        PfoList vertexPfos, nonVertexPfos;
        this->GetInputPfos(pVertex, vertexPfos, nonVertexPfos);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &vertexPfos, "vertexPfos", RED);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &nonVertexPfos, "nonVertexPfos", GREEN);
PandoraMonitoringApi::ViewEvent(this->GetPandora());
        PfoAssociationList pfoAssociationList;
        this->GetPfoAssociations(pVertex, vertexPfos, nonVertexPfos, pfoAssociationList);

        std::sort(pfoAssociationList.begin(), pfoAssociationList.end());
        const bool pfoMergeMade(this->ProcessPfoAssociations(pfoAssociationList));

        if (!pfoMergeMade)
            break;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoMergingAlgorithm::GetInputPfos(const Vertex *const pVertex, PfoList &vertexPfos, PfoList &nonVertexPfos) const
{
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

void VertexBasedPfoMergingAlgorithm::GetPfoAssociations(const Vertex *const pVertex, const PfoList &vertexPfos, const PfoList &nonVertexPfos,
    PfoAssociationList &pfoAssociationList) const
{
    for (PfoList::const_iterator iter1 = vertexPfos.begin(), iter1End = vertexPfos.end(); iter1 != iter1End; ++iter1)
    {
        Pfo *const pVertexPfo(*iter1);

        for (PfoList::const_iterator iter2 = nonVertexPfos.begin(), iter2End = nonVertexPfos.end(); iter2 != iter2End; ++iter2)
        {
            Pfo *const pDaughterPfo(*iter2);

            try
            {
                const PfoAssociation pfoAssociation(this->GetPfoAssociation(pVertex, pVertexPfo, pDaughterPfo));
                pfoAssociationList.push_back(pfoAssociation);
            }
            catch (StatusCodeException &) {}
        }
    }
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

VertexBasedPfoMergingAlgorithm::PfoAssociation VertexBasedPfoMergingAlgorithm::GetPfoAssociation(const Vertex *const pVertex, Pfo *const pVertexPfo,
    Pfo *const pDaughterPfo) const
{
    const ClusterList &vertexClusterList(pVertexPfo->GetClusterList());
    const ClusterList &daughterClusterList(pDaughterPfo->GetClusterList());

    if ((vertexClusterList.size() != daughterClusterList.size()) || (3 != vertexClusterList.size()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    HitTypeToAssociationMap hitTypeToAssociationMap;

    for (ClusterList::const_iterator iter1 = vertexClusterList.begin(), iter1End = vertexClusterList.end(); iter1 != iter1End; ++iter1)
    {
        Cluster *const pVertexCluster(*iter1);
        const HitType vertexHitType(LArClusterHelper::GetClusterHitType(pVertexCluster));
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), vertexHitType));

        for (ClusterList::const_iterator iter2 = daughterClusterList.begin(), iter2End = daughterClusterList.end(); iter2 != iter2End; ++iter2)
        {
            Cluster *const pDaughterCluster(*iter2);
            const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));

            if (vertexHitType != daughterHitType)
                continue;

            const ClusterAssociation clusterAssociation(this->GetClusterAssociation(vertexPosition2D, pVertexCluster, pDaughterCluster));

            if (!hitTypeToAssociationMap.insert(HitTypeToAssociationMap::value_type(vertexHitType, clusterAssociation)).second)
                throw STATUS_CODE_FAILURE;
        }
    }

    HitTypeToAssociationMap::const_iterator iterU = hitTypeToAssociationMap.find(TPC_VIEW_U);
    HitTypeToAssociationMap::const_iterator iterV = hitTypeToAssociationMap.find(TPC_VIEW_V);
    HitTypeToAssociationMap::const_iterator iterW = hitTypeToAssociationMap.find(TPC_VIEW_W);

    if ((hitTypeToAssociationMap.end() == iterU) || (hitTypeToAssociationMap.end() == iterV) || (hitTypeToAssociationMap.end() == iterW))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return PfoAssociation(pVertexPfo, pDaughterPfo, iterU->second, iterV->second, iterW->second);
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::ClusterAssociation VertexBasedPfoMergingAlgorithm::GetClusterAssociation(const CartesianVector &vertexPosition2D,
    Cluster *const pVertexCluster, Cluster *const pDaughterCluster) const
{
    const ConeParameters coneParameters(pVertexCluster, vertexPosition2D, 0.8f); // TODO config
    const float boundedFraction(coneParameters.GetBoundedFraction(pDaughterCluster, 4.f)); // TODO config

    return ClusterAssociation(pVertexCluster, pDaughterCluster, boundedFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::ClusterAssociation::ClusterAssociation(Cluster *pVertexCluster, Cluster *pDaughterCluster, const float boundedFraction) :
    m_pVertexCluster(pVertexCluster),
    m_pDaughterCluster(pDaughterCluster),
    m_boundedFraction(boundedFraction)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::PfoAssociation::PfoAssociation(Pfo *pVertexPfo, Pfo *pDaughterPfo, const ClusterAssociation &clusterAssociationU,
        const ClusterAssociation &clusterAssociationV, const ClusterAssociation &clusterAssociationW) :
    m_pVertexPfo(pVertexPfo),
    m_pDaughterPfo(pDaughterPfo),
    m_clusterAssociationU(clusterAssociationU),
    m_clusterAssociationV(clusterAssociationV),
    m_clusterAssociationW(clusterAssociationW)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexBasedPfoMergingAlgorithm::PfoAssociation::GetMeanBoundedFraction() const
{
    return ((this->GetClusterAssociationU().GetBoundedFraction() + this->GetClusterAssociationV().GetBoundedFraction() + this->GetClusterAssociationW().GetBoundedFraction()) / 3.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexBasedPfoMergingAlgorithm::PfoAssociation::GetMaxBoundedFraction() const
{
    return std::max(this->GetClusterAssociationU().GetBoundedFraction(), std::max(this->GetClusterAssociationV().GetBoundedFraction(), this->GetClusterAssociationW().GetBoundedFraction()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexBasedPfoMergingAlgorithm::PfoAssociation::GetMinBoundedFraction() const
{
    return std::min(this->GetClusterAssociationU().GetBoundedFraction(), std::min(this->GetClusterAssociationV().GetBoundedFraction(), this->GetClusterAssociationW().GetBoundedFraction()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexBasedPfoMergingAlgorithm::PfoAssociation::operator< (const PfoAssociation &rhs) const
{
    return (this->GetMeanBoundedFraction() > rhs.GetMeanBoundedFraction());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::ConeParameters::ConeParameters(const Cluster *const pCluster, const CartesianVector &vertexPosition2D, const float coneAngleCentile) :
    m_pCluster(pCluster),
    m_apex(vertexPosition2D),
    m_direction(0.f, 0.f, 0.f),
    m_coneLength(0.f),
    m_coneCosHalfAngle(0.f)
{
    m_direction = this->GetDirectionEstimate();
    m_coneLength = this->GetConeLength();
    m_coneCosHalfAngle = this->GetCosHalfAngleEstimate(coneAngleCentile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexBasedPfoMergingAlgorithm::ConeParameters::GetBoundedFraction(const Cluster *const pDaughterCluster, const float coneLengthMultiplier) const
{
    unsigned int nMatchedHits(0);
    const OrderedCaloHitList &orderedCaloHitList(pDaughterCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CartesianVector &positionVector((*hIter)->GetPositionVector());

            if (m_direction.GetCosOpeningAngle(positionVector - m_apex) < m_coneCosHalfAngle)
                continue;

            if (m_direction.GetDotProduct(positionVector - m_apex) > coneLengthMultiplier * m_coneLength)
                continue;

            ++nMatchedHits;
        }
    }

    if (0 == pDaughterCluster->GetNCaloHits())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return (static_cast<float>(nMatchedHits) / static_cast<float>(pDaughterCluster->GetNCaloHits()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector VertexBasedPfoMergingAlgorithm::ConeParameters::GetDirectionEstimate() const
{
    const OrderedCaloHitList &orderedCaloHitList(m_pCluster->GetOrderedCaloHitList());
    float sumDxDz(0.f), sumDxDx(0.f);

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CartesianVector apexDisplacement((*hitIter)->GetPositionVector() - m_apex);
            sumDxDz += apexDisplacement.GetX() * apexDisplacement.GetZ();
            sumDxDx += apexDisplacement.GetX() * apexDisplacement.GetX();
        }
    }

    if (sumDxDx < std::numeric_limits<float>::epsilon())
        return CartesianVector(0.f, 0.f, 1.f);

    return CartesianVector(1.f, 0.f, sumDxDz / sumDxDx).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexBasedPfoMergingAlgorithm::ConeParameters::GetConeLength() const
{
    float maxProjectedLength(0.f);
    const OrderedCaloHitList &orderedCaloHitList(m_pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
            maxProjectedLength = std::max(maxProjectedLength, m_direction.GetDotProduct((*hitIter)->GetPositionVector() - m_apex));
    }

    return maxProjectedLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexBasedPfoMergingAlgorithm::ConeParameters::GetCosHalfAngleEstimate(const float coneAngleCentile) const
{
    FloatVector halfAngleValues;
    const OrderedCaloHitList &orderedCaloHitList(m_pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
            halfAngleValues.push_back(m_direction.GetOpeningAngle((*hitIter)->GetPositionVector() - m_apex));
    }

    std::sort(halfAngleValues.begin(), halfAngleValues.end());

    if (halfAngleValues.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const unsigned int halfAngleBin(coneAngleCentile * halfAngleValues.size());
    return std::cos(halfAngleValues.at(halfAngleBin));
}

//------------------------------------------------------------------------------------------------------------------------------------------
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
