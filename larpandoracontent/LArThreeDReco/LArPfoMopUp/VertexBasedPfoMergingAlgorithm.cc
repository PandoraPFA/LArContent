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
#include "LArHelpers/LArPfoHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArPfoMopUp/VertexBasedPfoMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexBasedPfoMergingAlgorithm::VertexBasedPfoMergingAlgorithm() :
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexTransverseDistance(1.5f),
    m_minVertexAssociatedHitTypes(2),
    m_coneAngleCentile(0.8f),
    m_maxConeCosHalfAngle(0.95f),
    m_maxConeLengthMultiplier(3.f),
    m_directionTanAngle(1.732f),
    m_directionApexShift(0.333f),
    m_meanBoundedFractionCut(0.6f),
    m_maxBoundedFractionCut(0.7f),
    m_minBoundedFractionCut(0.3f),
    m_minConsistentDirections(2),
    m_minConsistentDirectionsTrack(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexBasedPfoMergingAlgorithm::Run()
{
    const VertexList *pVertexList = nullptr;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    const Vertex *const pSelectedVertex((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

    if (!pSelectedVertex)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexBasedPfoMerging: unable to find vertex in current list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    while (true)
    {
        PfoList vertexPfos, nonVertexPfos;
        this->GetInputPfos(pSelectedVertex, vertexPfos, nonVertexPfos);

        PfoAssociationList pfoAssociationList;
        this->GetPfoAssociations(pSelectedVertex, vertexPfos, nonVertexPfos, pfoAssociationList);

        std::sort(pfoAssociationList.begin(), pfoAssociationList.end());
        const bool pfoMergeMade(this->ProcessPfoAssociations(pfoAssociationList));

        if (!pfoMergeMade)
            break;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexBasedPfoMergingAlgorithm::IsVertexAssociated(const CartesianVector &vertex2D, const LArPointingCluster &pointingCluster) const
{
    return (LArPointingClusterHelper::IsNode(vertex2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(vertex2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::PfoAssociation VertexBasedPfoMergingAlgorithm::GetPfoAssociation(const Pfo *const pVertexPfo, const Pfo *const pDaughterPfo,
    HitTypeToAssociationMap &hitTypeToAssociationMap) const
{
    if ((pVertexPfo->GetClusterList().size() != pDaughterPfo->GetClusterList().size()) || (3 != pVertexPfo->GetClusterList().size()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return PfoAssociation(pVertexPfo, pDaughterPfo, hitTypeToAssociationMap.at(TPC_VIEW_U), hitTypeToAssociationMap.at(TPC_VIEW_V), hitTypeToAssociationMap.at(TPC_VIEW_W));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoMergingAlgorithm::GetInputPfos(const Vertex *const pVertex, PfoList &vertexPfos, PfoList &nonVertexPfos) const
{
    StringVector listNames;
    listNames.push_back(m_trackPfoListName);
    listNames.push_back(m_showerPfoListName);

    for (const std::string &listName : listNames)
    {
        const PfoList *pPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pPfoList))
            continue;

        for (const Pfo *const pPfo : *pPfoList)
        {
            PfoList &pfoTargetList(this->IsVertexAssociated(pPfo, pVertex) ? vertexPfos : nonVertexPfos);
            pfoTargetList.insert(pPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexBasedPfoMergingAlgorithm::IsVertexAssociated(const Pfo *const pPfo, const Vertex *const pVertex) const
{
    if (VERTEX_3D != pVertex->GetVertexType())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    HitTypeSet hitTypeSet;

    for (const Cluster *const pCluster : pPfo->GetClusterList())
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            continue;

        const CartesianVector vertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

        try
        {
            const LArPointingCluster pointingCluster(pCluster);

            if (this->IsVertexAssociated(vertex2D, pointingCluster))
                hitTypeSet.insert(hitType);
        }
        catch (StatusCodeException &) {}
    }

    const unsigned int nVertexAssociatedHitTypes(hitTypeSet.size());
    return (nVertexAssociatedHitTypes >= m_minVertexAssociatedHitTypes);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoMergingAlgorithm::GetPfoAssociations(const Vertex *const pVertex, const PfoList &vertexPfos, const PfoList &nonVertexPfos,
    PfoAssociationList &pfoAssociationList) const
{
    for (const Pfo *const pVertexPfo : vertexPfos)
    {
        for (const Pfo *const pDaughterPfo : nonVertexPfos)
        {
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

VertexBasedPfoMergingAlgorithm::PfoAssociation VertexBasedPfoMergingAlgorithm::GetPfoAssociation(const Vertex *const pVertex, const Pfo *const pVertexPfo,
    const Pfo *const pDaughterPfo) const
{
    if (pVertexPfo->GetClusterList().empty() || pDaughterPfo->GetClusterList().empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    HitTypeToAssociationMap hitTypeToAssociationMap;

    for (const Cluster *const pVertexCluster : pVertexPfo->GetClusterList())
    {
        const HitType vertexHitType(LArClusterHelper::GetClusterHitType(pVertexCluster));

        for (const Cluster *const pDaughterCluster : pDaughterPfo->GetClusterList())
        {
            const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));

            if (vertexHitType != daughterHitType)
                continue;

            const ClusterAssociation clusterAssociation(this->GetClusterAssociation(pVertex, pVertexCluster, pDaughterCluster));
            hitTypeToAssociationMap[vertexHitType] = clusterAssociation;
        }
    }

    return this->GetPfoAssociation(pVertexPfo, pDaughterPfo, hitTypeToAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::ClusterAssociation VertexBasedPfoMergingAlgorithm::GetClusterAssociation(const Vertex *const pVertex,
    const Cluster *const pVertexCluster, const Cluster *const pDaughterCluster) const
{
    const HitType vertexHitType(LArClusterHelper::GetClusterHitType(pVertexCluster));
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), vertexHitType));

    const ConeParameters coneParameters(pVertexCluster, vertexPosition2D, m_coneAngleCentile, m_maxConeCosHalfAngle);
    const float boundedFraction(coneParameters.GetBoundedFraction(pDaughterCluster, m_maxConeLengthMultiplier));

    const LArVertexHelper::ClusterDirection vertexClusterDirection(LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex,
        pVertexCluster, m_directionTanAngle, m_directionApexShift));
    const LArVertexHelper::ClusterDirection daughterClusterDirection(LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex,
        pDaughterCluster, m_directionTanAngle, m_directionApexShift));
    const bool isConsistentDirection(vertexClusterDirection == daughterClusterDirection);

    return ClusterAssociation(pVertexCluster, pDaughterCluster, boundedFraction, isConsistentDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexBasedPfoMergingAlgorithm::ProcessPfoAssociations(const PfoAssociationList &pfoAssociationList) const
{
    const PfoList *pTrackPfoList(nullptr);
    (void) PandoraContentApi::GetList(*this, m_trackPfoListName, pTrackPfoList);

    for (const PfoAssociation &pfoAssociation : pfoAssociationList)
    {
        if ((pfoAssociation.GetMeanBoundedFraction() < m_meanBoundedFractionCut) ||
            (pfoAssociation.GetMaxBoundedFraction() < m_maxBoundedFractionCut) ||
            (pfoAssociation.GetMinBoundedFraction() < m_minBoundedFractionCut) ||
            (pfoAssociation.GetNConsistentDirections() < m_minConsistentDirections))
        {
            continue;
        }

        if (pTrackPfoList)
        {
            if ((pTrackPfoList->count(pfoAssociation.GetVertexPfo()) > 0) && (pTrackPfoList->count(pfoAssociation.GetDaughterPfo()) > 0))
            {
                continue;
            }

            if (((pTrackPfoList->count(pfoAssociation.GetVertexPfo()) > 0) || (pTrackPfoList->count(pfoAssociation.GetDaughterPfo()) > 0)) &&
                (pfoAssociation.GetNConsistentDirections() < m_minConsistentDirectionsTrack))
            {
                continue;
            }
        }

        this->MergePfos(pfoAssociation);
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoMergingAlgorithm::MergePfos(const PfoAssociation &pfoAssociation) const
{
    const PfoList *pTrackPfoList(nullptr), *pShowerPfoList(nullptr);
    (void) PandoraContentApi::GetList(*this, m_trackPfoListName, pTrackPfoList);
    (void) PandoraContentApi::GetList(*this, m_showerPfoListName, pShowerPfoList);

    if (!pTrackPfoList && !pShowerPfoList)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const Pfo *const pVertexPfo(pfoAssociation.GetVertexPfo());
    const bool isvertexTrack(pTrackPfoList && (pTrackPfoList->count(pVertexPfo) > 0));
    const Pfo *pDaughterPfo(pfoAssociation.GetDaughterPfo());
    const bool isDaughterShower(pShowerPfoList && (pShowerPfoList->count(pDaughterPfo) > 0));

    this->MergeAndDeletePfos(pVertexPfo, pDaughterPfo);

    if (isvertexTrack && isDaughterShower)
    {
        PfoList vertexPfoList;
        vertexPfoList.insert(pVertexPfo);

        PandoraContentApi::ParticleFlowObject::Metadata metadata;
        metadata.m_particleId = E_MINUS;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AlterMetadata(*this, pVertexPfo, metadata));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, vertexPfoList));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoMergingAlgorithm::MergeAndDeletePfos(const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete) const
{
    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, this->GetListName(pPfoToDelete)));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoToEnlarge, pDaughterPfo));
    }

    for (const  Vertex *const pDaughterVertex : daughterVertices)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pDaughterVertex, this->GetListName(pDaughterVertex)));
    }

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *pParentCluster(this->GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster,
                this->GetListName(pParentCluster), this->GetListName(pDaughterCluster)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *VertexBasedPfoMergingAlgorithm::GetParentCluster(const ClusterList &clusterList, const HitType hitType) const
{
    unsigned int mostHits(0);
    const Cluster *pBestParentCluster(nullptr);

    for (const Cluster *const pParentCluster : clusterList)
    {
        if (hitType != LArClusterHelper::GetClusterHitType(pParentCluster))
            continue;

        const unsigned int nParentHits(pParentCluster->GetNCaloHits());

        if (nParentHits > mostHits)
        {
            mostHits = nParentHits;
            pBestParentCluster = pParentCluster;
        }
    }

    return pBestParentCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const std::string VertexBasedPfoMergingAlgorithm::GetListName(const T *const pT) const
{
    std::string currentListName;
    const std::MANAGED_CONTAINER<const T*> *pCurrentList(nullptr);
    (void) PandoraContentApi::GetCurrentList(*this, pCurrentList, currentListName);

    if (pCurrentList && (pCurrentList->count(pT)))
        return currentListName;

    for (const std::string &listName : m_daughterListNames)
    {
        const std::MANAGED_CONTAINER<const T*> *pList(nullptr);
        (void) PandoraContentApi::GetList(*this, listName, pList);

        if (pList && (pList->count(pT)))
            return listName;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::ClusterAssociation::ClusterAssociation() :
    m_pVertexCluster(nullptr),
    m_pDaughterCluster(nullptr),
    m_boundedFraction(0.f),
    m_isConsistentDirection(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::ClusterAssociation::ClusterAssociation(const Cluster *const pVertexCluster, const Cluster *const pDaughterCluster,
        const float boundedFraction, const bool isConsistentDirection) :
    m_pVertexCluster(pVertexCluster),
    m_pDaughterCluster(pDaughterCluster),
    m_boundedFraction(boundedFraction),
    m_isConsistentDirection(isConsistentDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::PfoAssociation::PfoAssociation(const Pfo *const pVertexPfo, const Pfo *const pDaughterPfo, const ClusterAssociation &clusterAssociationU,
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

unsigned int VertexBasedPfoMergingAlgorithm::PfoAssociation::GetNConsistentDirections() const
{
    unsigned int nConsistentDirections(0);

    if (this->GetClusterAssociationU().IsConsistentDirection())
        ++nConsistentDirections;

    if (this->GetClusterAssociationV().IsConsistentDirection())
        ++nConsistentDirections;

    if (this->GetClusterAssociationW().IsConsistentDirection())
        ++nConsistentDirections;

    return nConsistentDirections;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexBasedPfoMergingAlgorithm::PfoAssociation::operator< (const PfoAssociation &rhs) const
{
    if (std::fabs(this->GetMeanBoundedFraction() - rhs.GetMeanBoundedFraction()) > std::numeric_limits<float>::epsilon())
        return (this->GetMeanBoundedFraction() > rhs.GetMeanBoundedFraction());

    if (m_pVertexPfo != rhs.m_pVertexPfo)
        return LArPfoHelper::SortByNHits(m_pVertexPfo, rhs.m_pVertexPfo);

    return LArPfoHelper::SortByNHits(m_pDaughterPfo, rhs.m_pDaughterPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoMergingAlgorithm::ConeParameters::ConeParameters(const Cluster *const pCluster, const CartesianVector &vertexPosition2D,
        const float coneAngleCentile, const float maxCosHalfAngle) :
    m_pCluster(pCluster),
    m_apex(vertexPosition2D),
    m_direction(0.f, 0.f, 0.f),
    m_coneLength(0.f),
    m_coneCosHalfAngle(0.f)
{
    m_direction = this->GetDirectionEstimate();
    m_coneLength = this->GetSignedConeLength();

    // ATTN Compensate for coordinate shift when fitting with vertex constraint (cleaner way to do this?)
    if (m_coneLength < std::numeric_limits<float>::epsilon())
    {
        m_direction = m_direction * -1.f;
        m_coneLength = std::fabs(m_coneLength);
    }

    m_coneCosHalfAngle = std::min(maxCosHalfAngle, this->GetCosHalfAngleEstimate(coneAngleCentile));
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

float VertexBasedPfoMergingAlgorithm::ConeParameters::GetSignedConeLength() const
{
    float maxProjectedLength(0.f);
    const OrderedCaloHitList &orderedCaloHitList(m_pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const float projectedLength(m_direction.GetDotProduct((*hitIter)->GetPositionVector() - m_apex));

            if (std::fabs(projectedLength) > std::fabs(maxProjectedLength))
                maxProjectedLength = projectedLength;
        }
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

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "DaughterListNames", m_daughterListNames));

    m_daughterListNames.push_back(m_trackPfoListName);
    m_daughterListNames.push_back(m_showerPfoListName);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexAssociatedHitTypes", m_minVertexAssociatedHitTypes));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ConeAngleCentile", m_coneAngleCentile));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxConeCosHalfAngle", m_maxConeCosHalfAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxConeLengthMultiplier", m_maxConeLengthMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DirectionTanAngle", m_directionTanAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DirectionApexShift", m_directionApexShift));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MeanBoundedFractionCut", m_meanBoundedFractionCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxBoundedFractionCut", m_maxBoundedFractionCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinBoundedFractionCut", m_minBoundedFractionCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinConsistentDirections", m_minConsistentDirections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinConsistentDirectionsTrack", m_minConsistentDirectionsTrack));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
