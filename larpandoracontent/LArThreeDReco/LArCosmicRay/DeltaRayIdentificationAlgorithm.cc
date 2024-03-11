/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray identification algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DeltaRayIdentificationAlgorithm::DeltaRayIdentificationAlgorithm() :
    m_distanceForMatching(3.f),
    m_minParentLengthSquared(10.f * 10.f),
    m_maxDaughterLengthSquared(175.f * 175.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayIdentificationAlgorithm::Run()
{
    PfoVector parentPfos, daughterPfos;
    this->GetPfos(m_parentPfoListName, parentPfos);
    this->GetPfos(m_daughterPfoListName, daughterPfos);

    if (parentPfos.empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "DeltaRayIdentificationAlgorithm: pfo list " << m_parentPfoListName << " unavailable." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    // Build parent/daughter associations (currently using length and proximity)
    PfoAssociationMap pfoAssociationMap;
    this->BuildAssociationMap(parentPfos, daughterPfos, pfoAssociationMap);

    // Create the parent/daughter links
    PfoList newDaughterPfoList;
    this->BuildParentDaughterLinks(pfoAssociationMap, newDaughterPfoList);

    if (!newDaughterPfoList.empty())
    {
        PANDORA_THROW_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_parentPfoListName, m_daughterPfoListName, newDaughterPfoList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayIdentificationAlgorithm::GetPfos(const std::string &inputPfoListName, PfoVector &outputPfoVector) const
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, inputPfoListName, pPfoList));

    if (NULL == pPfoList)
        return;

    outputPfoVector.insert(outputPfoVector.end(), pPfoList->begin(), pPfoList->end());
    std::sort(outputPfoVector.begin(), outputPfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayIdentificationAlgorithm::BuildAssociationMap(const PfoVector &parentPfos, const PfoVector &daughterPfos, PfoAssociationMap &pfoAssociationMap) const
{
    PfoSet parentPfoList, daughterPfoList;
    parentPfoList.insert(parentPfos.begin(), parentPfos.end());
    daughterPfoList.insert(daughterPfos.begin(), daughterPfos.end());

    PfoVector allPfos(parentPfos.begin(), parentPfos.end());
    allPfos.insert(allPfos.end(), daughterPfos.begin(), daughterPfos.end());

    // Loop over possible daughter Pfos in primary list
    for (PfoVector::const_iterator iter1 = parentPfos.begin(), iterEnd1 = parentPfos.end(); iter1 != iterEnd1; ++iter1)
    {
        const ParticleFlowObject *const pDaughterPfo = *iter1;

        // Find the best parent Pfo using combined list
        const ParticleFlowObject *pBestParentPfo = NULL;
        float bestDisplacement(std::numeric_limits<float>::max());

        for (PfoVector::const_iterator iter2 = allPfos.begin(), iterEnd2 = allPfos.end(); iter2 != iterEnd2; ++iter2)
        {
            const ParticleFlowObject *const pThisParentPfo = *iter2;
            float thisDisplacement(std::numeric_limits<float>::max());

            if (pDaughterPfo == pThisParentPfo)
                continue;

            if (!this->IsAssociated(pDaughterPfo, pThisParentPfo, thisDisplacement))
                continue;

            if (thisDisplacement < bestDisplacement)
            {
                bestDisplacement = thisDisplacement;
                pBestParentPfo = pThisParentPfo;
            }
        }

        if (!pBestParentPfo)
            continue;

        // Case 1: candidate parent comes from primary list
        if (pBestParentPfo->GetParentPfoList().empty())
        {
            // Check: parent shouldn't live in the secondary list
            if (daughterPfoList.count(pBestParentPfo))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            pfoAssociationMap.insert(PfoAssociationMap::value_type(pDaughterPfo, pBestParentPfo));
        }

        // Case 2: candidate parent comes from secondary list
        else
        {
            // Check: parent shouldn't live in the primary list
            if (parentPfoList.count(pBestParentPfo))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            // Check: there should only be one parent
            if (pBestParentPfo->GetParentPfoList().size() != 1)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            // Check: get the new parent (and check there is no grand-parent)
            PfoList::const_iterator pIter = pBestParentPfo->GetParentPfoList().begin();
            const ParticleFlowObject *const pReplacementParentPfo = *pIter;
            if (pReplacementParentPfo->GetParentPfoList().size() != 0)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            pfoAssociationMap.insert(PfoAssociationMap::value_type(pDaughterPfo, pReplacementParentPfo));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayIdentificationAlgorithm::IsAssociated(
    const ParticleFlowObject *const pDaughterPfo, const ParticleFlowObject *const pParentPfo, float &displacement) const
{
    displacement = std::numeric_limits<float>::max();

    if (pDaughterPfo == pParentPfo)
        return false;

    const float daughterLengthSquared(LArPfoHelper::GetTwoDLengthSquared(pDaughterPfo));
    const float parentLengthSquared(LArPfoHelper::GetTwoDLengthSquared(pParentPfo));

    if (daughterLengthSquared > m_maxDaughterLengthSquared || parentLengthSquared < m_minParentLengthSquared || daughterLengthSquared > 0.5 * parentLengthSquared)
        return false;

    const float transitionLengthSquared(125.f);
    const float displacementCut((daughterLengthSquared > transitionLengthSquared)
            ? m_distanceForMatching
            : m_distanceForMatching * (2.f - daughterLengthSquared / transitionLengthSquared));

    try
    {
        displacement = this->GetTwoDSeparation(pDaughterPfo, pParentPfo);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    if (displacement > displacementCut)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayIdentificationAlgorithm::GetTwoDSeparation(const ParticleFlowObject *const pDaughterPfo, const ParticleFlowObject *const pParentPfo) const
{
    CartesianPointVector vertexVectorU, vertexVectorV, vertexVectorW;
    this->GetTwoDVertices(pDaughterPfo, TPC_VIEW_U, vertexVectorU);
    this->GetTwoDVertices(pDaughterPfo, TPC_VIEW_V, vertexVectorV);
    this->GetTwoDVertices(pDaughterPfo, TPC_VIEW_W, vertexVectorW);

    ClusterList clusterListU, clusterListV, clusterListW;
    LArPfoHelper::GetClusters(pParentPfo, TPC_VIEW_U, clusterListU);
    LArPfoHelper::GetClusters(pParentPfo, TPC_VIEW_V, clusterListV);
    LArPfoHelper::GetClusters(pParentPfo, TPC_VIEW_W, clusterListW);

    float sumViews(0.f);
    float sumDisplacementSquared(0.f);

    if (!vertexVectorU.empty())
    {
        const float thisDisplacement(this->GetClosestDistance(vertexVectorU, clusterListU));
        sumDisplacementSquared += thisDisplacement * thisDisplacement;
        sumViews += 1.f;
    }

    if (!vertexVectorV.empty())
    {
        const float thisDisplacement(this->GetClosestDistance(vertexVectorV, clusterListV));
        sumDisplacementSquared += thisDisplacement * thisDisplacement;
        sumViews += 1.f;
    }

    if (!vertexVectorW.empty())
    {
        const float thisDisplacement(this->GetClosestDistance(vertexVectorW, clusterListW));
        sumDisplacementSquared += thisDisplacement * thisDisplacement;
        sumViews += 1.f;
    }

    if (sumViews < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return std::sqrt(sumDisplacementSquared / sumViews);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayIdentificationAlgorithm::GetTwoDVertices(const ParticleFlowObject *const pPfo, const HitType &hitType, CartesianPointVector &vertexVector) const
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        CartesianVector firstCoordinate(0.f, 0.f, 0.f), secondCoordinate(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pCluster, firstCoordinate, secondCoordinate);

        vertexVector.push_back(firstCoordinate);
        vertexVector.push_back(secondCoordinate);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayIdentificationAlgorithm::GetClosestDistance(const CartesianPointVector &vertexVector, const ClusterList &clusterList) const
{
    if (vertexVector.empty() || clusterList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float bestDisplacement(std::numeric_limits<float>::max());

    for (CartesianPointVector::const_iterator iter1 = vertexVector.begin(), iterEnd1 = vertexVector.end(); iter1 != iterEnd1; ++iter1)
    {
        const CartesianVector &thisVertex = *iter1;

        for (ClusterList::const_iterator iter2 = clusterList.begin(), iterEnd2 = clusterList.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster *const pCluster = *iter2;
            const float thisDisplacement(LArClusterHelper::GetClosestDistance(thisVertex, pCluster));

            if (thisDisplacement < bestDisplacement)
                bestDisplacement = thisDisplacement;
        }
    }

    return bestDisplacement;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayIdentificationAlgorithm::BuildParentDaughterLinks(const PfoAssociationMap &pfoAssociationMap, PfoList &daughterPfoList) const
{
    PfoList pfoList;
    for (const auto &mapEntry : pfoAssociationMap)
        pfoList.push_back(mapEntry.first);
    pfoList.sort(LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pDaughterPfo : pfoList)
    {
        const ParticleFlowObject *const pParentPfo(this->GetParent(pfoAssociationMap, pDaughterPfo));

        if (!pParentPfo)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        if (!LArPfoHelper::IsTrack(pParentPfo))
            continue;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pParentPfo, pDaughterPfo));

        PandoraContentApi::ParticleFlowObject::Metadata metadata;
        metadata.m_particleId = E_MINUS;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pDaughterPfo, metadata));

        daughterPfoList.push_back(pDaughterPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *DeltaRayIdentificationAlgorithm::GetParent(const PfoAssociationMap &pfoAssociationMap, const ParticleFlowObject *const pPfo) const
{
    const ParticleFlowObject *pParentPfo = nullptr;
    const ParticleFlowObject *pDaughterPfo = pPfo;

    while (1)
    {
        PfoAssociationMap::const_iterator iter = pfoAssociationMap.find(pDaughterPfo);
        if (pfoAssociationMap.end() == iter)
            break;

        pParentPfo = iter->second;
        pDaughterPfo = pParentPfo;
    }

    return pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayIdentificationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParentPfoListName", m_parentPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DaughterPfoListName", m_daughterPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DistanceForMatching", m_distanceForMatching));

    float minParentLength = std::sqrt(m_minParentLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinParentLength", minParentLength));
    m_minParentLengthSquared = minParentLength * minParentLength;

    float maxDaughterLength = std::sqrt(m_maxDaughterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDaughterLength", maxDaughterLength));
    m_maxDaughterLengthSquared = maxDaughterLength * maxDaughterLength;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
