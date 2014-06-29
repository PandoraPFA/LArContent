/**
 *  @file   LArContent/src/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray identification algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode DeltaRayIdentificationAlgorithm::Run()
{
    PfoList inputPfos, outputPfos;
    this->GetPfos(m_inputPfoListName, inputPfos);
    this->GetPfos(m_outputPfoListName, outputPfos);

    if (inputPfos.empty())
    {
        std::cout << "DeltaRayIdentificationAlgorithm: could not find pfo list " << m_inputPfoListName << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    // Build parent/daughter associations (currently using length and proximity)
    PfoAssociationMap pfoAssociationMap;
    this->BuildAssociationMap(inputPfos, outputPfos, pfoAssociationMap);

    // Create the parent/daughter links
    PfoList daughterPfoList;
    this->BuildParentDaughterLinks(pfoAssociationMap, daughterPfoList);

    if (!daughterPfoList.empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_inputPfoListName, m_outputPfoListName,
            daughterPfoList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayIdentificationAlgorithm::GetPfos(const std::string inputPfoListName, PfoList &outputPfoList) const
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputPfoListName, pPfoList));

    if (NULL == pPfoList)
        return;

    outputPfoList.insert(pPfoList->begin(), pPfoList->end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayIdentificationAlgorithm::BuildAssociationMap(const PfoList &inputPfos, const PfoList &outputPfos,
    PfoAssociationMap &pfoAssociationMap) const
{
    PfoList allPfos;
    allPfos.insert(inputPfos.begin(), inputPfos.end());
    allPfos.insert(outputPfos.begin(), outputPfos.end());

    // Loop over possible daughter Pfos in primary list
    for (PfoList::const_iterator iter1 = inputPfos.begin(), iterEnd1 = inputPfos.end(); iter1 != iterEnd1; ++iter1)
    {
        const ParticleFlowObject *pDaughterPfo = *iter1;

        // Find the best parent Pfo using combined list
        ParticleFlowObject *pBestParentPfo = NULL;
        float bestDisplacement(std::numeric_limits<float>::max());

        for (PfoList::const_iterator iter2 = allPfos.begin(), iterEnd2 = allPfos.end(); iter2 != iterEnd2; ++iter2)
        {
            ParticleFlowObject *pThisParentPfo = *iter2;
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
            if (outputPfos.count(pBestParentPfo))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            pfoAssociationMap.insert(PfoAssociationMap::value_type(pDaughterPfo, pBestParentPfo));
        }

        // Case 2: candidate parent comes from secondary list
        else
        {
            // Check: parent shouldn't live in the primary list
            if (inputPfos.count(pBestParentPfo))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            // Check: there should only be one parent
            if (pBestParentPfo->GetParentPfoList().size() != 1)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            // Check: get the new parent (and check there is no grand-parent)
            PfoList::iterator pIter = pBestParentPfo->GetParentPfoList().begin();
            ParticleFlowObject *pReplacementParentPfo = *pIter;
            if (pReplacementParentPfo->GetParentPfoList().size() != 0)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            pfoAssociationMap.insert(PfoAssociationMap::value_type(pDaughterPfo, pReplacementParentPfo));
        }

    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayIdentificationAlgorithm::IsAssociated(const ParticleFlowObject *const pDaughterPfo, const ParticleFlowObject *const pParentPfo,
    float &displacementSquared) const
{
    displacementSquared = std::numeric_limits<float>::max();

    if (pDaughterPfo == pParentPfo)
        return false;

    float sumViews(0.f);
    float sumDeltaSquared(0.f);

    const ClusterList &daughterList(pDaughterPfo->GetClusterList());
    const ClusterList &parentList(pParentPfo->GetClusterList());

    for (ClusterList::const_iterator cIter1 = daughterList.begin(), cIterEnd1 = daughterList.end(); cIter1 != cIterEnd1; ++cIter1)
    {
        const Cluster *pDaughterCluster = *cIter1;

        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        if (TPC_3D == daughterHitType)
            continue;

        bool foundMatch(false);

        for (ClusterList::const_iterator cIter2 = parentList.begin(), cIterEnd2 = parentList.end(); cIter2 != cIterEnd2; ++cIter2)
        {
            const Cluster *pParentCluster = *cIter2;

            if (daughterHitType != LArClusterHelper::GetClusterHitType(pParentCluster))
                continue;

            if (LArClusterHelper::GetLengthSquared(pDaughterCluster) > 0.5f * LArClusterHelper::GetLengthSquared(pParentCluster))
                break;

            const float delta(LArClusterHelper::GetClosestDistance(pParentCluster, pDaughterCluster));
            const float deltaSquared(delta * delta);

            if (deltaSquared > 5.f * m_maxDisplacementSquared)
                break;

            sumDeltaSquared += deltaSquared;
            sumViews += 1.f;
            foundMatch = true;
        }

        if (!foundMatch)
            return false;
    }

    if (sumViews < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (sumDeltaSquared / sumViews > m_maxDisplacementSquared)
        return false;

    displacementSquared = sumDeltaSquared / sumViews;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayIdentificationAlgorithm::BuildParentDaughterLinks(const PfoAssociationMap &pfoAssociationMap, PfoList &outputPfoList) const
{
    for (PfoAssociationMap::const_iterator iter = pfoAssociationMap.begin(), iterEnd = pfoAssociationMap.end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *pPfo = iter->first;

        ParticleFlowObject *pDaughterPfo = const_cast<ParticleFlowObject*>(pPfo);
        ParticleFlowObject *pParentPfo(this->GetParent(pfoAssociationMap, pDaughterPfo));

        if (NULL == pParentPfo)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pParentPfo, pDaughterPfo));
        outputPfoList.insert(pDaughterPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleFlowObject *DeltaRayIdentificationAlgorithm::GetParent(const PfoAssociationMap &pfoAssociationMap,
    const ParticleFlowObject *const pPfo) const
{
    ParticleFlowObject *pParentPfo = NULL;
    ParticleFlowObject *pDaughterPfo = const_cast<ParticleFlowObject*>(pPfo);

    while(1)
    {
        PfoAssociationMap::const_iterator iter = pfoAssociationMap.find(pDaughterPfo);
        if (pfoAssociationMap.end() == iter)
            break;

        pParentPfo = const_cast<ParticleFlowObject*>(iter->second);
        pDaughterPfo = pParentPfo;
    }

    return pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayIdentificationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    float maxDisplacement = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDisplacement", maxDisplacement));
    m_maxDisplacementSquared = maxDisplacement * maxDisplacement;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
