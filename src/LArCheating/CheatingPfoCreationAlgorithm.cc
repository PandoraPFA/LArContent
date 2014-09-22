/**
 *  @file   LArContent/src/LArCheating/CheatingPfoCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCheating/CheatingPfoCreationAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingPfoCreationAlgorithm::CheatingPfoCreationAlgorithm() :
    m_idOffsetU(0),
    m_idOffsetV(0),
    m_idOffsetW(0),
    m_useOnlyAvailableClusters(true),
    m_minGoodHitTypes(0),
    m_nHitsForGoodHitType(10)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoCreationAlgorithm::Run()
{
    const ClusterList *pInputClusterListU(NULL), *pInputClusterListV(NULL), *pInputClusterListW(NULL);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputClusterListNameU, pInputClusterListU));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputClusterListNameV, pInputClusterListV));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_inputClusterListNameW, pInputClusterListW));

    IdToClusterListMap idToClusterListMap;
    this->GetIdToClusterListMap(pInputClusterListU, m_idOffsetU, idToClusterListMap);
    this->GetIdToClusterListMap(pInputClusterListV, m_idOffsetV, idToClusterListMap);
    this->GetIdToClusterListMap(pInputClusterListW, m_idOffsetW, idToClusterListMap);

    IdToMCParticleMap idToMCParticleMap;
    this->GetIdToMCParticleMap(idToMCParticleMap);

    this->CreatePfos(idToClusterListMap, idToMCParticleMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPfoCreationAlgorithm::GetIdToClusterListMap(const ClusterList *const pClusterList, const int idOffset,
    IdToClusterListMap &idToClusterListMap) const
{
    if (NULL == pClusterList)
        return;

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            Cluster *pCluster(*iter);

            if (m_useOnlyAvailableClusters && !pCluster->IsAvailable())
                continue;

            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

            if (!m_particleIdList.empty() && !m_particleIdList.count(pMCParticle->GetParticleId()))
                continue;

            const int id(intptr_t(pMCParticle->GetUid()) - idOffset);
            idToClusterListMap[id].insert(pCluster);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPfoCreationAlgorithm::GetIdToMCParticleMap(IdToMCParticleMap &idToMCParticleMap) const
{
    const MCParticleList *pMCParticleList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticle3DListName, pMCParticleList));

    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        MCParticle *pMCParticle(*iter);
        const int id(intptr_t(pMCParticle->GetUid()));

        if (!idToMCParticleMap.insert(IdToMCParticleMap::value_type(id, pMCParticle)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPfoCreationAlgorithm::CreatePfos(const IdToClusterListMap &idToClusterListMap, const IdToMCParticleMap &idToMCParticleMap) const
{
    if (idToClusterListMap.empty())
        return;

    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (IdToClusterListMap::const_iterator iter = idToClusterListMap.begin(), iterEnd = idToClusterListMap.end(); iter != iterEnd; ++iter)
    {
        const int id(iter->first);
        const ClusterList &clusterList(iter->second);

        if (this->GetNHitTypesAboveThreshold(clusterList, m_nHitsForGoodHitType) < m_minGoodHitTypes)
            continue;

        IdToMCParticleMap::const_iterator mcIter = idToMCParticleMap.find(id);

        if (idToMCParticleMap.end() == mcIter)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        MCParticle *pMCParticle = mcIter->second;

        try
        {
            PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
            pfoParameters.m_particleId = pMCParticle->GetParticleId();
            pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
            pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
            pfoParameters.m_energy = pMCParticle->GetEnergy();
            pfoParameters.m_momentum = pMCParticle->GetMomentum();
            pfoParameters.m_clusterList.insert(clusterList.begin(), clusterList.end());

            ParticleFlowObject *pPfo(NULL);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
        }
        catch (StatusCodeException &)
        {
            std::cout << "CheatingPfoCreationAlgorithm: Could not create PFO for MCParticle with pdg code " << pMCParticle->GetParticleId() << std::endl;
        }
    }

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int CheatingPfoCreationAlgorithm::GetNHitTypesAboveThreshold(const ClusterList &clusterList, const unsigned int nHitsThreshold) const
{
    HitTypeMap hitTypeMap;

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        hitTypeMap[LArClusterHelper::GetClusterHitType(*cIter)] += (*cIter)->GetNCaloHits();
    }

    unsigned int nGoodViews(0);

    for (HitTypeMap::const_iterator hIter = hitTypeMap.begin(), hIterEnd = hitTypeMap.end(); hIter != hIterEnd; ++hIter)
    {
        if (hIter->second > nHitsThreshold)
            ++nGoodViews;
    }

    return nGoodViews;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "IdOffsetU", m_idOffsetU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "IdOffsetV", m_idOffsetV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "IdOffsetW", m_idOffsetW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticle3DListName", m_mcParticle3DListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseOnlyAvailableClusters", m_useOnlyAvailableClusters));

    IntVector particleIdVector;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ParticleIdList", particleIdVector));

    m_particleIdList.insert(particleIdVector.begin(), particleIdVector.end());

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinGoodHitTypes", m_minGoodHitTypes));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NHitsForGoodHitType", m_nHitsForGoodHitType));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
