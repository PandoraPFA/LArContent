/**
 *  @file   larpandoracontent/LArStitching/MultiPandoraApiImpl.cc
 *
 *  @brief  Implementation of the MultiPandoraApiImpl class.
 *
 *  $Log: $
 */

#include "Pandora/Pandora.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArStitching/MultiPandoraApiImpl.h"

const PandoraInstanceMap &MultiPandoraApiImpl::GetPandoraInstanceMap() const
{
    return m_pandoraInstanceMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PandoraInstanceList &MultiPandoraApiImpl::GetDaughterPandoraInstanceList(const pandora::Pandora *const pPrimaryPandora) const
{
    PandoraInstanceMap::const_iterator iter = m_pandoraInstanceMap.find(pPrimaryPandora);

    if (m_pandoraInstanceMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MultiPandoraApiImpl::GetDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const int idNumber) const
{
    PrimaryToDaughterMap::const_iterator iter = m_primaryToDaughterMap.find(pPrimaryPandora);

    if (m_primaryToDaughterMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    VolumeIdToDaughterMap::const_iterator idIter = iter->second.find(idNumber);

    if (iter->second.end() == idIter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return idIter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MultiPandoraApiImpl::GetPrimaryPandoraInstance(const pandora::Pandora *const pDaughterPandora) const
{
    PandoraRelationMap::const_iterator iter = m_pandoraRelationMap.find(pDaughterPandora);

    if (m_pandoraRelationMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VolumeInfo &MultiPandoraApiImpl::GetVolumeInfo(const pandora::Pandora *const pPandora) const
{
    VolumeInfoMap::const_iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return *(iter->second);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VolumeInfo &MultiPandoraApiImpl::GetVolumeInfo(const pandora::Pandora *const pPrimaryPandora, const int idNumber) const
{
    PrimaryToInfoMap::const_iterator iter = m_primaryToInfoMap.find(pPrimaryPandora);

    if (m_primaryToInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    VolumeIdToInfoMap::const_iterator idIter = iter->second.find(idNumber);

    if (iter->second.end() == idIter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return *(idIter->second);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VolumeIdList &MultiPandoraApiImpl::GetVolumeIdList(const pandora::Pandora *const pPrimaryPandora) const
{
    PandoraVolumeIdMap::const_iterator iter = m_pandoraVolumeIdMap.find(pPrimaryPandora);

    if (m_pandoraVolumeIdMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::SetVolumeInfo(const pandora::Pandora *const pPandora, VolumeInfo *const pVolumeInfo)
{
    if (m_volumeInfoMap.count(pPandora))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!m_volumeInfoMap.insert(VolumeInfoMap::value_type(pPandora, pVolumeInfo)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const pandora::Pandora *const pPrimaryPandora((m_pandoraInstanceMap.count(pPandora)) ? pPandora : this->GetPrimaryPandoraInstance(pPandora));

    // Populate primary --> id ---> daughter map
    PrimaryToDaughterMap::iterator iter1 = m_primaryToDaughterMap.find(pPrimaryPandora);

    if (m_primaryToDaughterMap.end() == iter1)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    if (!iter1->second.insert(VolumeIdToDaughterMap::value_type(pVolumeInfo->GetIdNumber(), pPandora)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    // Populate primary --> id --> info map
    PrimaryToInfoMap::iterator iter2 = m_primaryToInfoMap.find(pPrimaryPandora);

    if (m_primaryToInfoMap.end() == iter2)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    if (!iter2->second.insert(VolumeIdToInfoMap::value_type(pVolumeInfo->GetIdNumber(), pVolumeInfo)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    // Populate primary -> id list
    PandoraVolumeIdMap::iterator iter3 = m_pandoraVolumeIdMap.find(pPrimaryPandora);

    if (m_pandoraVolumeIdMap.end() == iter3)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter3->second.push_back(pVolumeInfo->GetIdNumber());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::SetParticleX0(const pandora::Pandora *const pPandora, const pandora::ParticleFlowObject *const pPfo, const float x0)
{
    VolumeInfoMap::iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second->SetParticleX0(pPfo, x0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::ClearParticleX0Map(const pandora::Pandora *const pPandora)
{
    VolumeInfoMap::iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second->ClearParticleX0Map();
}

//------------------------------------------------------------------------------------------------------------------------------------------

MultiPandoraApiImpl::MultiPandoraApiImpl()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MultiPandoraApiImpl::~MultiPandoraApiImpl()
{
    PandoraInstanceMap pandoraInstanceMap(m_pandoraInstanceMap);

    for (const auto &mapElement : pandoraInstanceMap)
        this->DeletePandoraInstances(mapElement.first);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::AddPrimaryPandoraInstance(const pandora::Pandora *const pPrimaryPandora)
{
    if (m_pandoraInstanceMap.count(pPrimaryPandora))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!m_pandoraInstanceMap.insert(PandoraInstanceMap::value_type(pPrimaryPandora, PandoraInstanceList())).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    if (!m_pandoraVolumeIdMap.insert(PandoraVolumeIdMap::value_type(pPrimaryPandora, VolumeIdList())).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    if (!m_primaryToDaughterMap.insert(PrimaryToDaughterMap::value_type(pPrimaryPandora, VolumeIdToDaughterMap())).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    if (!m_primaryToInfoMap.insert(PrimaryToInfoMap::value_type(pPrimaryPandora, VolumeIdToInfoMap())).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::AddDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const pandora::Pandora *const pDaughterPandora)
{
    PandoraInstanceMap::iterator iter = m_pandoraInstanceMap.find(pPrimaryPandora);

    if (m_pandoraInstanceMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second.push_back(pDaughterPandora);

    if (!m_pandoraRelationMap.insert(PandoraRelationMap::value_type(pDaughterPandora, pPrimaryPandora)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApiImpl::DeletePandoraInstances(const pandora::Pandora *const pPrimaryPandora)
{
    PandoraInstanceList pandoraInstanceList;

    try
    {
        pandoraInstanceList = this->GetDaughterPandoraInstanceList(pPrimaryPandora);
    }
    catch (const pandora::StatusCodeException &)
    {
        std::cout << "MultiPandoraApiImpl::DeletePandoraInstances - unable to find daughter instances associated with primary " << pPrimaryPandora << std::endl;
    }

    pandoraInstanceList.push_back(pPrimaryPandora);

    m_pandoraInstanceMap.erase(pPrimaryPandora);
    m_pandoraVolumeIdMap.erase(pPrimaryPandora);
    m_primaryToDaughterMap.erase(pPrimaryPandora);
    m_primaryToInfoMap.erase(pPrimaryPandora);

    for (const pandora::Pandora *const pPandora : pandoraInstanceList)
    {
        VolumeInfoMap::iterator volumeInfoIter = m_volumeInfoMap.find(pPandora);

        if (m_volumeInfoMap.end() != volumeInfoIter)
        {
            delete volumeInfoIter->second;
            m_volumeInfoMap.erase(volumeInfoIter);
        }

        m_pandoraRelationMap.erase(pPandora);
        delete pPandora;
    }
}
