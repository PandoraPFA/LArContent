/**
 *  @file   larpandoracontent/LArStitching/MultiPandoraApi.cc
 *
 *  @brief  Implementation of the MultiPandoraApi class.
 *
 *  $Log: $
 */

#include "Pandora/Pandora.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"
#include "larpandoracontent/LArStitching/MultiPandoraApiImpl.h"

MultiPandoraApiImpl MultiPandoraApi::m_multiPandoraApiImpl;

//------------------------------------------------------------------------------------------------------------------------------------------

const PandoraInstanceMap &MultiPandoraApi::GetPandoraInstanceMap()
{
    return m_multiPandoraApiImpl.GetPandoraInstanceMap();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PandoraInstanceList &MultiPandoraApi::GetDaughterPandoraInstanceList(const pandora::Pandora *const pPrimaryPandora)
{
    return m_multiPandoraApiImpl.GetDaughterPandoraInstanceList(pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MultiPandoraApi::GetDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const int idNumber)
{
    return m_multiPandoraApiImpl.GetDaughterPandoraInstance(pPrimaryPandora, idNumber);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MultiPandoraApi::GetPrimaryPandoraInstance(const pandora::Pandora *const pDaughterPandora)
{
    return m_multiPandoraApiImpl.GetPrimaryPandoraInstance(pDaughterPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VolumeInfo &MultiPandoraApi::GetVolumeInfo(const pandora::Pandora *const pPandora)
{
    return m_multiPandoraApiImpl.GetVolumeInfo(pPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VolumeInfo &MultiPandoraApi::GetVolumeInfo(const pandora::Pandora *const pPrimaryPandora, const int idNumber)
{
    return m_multiPandoraApiImpl.GetVolumeInfo(pPrimaryPandora, idNumber);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VolumeIdList &MultiPandoraApi::GetVolumeIdList(const pandora::Pandora *const pPrimaryPandora)
{
    return m_multiPandoraApiImpl.GetVolumeIdList(pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApi::AddPrimaryPandoraInstance(const pandora::Pandora *const pPrimaryPandora)
{
    m_multiPandoraApiImpl.AddPrimaryPandoraInstance(pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApi::AddDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const pandora::Pandora *const pDaughterPandora)
{
    m_multiPandoraApiImpl.AddDaughterPandoraInstance(pPrimaryPandora, pDaughterPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApi::DeletePandoraInstances(const pandora::Pandora *const pPrimaryPandora)
{
    m_multiPandoraApiImpl.DeletePandoraInstances(pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApi::SetVolumeInfo(const pandora::Pandora *const pPandora, VolumeInfo *const pVolumeInfo)
{
    m_multiPandoraApiImpl.SetVolumeInfo(pPandora, pVolumeInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

VolumeInfo::VolumeInfo(const int idNumber) :
    m_idNumber(idNumber)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

int VolumeInfo::GetIdNumber() const
{
    return m_idNumber;
}
