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

const VolumeIdList &MultiPandoraApi::GetVolumeIdList()
{
    return m_multiPandoraApiImpl.GetVolumeIdList();
}

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

const VolumeInfo &MultiPandoraApi::GetVolumeInfo(const int volumeID)
{
    return m_multiPandoraApiImpl.GetVolumeInfo(volumeID);
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

void MultiPandoraApi::SetParticleX0(const pandora::Pandora *const pPandora, const pandora::ParticleFlowObject *const pPfo, const float x0)
{
    m_multiPandoraApiImpl.SetParticleX0(pPandora, pPfo, x0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MultiPandoraApi::ClearParticleX0Map(const pandora::Pandora *const pPandora)
{
    m_multiPandoraApiImpl.ClearParticleX0Map(pPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

VolumeInfo::VolumeInfo(const int idNumber, const std::string &idString, const float centerX, const float centerY, const float centerZ,
    const float widthX, const float widthY, const float widthZ, const bool isDriftInPositiveX) :
    m_idNumber(idNumber),
    m_idString(idString),
    m_centerX(centerX),
    m_centerY(centerY),
    m_centerZ(centerZ),
    m_widthX(widthX),
    m_widthY(widthY),
    m_widthZ(widthZ),
    m_isDriftInPositiveX(isDriftInPositiveX)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

int VolumeInfo::GetIdNumber() const
{
    return m_idNumber;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string &VolumeInfo::GetIdString() const
{
    return m_idString;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VolumeInfo::GetCenterX() const
{
    return m_centerX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VolumeInfo::GetCenterY() const
{
    return m_centerY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VolumeInfo::GetCenterZ() const
{
    return m_centerZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VolumeInfo::GetWidthX() const
{
    return m_widthX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VolumeInfo::GetWidthY() const
{
    return m_widthY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VolumeInfo::GetWidthZ() const
{
    return m_widthZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VolumeInfo::IsDriftInPositiveX() const
{
    return m_isDriftInPositiveX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VolumeInfo::GetParticleX0(const pandora::ParticleFlowObject *const pPfo) const
{
    ParticleX0Map::const_iterator iter = m_particleX0Map.find(pPfo);

    if (m_particleX0Map.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VolumeInfo::SetParticleX0(const pandora::ParticleFlowObject *const pPfo, const float x0)
{
    if (m_particleX0Map.count(pPfo))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!m_particleX0Map.insert(ParticleX0Map::value_type(pPfo, x0)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VolumeInfo::ClearParticleX0Map()
{
    m_particleX0Map.clear();
}
