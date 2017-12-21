/**
 *  @file   larpandoracontent/LArControlFlow/MultiPandoraApi.cc
 *
 *  @brief  Implementation of the MultiPandoraApi class.
 *
 *  $Log: $
 */

#include "Pandora/Pandora.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApiImpl.h"

MultiPandoraApiImpl MultiPandoraApi::m_multiPandoraApiImpl;

//------------------------------------------------------------------------------------------------------------------------------------------

const PandoraInstanceMap &MultiPandoraApi::GetPandoraInstanceMap()
{
    return m_multiPandoraApiImpl.GetPandoraInstanceMap();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MultiPandoraApi::GetPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const unsigned int volumeId)
{
    return m_multiPandoraApiImpl.GetPandoraInstance(pPrimaryPandora, volumeId);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PandoraInstanceList &MultiPandoraApi::GetDaughterPandoraInstanceList(const pandora::Pandora *const pPrimaryPandora)
{
    return m_multiPandoraApiImpl.GetDaughterPandoraInstanceList(pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *MultiPandoraApi::GetPrimaryPandoraInstance(const pandora::Pandora *const pDaughterPandora)
{
    return m_multiPandoraApiImpl.GetPrimaryPandoraInstance(pDaughterPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int MultiPandoraApi::GetVolumeId(const pandora::Pandora *const pPandora)
{
    return m_multiPandoraApiImpl.GetVolumeId(pPandora);
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

void MultiPandoraApi::SetVolumeId(const pandora::Pandora *const pPandora, const unsigned int volumeId)
{
    m_multiPandoraApiImpl.SetVolumeId(pPandora, volumeId);
}
