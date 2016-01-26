/**
 *  @file   LArContent/include/LArStitching/MultiPandora.h
 * 
 *  @brief  Header file for the MultiPandora class.
 * 
 *  $Log: $
 */
#ifndef MULTI_PANDORA_H
#define MULTI_PANDORA_H 1

#include "Objects/CartesianVector.h"

#include "Pandora/StatusCodes.h"

#include <string>
#include <unordered_map>

class MultiPandora;
class VolumeInfo;

namespace pandora { class Pandora; }

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MultiPandoraApi class
 */
class MultiPandoraApi
{
public:
    /**
     *  @brief  Get the multi pandora instance
     * 
     *  @return the multipandora instance
     */
    static const MultiPandora &GetMultiPandora();

    /**
     *  @brief  Declare a new primary pandora instance and receive the relevant multi pandora book-keeping instance
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     * 
     *  @return the multipandora instance
     */
    static void AddPrimaryPandoraInstance(const pandora::Pandora *const pPrimaryPandora);

    /**
     *  @brief  Add a pandora daughter instance, associated to a primary pandora instance
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *  @param  pDaughterPandora the address of the daughter pandora instance
     */
    static void AddDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const pandora::Pandora *const pDaughterPandora);

    /**
     *  @brief  Delete all pandora instances associated with (and including) a specified primary pandora instance
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     */
    static void DeletePandoraInstances(const pandora::Pandora *const pPrimaryPandora);

    /**
     *  @brief  Set the volume info block associated with a given pandora instance
     * 
     *  @param  pPandora the address of the pandora instance
     *  @param  pVolumeInfo the address of the volume info block
     */
    static void SetVolumeInfo(const pandora::Pandora *const pPandora, VolumeInfo *const pVolumeInfo);

    /**
     *  @brief  Set the x0 value for a specified particle
     * 
     *  @param  pPandora the address of the pandora instance
     *  @param  pPfo the address of the particle
     *  @param  x0 the x0 value for the particle
     */
    static void SetParticleX0(const pandora::Pandora *const pPandora, const pandora::ParticleFlowObject *const pPfo, const float x0);

    /**
     *  @brief  Clear the particle x0 map
     * 
     *  @param  pPandora the address of the pandora instance
     */
    static void ClearParticleX0Map(const pandora::Pandora *const pPandora);

private:
    static MultiPandora             m_multiPandora;             ///< The multi pandora
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MultiPandora class
 */
class MultiPandora
{
public:
    typedef std::vector<const pandora::Pandora *> PandoraInstanceList;
    typedef std::unordered_map<const pandora::Pandora *, PandoraInstanceList> PandoraInstanceMap;

    /**
     *  @brief  Get the pandora instance map
     * 
     *  @return the pandora instance map
     */
    const PandoraInstanceMap &GetPandoraInstanceMap() const;

    /**
     *  @brief  Get the list of daughter pandora instances associated with a given primary pandora instance
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     * 
     *  @return the pandora instance list
     */
    const PandoraInstanceList &GetDaughterPandoraInstances(const pandora::Pandora *const pPrimaryPandora) const;

    /**
     *  @brief  Get the volume info block associated with a given pandora instance
     * 
     *  @param  pPandora the address of the pandora instance
     * 
     *  @return the volume info block
     */
    const VolumeInfo &GetVolumeInfo(const pandora::Pandora *const pPandora) const;

private:
    /**
     *  @brief  Default constructor;
     */
    MultiPandora();

    /**
     *  @brief  Declare a new primary pandora instance
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     */
    void AddPrimaryPandoraInstance(const pandora::Pandora *const pPrimaryPandora);

    /**
     *  @brief  Add a pandora daughter instance, associated to a primary pandora instance
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *  @param  pDaughterPandora the address of the daughter pandora instance
     */
    void AddDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const pandora::Pandora *const pDaughterPandora);

    /**
     *  @brief  Delete all pandora instances associated with (and including) a specified primary pandora instance
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     */
    void DeletePandoraInstances(const pandora::Pandora *const pPrimaryPandora);

    /**
     *  @brief  Set the volume info block associated with a given pandora instance
     * 
     *  @param  pPandora the address of the pandora instance
     *  @param  pVolumeInfo the address of the volume info block
     */
    void SetVolumeInfo(const pandora::Pandora *const pPandora, VolumeInfo *const pVolumeInfo);

    /**
     *  @brief  Set the x0 value for a specified particle
     * 
     *  @param  pPandora the address of the pandora instance
     *  @param  pPfo the address of the particle
     *  @param  x0 the x0 value for the particle
     */
    void SetParticleX0(const pandora::Pandora *const pPandora, const pandora::ParticleFlowObject *const pPfo, const float x0);

    /**
     *  @brief  Clear the particle x0 map
     * 
     *  @param  pPandora the address of the pandora instance
     */
    void ClearParticleX0Map(const pandora::Pandora *const pPandora);

    typedef std::unordered_map<const pandora::Pandora *, VolumeInfo*> VolumeInfoMap;

    VolumeInfoMap                   m_volumeInfoMap;        ///< The volume info map
    PandoraInstanceMap              m_pandoraInstanceMap;   ///< The pandora instance map

    friend class MultiPandoraApi;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  VolumeInfo class
 */
class VolumeInfo
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  identifier the volume identifier
     *  @param  center the center of the drift volume 3D coordinate system, in the World volume
     *  @param  isDriftInPositiveX whether the drift direction for electrons corresponds to positive x direction
     */
    VolumeInfo(const std::string &identifier, const pandora::CartesianVector &center, const bool isDriftInPositiveX);

    /**
     *  @brief  Get the volume identifier
     * 
     *  @return the volume identifier
     */
    const std::string &GetIdentifier() const;

    /**
     *  @brief  Get the center of the drift volume 3D coordinate system, in the World volume
     * 
     *  @return the center of the drift volume 3D coordinate system, in the World volume
     */
    const pandora::CartesianVector &GetCenter() const;

    /**
     *  @brief  Whether the drift direction for electrons corresponds to positive x direction
     * 
     *  @return boolean
     */
    bool IsDriftInPositiveX() const;

    /**
     *  @brief  Get the x0 value for a specified particle
     * 
     *  @param  pPfo the address of the particle
     * 
     *  @return the x0 value for the particle
     */
    float GetParticleX0(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Set the x0 value for a specified particle
     * 
     *  @param  pPfo the address of the particle
     *  @param  x0 the x0 value for the particle
     */
    void SetParticleX0(const pandora::ParticleFlowObject *const pPfo, const float x0);

    /**
     *  @brief  Clear the particle x0 map
     */
    void ClearParticleX0Map();

private:
    typedef std::map<const pandora::ParticleFlowObject*, float> ParticleX0Map;

    std::string                 m_identifier;           ///< The volume identifier
    pandora::CartesianVector    m_center;               ///< The center of the drift volume 3D coordinate system, in the World volume
    bool                        m_isDriftInPositiveX;   ///< Whether the drift direction for electrons corresponds to positive x direction
    ParticleX0Map               m_particleX0Map;        ///< The particle x0 map (x0 is t0 times drift velocity)
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const MultiPandora &MultiPandoraApi::GetMultiPandora()
{
    return m_multiPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandoraApi::AddPrimaryPandoraInstance(const pandora::Pandora *const pPrimaryPandora)
{
    m_multiPandora.AddPrimaryPandoraInstance(pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandoraApi::AddDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const pandora::Pandora *const pDaughterPandora)
{
    m_multiPandora.AddDaughterPandoraInstance(pPrimaryPandora, pDaughterPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandoraApi::DeletePandoraInstances(const pandora::Pandora *const pPrimaryPandora)
{
    if (!m_multiPandora.GetPandoraInstanceMap().count(pPrimaryPandora))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    m_multiPandora.DeletePandoraInstances(pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandoraApi::SetVolumeInfo(const pandora::Pandora *const pPandora, VolumeInfo *const pVolumeInfo)
{
    m_multiPandora.SetVolumeInfo(pPandora, pVolumeInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandoraApi::SetParticleX0(const pandora::Pandora *const pPandora, const pandora::ParticleFlowObject *const pPfo, const float x0)
{
    m_multiPandora.SetParticleX0(pPandora, pPfo, x0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandoraApi::ClearParticleX0Map(const pandora::Pandora *const pPandora)
{
    m_multiPandora.ClearParticleX0Map(pPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const MultiPandora::PandoraInstanceMap &MultiPandora::GetPandoraInstanceMap() const
{
    return m_pandoraInstanceMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const MultiPandora::PandoraInstanceList &MultiPandora::GetDaughterPandoraInstances(const pandora::Pandora *const pPrimaryPandora) const
{
    PandoraInstanceMap::const_iterator iter = m_pandoraInstanceMap.find(pPrimaryPandora);

    if (m_pandoraInstanceMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const VolumeInfo &MultiPandora::GetVolumeInfo(const pandora::Pandora *const pPandora) const
{
    VolumeInfoMap::const_iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return *(iter->second);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandora::SetVolumeInfo(const pandora::Pandora *const pPandora, VolumeInfo *const pVolumeInfo)
{
    if (m_volumeInfoMap.count(pPandora))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!m_volumeInfoMap.insert(VolumeInfoMap::value_type(pPandora, pVolumeInfo)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandora::SetParticleX0(const pandora::Pandora *const pPandora, const pandora::ParticleFlowObject *const pPfo, const float x0)
{
    VolumeInfoMap::iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second->SetParticleX0(pPfo, x0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandora::ClearParticleX0Map(const pandora::Pandora *const pPandora)
{
    VolumeInfoMap::iterator iter = m_volumeInfoMap.find(pPandora);

    if (m_volumeInfoMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second->ClearParticleX0Map();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline MultiPandora::MultiPandora()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandora::AddPrimaryPandoraInstance(const pandora::Pandora *const pPrimaryPandora)
{
    if (m_pandoraInstanceMap.count(pPrimaryPandora))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!m_pandoraInstanceMap.insert(PandoraInstanceMap::value_type(pPrimaryPandora, PandoraInstanceList())).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandora::AddDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const pandora::Pandora *const pDaughterPandora)
{
    PandoraInstanceMap::iterator iter = m_pandoraInstanceMap.find(pPrimaryPandora);

    if (m_pandoraInstanceMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    iter->second.push_back(pDaughterPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MultiPandora::DeletePandoraInstances(const pandora::Pandora *const pPrimaryPandora)
{
    PandoraInstanceMap::const_iterator iter = m_pandoraInstanceMap.find(pPrimaryPandora);

    if (m_pandoraInstanceMap.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    PandoraInstanceList pandoraInstanceList(iter->second);
    pandoraInstanceList.push_back(pPrimaryPandora);
    m_pandoraInstanceMap.erase(iter);
    
    for (const pandora::Pandora *const pPandora : pandoraInstanceList)
    {
        VolumeInfoMap::iterator volumeInfoIter = m_volumeInfoMap.find(pPandora);

        if (m_volumeInfoMap.end() != volumeInfoIter)
        {
            delete volumeInfoIter->second;
            m_volumeInfoMap.erase(volumeInfoIter);
        }
        
        delete pPandora;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline VolumeInfo::VolumeInfo(const std::string &identifier, const pandora::CartesianVector &center, const bool isDriftInPositiveX) :
    m_identifier(identifier),
    m_center(center),
    m_isDriftInPositiveX(isDriftInPositiveX)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &VolumeInfo::GetIdentifier() const
{
    return m_identifier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &VolumeInfo::GetCenter() const
{
    return m_center;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool VolumeInfo::IsDriftInPositiveX() const
{
    return m_isDriftInPositiveX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VolumeInfo::GetParticleX0(const pandora::ParticleFlowObject *const pPfo) const
{
    ParticleX0Map::const_iterator iter = m_particleX0Map.find(pPfo);

    if (m_particleX0Map.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VolumeInfo::SetParticleX0(const pandora::ParticleFlowObject *const pPfo, const float x0)
{
    if (m_particleX0Map.count(pPfo))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!m_particleX0Map.insert(ParticleX0Map::value_type(pPfo, x0)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VolumeInfo::ClearParticleX0Map()
{
    m_particleX0Map.clear();
}

#endif // #ifndef MULTI_PANDORA_H
