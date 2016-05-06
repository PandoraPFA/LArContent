/**
 *  @file   LArContent/LArStitching/MultiPandoraApi.h
 * 
 *  @brief  Header file for the MultiPandoraApi class.
 * 
 *  $Log: $
 */
#ifndef MULTI_PANDORA_API_H
#define MULTI_PANDORA_API_H 1

#include "Objects/CartesianVector.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace pandora
{
    class Pandora;
    class ParticleFlowObject;
}

class MultiPandoraApiImpl;
class VolumeInfo;

typedef std::vector<const pandora::Pandora *> PandoraInstanceList;
typedef std::unordered_map<const pandora::Pandora *, PandoraInstanceList> PandoraInstanceMap;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MultiPandoraApi class
 */
class MultiPandoraApi
{
public:
    /**
     *  @brief  Get the pandora instance map
     * 
     *  @return the pandora instance map
     */
    static const PandoraInstanceMap &GetPandoraInstanceMap();

    /**
     *  @brief  Get the list of daughter pandora instances associated with a given primary pandora instance
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     * 
     *  @return the daughter pandora instance list
     */
    static const PandoraInstanceList &GetDaughterPandoraInstanceList(const pandora::Pandora *const pPrimaryPandora);

    /**
     *  @brief  Get the address of the daughter pandora instance associated with a given primary pandora instance and volume id number
     * 
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *  @param  idNumber the volume identifier number
     * 
     *  @return the address of the daughter pandora instance
     */
    static const pandora::Pandora *GetDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const int idNumber);

    /**
     *  @brief  Get the address of the primary pandora instance associated with a given daughter pandora instance
     * 
     *  @param  pDaughterPandora the address of the daughter pandora instance
     * 
     *  @return the address of the primary pandora instance
     */
    static const pandora::Pandora *GetPrimaryPandoraInstance(const pandora::Pandora *const pDaughterPandora);

    /**
     *  @brief  Get the volume info block associated with a given pandora instance
     * 
     *  @param  pPandora the address of the pandora instance
     * 
     *  @return the volume info block
     */
    static const VolumeInfo &GetVolumeInfo(const pandora::Pandora *const pPandora);

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
    static MultiPandoraApiImpl      m_multiPandoraApiImpl;          ///< The multi pandora api implementation
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
     *  @param  idNumber the volume identifier number
     *  @param  idString the volume identifier string or name
     *  @param  center the center of the drift volume 3D coordinate system, in the World volume
     *  @param  isDriftInPositiveX whether the drift direction for electrons corresponds to positive x direction
     */
    VolumeInfo(const int idNumber, const std::string &idString, const pandora::CartesianVector &center, const bool isDriftInPositiveX);

    /**
     *  @brief  Get the volume identifier number
     * 
     *  @return the volume identifier number
     */
    int GetIdNumber() const;

    /**
     *  @brief  Get the volume identifier string or name
     * 
     *  @return the volume identifier string or name
     */
    const std::string &GetIdString() const;

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
    typedef std::unordered_map<const pandora::ParticleFlowObject*, float> ParticleX0Map;

    const int                       m_idNumber;             ///< The volume identifier number
    const std::string               m_idString;             ///< The volume identifier string or name
    const pandora::CartesianVector  m_center;               ///< The center of the drift volume 3D coordinate system, in the World volume
    const bool                      m_isDriftInPositiveX;   ///< Whether the drift direction for electrons corresponds to positive x direction
    ParticleX0Map                   m_particleX0Map;        ///< The particle x0 map (x0 is t0 times drift velocity)
};

#endif // #ifndef MULTI_PANDORA_API_H
