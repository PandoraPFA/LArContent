/**
 *  @file   larpandoracontent/LArStitching/MultiPandoraApi.h
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

typedef std::vector<int> VolumeIdList;
typedef std::unordered_map<const pandora::Pandora *, VolumeIdList> PandoraVolumeIdMap;

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
     *  @brief  Get the volume info block associated with a given pandora stitching instance and volume ID number
     *
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *  @param  idNumber the volume identifier number
     *
     *  @return the volume info block
     */
    static const VolumeInfo &GetVolumeInfo(const pandora::Pandora *const pPrimaryPandora, const int idNumber);

    /**
     *  @brief  Get the volume info block associated with a given pandora instance
     *
     *  @param  pPandora the address of the pandora instance
     *
     *  @return the volume info block
     */
    static const VolumeInfo &GetVolumeInfo(const pandora::Pandora *const pPandora);

    /**
     *  @brief  Get the list of volume IDs associated with a given primary pandora instance
     *
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *
     *  @return the volume ID list
     */
    static const VolumeIdList &GetVolumeIdList(const pandora::Pandora *const pPrimaryPandora);

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

private:
    static MultiPandoraApiImpl      m_multiPandoraApiImpl;          ///< The multi pandora api implementation
};

//------------------------------------------------------------------------------------------------------------------------------------------
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
     */
    VolumeInfo(const int idNumber);

    /**
     *  @brief  Get the volume identifier number
     *
     *  @return the volume identifier number
     */
    int GetIdNumber() const;

private:
    const int                       m_idNumber;                     ///< The volume identifier number
};

#endif // #ifndef MULTI_PANDORA_API_H
