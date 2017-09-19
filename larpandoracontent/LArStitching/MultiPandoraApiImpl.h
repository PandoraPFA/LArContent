/**
 *  @file   larpandoracontent/LArStitching/MultiPandoraApiImpl.h
 *
 *  @brief  Header file for the MultiPandoraApiImpl class.
 *
 *  $Log: $
 */
#ifndef MULTI_PANDORA_API_IMPL_H
#define MULTI_PANDORA_API_IMPL_H 1

#include "Objects/CartesianVector.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

#include <map>
#include <unordered_map>

/**
 *  @brief  MultiPandoraApiImpl class
 */
class MultiPandoraApiImpl
{
private:
    /**
     *  @brief  Default constructor;
     */
    MultiPandoraApiImpl();

    /**
     *  @brief  Destructor;
     */
    ~MultiPandoraApiImpl();

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
     *  @return the daughter pandora instance list
     */
    const PandoraInstanceList &GetDaughterPandoraInstanceList(const pandora::Pandora *const pPrimaryPandora) const;

    /**
     *  @brief  Get the address of the daughter pandora instance associated with a given primary pandora instance and volume id number
     *
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *  @param  idNumber the volume identifier number
     *
     *  @return the address of the daughter pandora instance
     */
    const pandora::Pandora *GetDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const int idNumber) const;

    /**
     *  @brief  Get the address of the primary pandora instance associated with a given daughter pandora instance
     *
     *  @param  pDaughterPandora the address of the daughter pandora instance
     *
     *  @return the address of the primary pandora instance
     */
    const pandora::Pandora *GetPrimaryPandoraInstance(const pandora::Pandora *const pDaughterPandora) const;

    /**
     *  @brief  Get the volume info block associated with a given pandora instance
     *
     *  @param  pPandora the address of the pandora instance
     *
     *  @return the volume info block
     */
    const VolumeInfo &GetVolumeInfo(const pandora::Pandora *const pPandora) const;

    /**
     *  @brief  Get the volume info block associated with a given pandora stitching instance and volume ID number
     *
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *  @param  idNumber the volume identifier number
     *
     *  @return the volume info block
     */
    const VolumeInfo &GetVolumeInfo(const pandora::Pandora *const pPrimaryPandora, const int idNumber) const;

    /**
     *  @brief  Get the list of volume IDs associated with a given primary pandora instance
     *
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *
     *  @return the volume ID list
     */
    const VolumeIdList &GetVolumeIdList(const pandora::Pandora *const pPrimaryPandora) const;

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

    PandoraInstanceMap              m_pandoraInstanceMap;   ///< The map from primary pandora instance to list of daughter pandora instances
    PandoraVolumeIdMap              m_pandoraVolumeIdMap;   ///< The map from primary pandora instance to list of volume IDs

    typedef std::unordered_map<const pandora::Pandora *, const pandora::Pandora *> PandoraRelationMap;
    PandoraRelationMap              m_pandoraRelationMap;   ///< The map from daughter pandora instance to primary pandora instance

    typedef std::unordered_map<const pandora::Pandora *, VolumeInfo *> VolumeInfoMap;
    VolumeInfoMap                   m_volumeInfoMap;        ///< The map from pandora instance to volume info object

    typedef std::map<int, VolumeInfo *> VolumeIdToInfoMap;
    typedef std::unordered_map<const pandora::Pandora *, VolumeIdToInfoMap> PrimaryToInfoMap;
    PrimaryToInfoMap                m_primaryToInfoMap; ///< The primary --> id --> volume info mapping

    typedef std::map<int, const pandora::Pandora *> VolumeIdToDaughterMap;
    typedef std::unordered_map<const pandora::Pandora *, VolumeIdToDaughterMap> PrimaryToDaughterMap;
    PrimaryToDaughterMap            m_primaryToDaughterMap; ///< The primary --> id --> daughter mapping

    friend class MultiPandoraApi;
};

#endif // #ifndef MULTI_PANDORA_API_IMPL_H
