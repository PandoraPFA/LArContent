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
     *  @brief  Get the address of the pandora instance associated with a given primary pandora instance and volume id number
     *
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *  @param  volumeId the volume identifier number
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *GetPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const unsigned int volumeId) const;

    /**
     *  @brief  Get the list of daughter pandora instances associated with a given primary pandora instance
     *
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *
     *  @return the daughter pandora instance list
     */
    const PandoraInstanceList &GetDaughterPandoraInstanceList(const pandora::Pandora *const pPrimaryPandora) const;

    /**
     *  @brief  Get the address of the primary pandora instance associated with a given daughter pandora instance
     *
     *  @param  pDaughterPandora the address of the daughter pandora instance
     *
     *  @return the address of the primary pandora instance
     */
    const pandora::Pandora *GetPrimaryPandoraInstance(const pandora::Pandora *const pDaughterPandora) const;

    /**
     *  @brief  Get the volume id associated with a given pandora instance
     *
     *  @param  pPandora the address of the pandora instance
     *
     *  @return the volume id
     */
    unsigned int GetVolumeId(const pandora::Pandora *const pPandora) const;

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
     *  @brief  Set the volume id associated with a given pandora instance
     *
     *  @param  pPandora the address of the pandora instance
     *  @param  volumeId the volume id
     */
    void SetVolumeId(const pandora::Pandora *const pPandora, const unsigned int volumeId);

    typedef std::unordered_map<const pandora::Pandora *, const pandora::Pandora *> PandoraRelationMap;
    typedef std::unordered_map<const pandora::Pandora *, unsigned int> PandoraToVolumeIdMap;

    PandoraInstanceMap              m_primaryToDaughtersMap;    ///< The map from primary pandora instance to list of daughter pandora instances
    PandoraRelationMap              m_daughterToPrimaryMap;     ///< The map from daughter pandora instance to primary pandora instance
    PandoraToVolumeIdMap            m_pandoraToVolumeIdMap;     ///< The map from pandora instance to volume id

    friend class MultiPandoraApi;
};

#endif // #ifndef MULTI_PANDORA_API_IMPL_H
