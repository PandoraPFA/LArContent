/**
 *  @file   larpandoracontent/LArStitching/MultiPandoraApi.h
 *
 *  @brief  Header file for the MultiPandoraApi class.
 *
 *  $Log: $
 */
#ifndef MULTI_PANDORA_API_H
#define MULTI_PANDORA_API_H 1

#include <unordered_map>
#include <vector>

namespace pandora
{
    class Pandora;
    class ParticleFlowObject;
}

class MultiPandoraApiImpl;

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
     *  @param  volumeId the volume identifier number
     *
     *  @return the address of the daughter pandora instance
     */
    static const pandora::Pandora *GetDaughterPandoraInstance(const pandora::Pandora *const pPrimaryPandora, const unsigned int volumeId);

    /**
     *  @brief  Get the address of the primary pandora instance associated with a given daughter pandora instance
     *
     *  @param  pDaughterPandora the address of the daughter pandora instance
     *
     *  @return the address of the primary pandora instance
     */
    static const pandora::Pandora *GetPrimaryPandoraInstance(const pandora::Pandora *const pDaughterPandora);

    /**
     *  @brief  Get the volume id associated with a given pandora instance
     *
     *  @param  pPandora the address of the pandora instance
     *
     *  @return the volume id
     */
    static unsigned int GetVolumeId(const pandora::Pandora *const pPandora);

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
     *  @brief  Set the volume id associated with a given pandora instance
     *
     *  @param  pPandora the address of the pandora instance
     *  @param  volumeId the volume id
     */
    static void SetVolumeId(const pandora::Pandora *const pPandora, const unsigned int volumeId);

private:
    static MultiPandoraApiImpl      m_multiPandoraApiImpl;          ///< The multi pandora api implementation
};

#endif // #ifndef MULTI_PANDORA_API_H
