/**
 *  @file   larpandoracontent/LArContent.h
 *
 *  @brief  Header file detailing content for use with particle flow reconstruction at liquid argon time projection chambers
 *
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

namespace pandora
{
class Pandora;
}

/**
 *  @brief  LArContent class
 */
class LArContent
{
public:
    /**
     *  @brief  Register all the lar content algorithms and tools with pandora
     *
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);

    /**
     *  @brief  Register the basic lar content plugins with pandora
     *
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterBasicPlugins(const pandora::Pandora &pandora);
};

#endif // #ifndef LAR_CONTENT_H
