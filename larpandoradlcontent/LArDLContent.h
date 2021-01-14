/**
 *  @file   larpandoradlcontent/LArDLContent.h
 *
 *  @brief  Header file detailing content for use with particle flow reconstruction at liquid argon time projection chambers
 *
 *  $Log: $
 */
#ifndef LAR_DL_CONTENT_H
#define LAR_DL_CONTENT_H 1

namespace pandora { class Pandora; }

/**
 *  @brief  LArDLContent class
 */
class LArDLContent
{
public:
    /**
     *  @brief  Register all the lar dl content algorithms and tools with pandora
     *
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);
};

#endif // #ifndef LAR_DL_CONTENT_H

