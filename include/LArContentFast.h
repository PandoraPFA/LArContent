/**
 *  @file   LArContent/include/LArContentFast.h
 * 
 *  @brief  Header file detailing faster versions of algorithms in the LArContent library, using e.g. KD-trees and *relying on c++11*
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_FAST_H
#define LAR_CONTENT_FAST_H 1

#include "LArContentFast/ListPreparationAlgorithmFast.h"
#include "LArContentFast/LongitudinalAssociationAlgorithmFast.h"
#include "LArContentFast/TransverseAssociationAlgorithmFast.h"
#include "LArContentFast/VertexSelectionAlgorithmFast.h"

/**
 *  @brief  LArContentFast class
 */
class LArContentFast
{
public:
    #define LAR_ALGORITHM_FAST_LIST(d)                                                                                          \
        d("LArListPreparationFast",                    lar_content_fast::ListPreparationAlgorithm::Factory)                     \
        d("LArLongitudinalAssociationFast",            lar_content_fast::LongitudinalAssociationAlgorithm::Factory)             \
        d("LArTransverseAssociationFast",              lar_content_fast::TransverseAssociationAlgorithm::Factory)               \
        d("LArVertexSelectionFast",                    lar_content_fast::VertexSelectionAlgorithm::Factory)

    /**
     *  @brief  Register all the lar content fast algorithms with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContentFast::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LAR_ALGORITHM_FAST_LIST(PANDORA_REGISTER_ALGORITHM);

    return pandora::STATUS_CODE_SUCCESS;
}

#endif // #ifndef LAR_CONTENT_FAST_H
