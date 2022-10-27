/**
 *  @file   larpandoracontent/LArObjects/LArSlice.h
 *
 *  @brief  Header file for a simple class representing a slice.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_H
#define LAR_SLICE_H 1

#include <vector>

#include "Objects/CaloHit.h"

namespace lar_content
{

/**
 *  @brief  Slice class
 */
class Slice
{
public:
    pandora::CaloHitList m_caloHitListU; ///< The u calo hit list
    pandora::CaloHitList m_caloHitListV; ///< The v calo hit list
    pandora::CaloHitList m_caloHitListW; ///< The w calo hit list
};

typedef std::vector<Slice> SliceList;

}

#endif // #ifndef LAR_SLICE_H
