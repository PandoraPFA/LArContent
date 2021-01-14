/**
 *  @file   larpandoracontent/LArHelpers/LArRandomHelper.h
 *
 *  @brief  Header file for the random helper class.
 *
 *  $Log: $
 */

#ifndef LAR_RANDOM_HELPER_H
#define LAR_RANDOM_HELPER_H 1

#include <random>

namespace lar_content
{

/**
 *  @brief  Specific implementation of std::uniform_int_distribution. The default can vary between compilers, which would lead
 *          to inconsistent randomness between gcc and clang. This gives consistent results for both gcc and clang, as
 *          std::mt19937 will return the same value for both.
 *
 *  @param  low  Lowest number in range to pick.
 *  @param  high Highest number in range to pick.
 *  @param  eng The MT19937 to give a random number to scale.
 *
 *  @return int An integer between the limits, inclusive.
 */
int GetIntsInRange(const int low, const int high, std::mt19937 &eng);

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_GEOMETRY_HELPER_H
