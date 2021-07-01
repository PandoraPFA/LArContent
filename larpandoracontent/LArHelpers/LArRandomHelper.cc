/**
 *  @file   larpandoracontent/LArHelpers/LArRandomHelper.cc
 *
 *  @brief  Implementation of the random helper class.
 *
 *  $Log: $
 */

#include <random>

namespace lar_content
{

int GetIntsInRange(const int low, const int high, std::mt19937 &eng)
{
    const double answer = eng() / (1.0 + eng.max());
    const int totalRange = high - low + 1;
    return (int) (answer * totalRange) + low;
}

} // namespace lar_content
