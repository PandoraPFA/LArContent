/**
 *  @file   larpandoracontent/LArHelpers/LArUtilityHelper.h
 *
 *  @brief  Header file for the cluster helper class.
 *
 *  $Log: $
 */
#ifndef LAR_UTILITY_HELPER_H
#define LAR_UTILITY_HELPER_H 1

#include <algorithm>
#include <numeric>
#include <vector>

namespace lar_content
{

/**
 *  @brief  LArUtilityHelper class
 */
class LArUtilityHelper
{
public:
    /**
     *  @brief  Determine the permutation that would apply to the elements of a vector if sorted in ascending order
     *
     *  @param  input the vector for which an order is to be determined
     *  @param  compare the function to use for element comparisons
     *
     *  @return The order of indices determined by the sort operation
     */
    template <typename T, typename Comparison>
    static std::vector<std::size_t> GetSortIndices(const std::vector<T> &input, Comparison &compare);

    /**
     *  @brief  Sort a vector in place based on a supplied index ordering
     *
     *  @param  order the index ordering that should be applied to the vector
     *  @param  vector the vector to be sorted in place
     */
    template <typename T>
    static void SortByIndices(const std::vector<std::size_t> &order, std::vector<T> &vector);
};

// ATTN templated static functions need to be defined in the header file or you get a "used by not defined" error
template <typename T, typename Comparison>
std::vector<std::size_t> LArUtilityHelper::GetSortIndices(const std::vector<T> &input, Comparison &compare)
{
    std::vector<std::size_t> order(input.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j) { return compare(input[i], input[j]); });

    return order;
}

template <typename T>
void LArUtilityHelper::SortByIndices(const std::vector<std::size_t> &order, std::vector<T> &vector)
{
    std::vector<bool> done(vector.size());
    for (std::size_t i = 0; i < vector.size(); ++i)
    {
        if (done[i])
            continue;
        done[i] = true;
        std::size_t idx1{i};
        std::size_t idx2{order[i]};
        // keep swapping indices until we move the matching index into the current index
        while (i != idx2)
        {
            std::swap(vector[idx1], vector[idx2]);
            done[idx2] = true;
            idx1 = idx2;
            idx2 = order[idx2];
        }
    }
}

} // namespace lar_content

#endif // #ifndef LAR_UTILITY_HELPER_H
