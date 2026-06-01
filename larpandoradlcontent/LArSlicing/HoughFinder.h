/**
 * @file   larpandoradlcontent/LArSlicing/HoughFinder.h
 * @brief  Header file for the Fast Hough Transform vertex finder
 */
#ifndef LAR_FAST_HOUGH_FINDER_H
#define LAR_FAST_HOUGH_FINDER_H 1

#include "Objects/CartesianVector.h"
#include <vector>

namespace lar_content
{

class FastHoughFinder
{
public:
    /**
     * @brief Constructor
     *
     * @param thresholds List of distance bin edge thresholds
     * @param scalingFactor The scaling factor to convert network output back to physical cm
     * @param tolerance The distance tolerance (in cm) for a vote to count
     * @param minVotes Minimum number of votes required to keep a vertex candidate
     * @param nmsRadius Radius (in cm) for Non-Maximum Suppression (NMS)
     */
    FastHoughFinder(const std::vector<float> &thresholds, const float scalingFactor = 400.0f, const float tolerance = 25.0f,
        const int minVotes = 3, const float nmsRadius = 35.0f);

    /**
     * @brief Finds vertices based on hit positions and network distance logits
     *
     * @param hitPositions Vector of physical hit positions
     * @param logits Vector of network logits for each hit (flattened, size = N_hits * N_classes)
     *
     * @return Vector of found vertex positions
     */
    std::vector<pandora::CartesianVector> Fit(const std::vector<pandora::CartesianVector> &hitPositions, const std::vector<float> &logits) const;

private:
    std::vector<float> m_thresholds;
    std::vector<float> m_binCenters;
    float m_scalingFactor;
    float m_tolerance;
    int m_minVotes;
    float m_nmsRadiusSq; ///< Stored as squared radius for faster comparison
};

} // namespace lar_content

#endif // LAR_FAST_HOUGH_FINDER_H
