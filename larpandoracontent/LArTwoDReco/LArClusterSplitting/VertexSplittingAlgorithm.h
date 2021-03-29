/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h
 *
 *  @brief  Header file for the vertex splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_SPLITTING_ALGORITHM_H
#define LAR_VERTEX_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAlgorithm.h"

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexSplittingAlgorithm class
 */
class VertexSplittingAlgorithm : public TwoDSlidingFitSplittingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, pandora::CartesianVector &splitPosition) const;

    float m_splitDisplacementSquared;  ///< Maximum displacement squared
    float m_vertexDisplacementSquared; ///< Maximum displacement squared
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SPLITTING_ALGORITHM_H
