/**
 *  @file   larpandoradlcontent/LArTwoDReco/DlClusterAlgorithm.h
 *
 *  @brief  Header file for the deep learning clustering algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_CLUSTER_ALGORITHM_H
#define LAR_DL_CLUSTER_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DlClusterAlgorithm class
 */
class DlClusterAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlClusterAlgorithm() = default;

    virtual ~DlClusterAlgorithm() = default;

private:
    /**
     *  @brief  Produce a training sample and output to CSV file.
     *
     *  @return The StatusCode
     */
    pandora::StatusCode PrepareTrainingSample();

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_caloHitListNames;  ///< The names of the calo hit lists
    std::string m_outputFilePrefix;            ///< The prefix to use for output CSV filenames
    float m_maxEdgeLengthSquared;              ///< The maximum length of an edge
    float m_maxEdgeRatioSquared;               ///< The maximum ratio of the longest to shortest edge of a triangle
    bool m_prune;                              ///< Whether or not to prune the graph
};

} // namespace lar_dl_content

#endif // LAR_DL_CLUSTER_ALGORITHM_H
