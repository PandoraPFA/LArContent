/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/ExampleClusterMopUpAlgorithm.h
 * 
 *  @brief  Header file for the example cluster mop up algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_EXAMPLE_CLUSTER_MOP_UP_ALGORITHM_H
#define LAR_EXAMPLE_CLUSTER_MOP_UP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ExampleClusterMopUpAlgorithm class
 */
class ExampleClusterMopUpAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ExampleClusterMopUpAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_maxMergeDistance;         ///< The maximum distance between clusters for merging to occur
};

} // namespace lar_content

#endif // #ifndef LAR_EXAMPLE_CLUSTER_MOP_UP_ALGORITHM_H
