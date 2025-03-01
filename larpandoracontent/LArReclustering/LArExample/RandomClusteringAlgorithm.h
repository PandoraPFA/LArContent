/**
 *  @file   larpandoracontent/LArReclustering/RandomClusteringAlgorithm.h
 *
 *  @brief  Header file for the random clustering algorithm class. The class exists as a test for
 *          the reclustering framework.
 *
 *  $Log: $
 */

#ifndef RANDOM_CLUSTERING_ALGORITHM_H
#define RANDOM_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
  *  @brief  RandomClusteringAlgorithm class
  */
class RandomClusteringAlgorithm : public pandora::Algorithm
{
public:

    /**
     *  @brief  Default constructor
     */
    RandomClusteringAlgorithm();

    /**
    *  @brief  Destructor
    */
    ~RandomClusteringAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_numNewClusters; ///< Number of new random clusters to generate from original cluster list
};

} // namespace lar_content

#endif // #ifndef RANDOM_CLUSTERING_ALGORITHM_H
