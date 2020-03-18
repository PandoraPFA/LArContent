/**
 *  @file   TestAlgorithm.h
 *
 *  @brief  Header file for the test algorithm class.
 *
 *  $Log: $
 */

#ifndef LAR_TEST_ALGORITHM_H
#define LAR_TEST_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTest/TestAlgorithm.h"

namespace lar_content
{

/**
 *  @brief TestAlgorithm class
 */
  class TestAlgorithm : public pandora::Algorithm
{

public:

    /**
     *  @brief  Default constructor
     */
    TestAlgorithm();

private:

    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_clusterListName;
    std::string m_caloHitListName;

};

} //namespace lar_content

#endif //LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H
