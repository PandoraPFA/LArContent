/**
 *  @file   larpandoracontent/LArReclustering/ThreeDMultiReclusteringAlgorithm.h
 *
 *  @brief  Header file for the reclustering algorithm class.
 *
 *  $Log: $
 */

#ifndef LAR_THREE_D_MULTI_RECLUSTERING_ALGORITHM_H
#define LAR_THREE_D_MULTI_RECLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
  *  @brief  ThreeDMultiReclusteringAlgorithm class
  */
class ThreeDMultiReclusteringAlgorithm : public pandora::Algorithm
{
public:

    /**
     *  @brief  Default constructor
     */
    ThreeDMultiReclusteringAlgorithm();

    /**
    *  @brief  Destructor
    */
    ~ThreeDMultiReclusteringAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string               m_pfoListName;     ///< Name of list of pfos to consider for reclustering
    std::string               m_clusterListName; ///< Name of list of 3D clusters that comprise the pfos
    std::vector<std::string>  m_clusteringAlgs;  ///< The ordered list of clustering algorithms to be used
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_MULTI_RECLUSTERING_ALGORITHM_H
