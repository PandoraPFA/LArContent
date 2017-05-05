/**
 *  @file   larpandoracontent/LArUtility/ListMovingAlgorithm.h
 * 
 *  @brief  Header file for the list moving algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LIST_MOVING_ALGORITHM_H
#define LAR_LIST_MOVING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ListMovingAlgorithm class
 */
class ListMovingAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_pfoListNames;         ///< The list of pfo list names
    pandora::StringVector   m_clusterListNames;     ///< The list of cluster list names
    pandora::StringVector   m_vertexListNames;      ///< The list of vertex list names

    std::string             m_prefix;               ///< The prefix to add to the list names, in order to yield the output list names
};

} // namespace lar_content

#endif // #ifndef LAR_LIST_MOVING_ALGORITHM_H
