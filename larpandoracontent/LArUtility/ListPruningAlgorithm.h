/**
 *  @file   larpandoracontent/LArUtility/ListPruningAlgorithm.h
 *
 *  @brief  Header file for the list pruning algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LIST_PRUNING_ALGORITHM_H
#define LAR_LIST_PRUNING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ListPruningAlgorithm class
 */
class ListPruningAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ListPruningAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;     ///< The pfo list names
    pandora::StringVector m_clusterListNames; ///< The cluster list names
    pandora::StringVector m_vertexListNames;  ///< The vertex list names
    bool m_warnIfObjectsUnavailable;          ///< Whether to print warning if attempt made to delete unavailable objects
};

} // namespace lar_content

#endif // #ifndef LAR_LIST_PRUNING_ALGORITHM_H
