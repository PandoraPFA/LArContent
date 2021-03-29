/**
 *  @file   larpandoracontent/LArControlFlow/PostProcessingAlgorithm.h
 *
 *  @brief  Header file for the post processing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_POST_PROCESSING_ALGORITHM_H
#define LAR_POST_PROCESSING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PostProcessingAlgorithm class
 */
class PostProcessingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PostProcessingAlgorithm();

private:
    pandora::StatusCode Reset();
    pandora::StatusCode Run();

    /**
     *  @brief  Rename a list of relevant type with specified name - the new name will be the old name with appended list counter
     *
     *  @param  oldListName the old list name
     */
    template <typename T>
    pandora::StatusCode RenameList(const std::string &oldListName) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;     ///< The list of pfo list names
    pandora::StringVector m_clusterListNames; ///< The list of cluster list names
    pandora::StringVector m_vertexListNames;  ///< The list of vertex list names
    pandora::StringVector m_caloHitListNames; ///< The list of calo hit list names

    std::string m_currentPfoListReplacement; ///< The name of the pfo list to replace the current list

    unsigned int m_listCounter; ///< The counter appended to output (and replacement current) list names and reset each event
};

} // namespace lar_content

#endif // #ifndef LAR_POST_PROCESSING_ALGORITHM_H
