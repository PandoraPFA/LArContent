/**
 *  @file   LArContent/include/LArUtility/ListDissolutionAlgorithm.h
 * 
 *  @brief  Header file for the list dissolution algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LIST_DISSOLUTION_ALGORITHM_H
#define LAR_LIST_DISSOLUTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ListDissolutionAlgorithm class
 */
class ListDissolutionAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    ListDissolutionAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_pfoListNames;                 ///< The pfo list names
    pandora::StringVector   m_clusterListNames;             ///< The cluster list names
    bool                    m_warnIfClustersUnavailable;    ///< Whether to print warning if attempt made to delete unavailable clusters
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ListDissolutionAlgorithm::Factory::CreateAlgorithm() const
{
    return new ListDissolutionAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_LIST_DISSOLUTION_ALGORITHM_H
