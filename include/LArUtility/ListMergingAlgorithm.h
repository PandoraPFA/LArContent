/**
 *  @file   LArContent/include/LArUtility/ListMergingAlgorithm.h
 * 
 *  @brief  Header file for the list merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LIST_MERGING_ALGORITHM_H
#define LAR_LIST_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ListMergingAlgorithm class
 */
class ListMergingAlgorithm : public pandora::Algorithm
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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_sourceClusterListNames;   ///< The source cluster list names
    pandora::StringVector   m_targetClusterListNames;   ///< The target cluster list names

    pandora::StringVector   m_sourcePfoListNames;       ///< The source pfo list names
    pandora::StringVector   m_targetPfoListNames;       ///< The target pfo list names
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ListMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ListMergingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_LIST_MERGING_ALGORITHM_H
