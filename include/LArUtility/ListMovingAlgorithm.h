/**
 *  @file   LArContent/include/LArUtility/ListMovingAlgorithm.h
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

    pandora::StringVector   m_clusterListNames;     ///< The list of cluster list names
    pandora::StringVector   m_vertexListNames;      ///< The list of vertex list names
    pandora::StringVector   m_pfoListNames;         ///< The list of pfo list names
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ListMovingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ListMovingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_LIST_MOVING_ALGORITHM_H
