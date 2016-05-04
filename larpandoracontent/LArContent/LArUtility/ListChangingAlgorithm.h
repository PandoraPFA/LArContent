/**
 *  @file   LArContent/include/LArUtility/ListChangingAlgorithm.h
 * 
 *  @brief  Header file for the list changing algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LIST_CHANGING_ALGORITHM_H
#define LAR_LIST_CHANGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ListChangingAlgorithm::Algorithm class
 */
class ListChangingAlgorithm : public pandora::Algorithm
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

    std::string     m_caloHitListName;  ///< The calo hit list name to set as the current calo hit list
    std::string     m_clusterListName;  ///< The cluster list name to set as the current cluster list
    std::string     m_vertexListName;   ///< The vertex list name to set as the current vertex list
    std::string     m_pfoListName;      ///< The pfo list name to set as the current pfo list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ListChangingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ListChangingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_LIST_CHANGING_ALGORITHM_H
