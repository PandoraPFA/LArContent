/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoParentAlgorithm.h
 * 
 *  @brief  Header file for the neutrino parent algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_PARENT_ALGORITHM_H
#define LAR_NEUTRINO_PARENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoParentAlgorithm class
 */
class NeutrinoParentAlgorithm : public pandora::Algorithm
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

    typedef std::vector<pandora::ClusterList> ClusterCollectionList;

    pandora::StringVector       m_inputClusterListNames;        ///< The names of the input cluster lists
    std::string                 m_outputClusterListName;        ///< The name of the final output 2D cluster list
    std::string                 m_workingClusterListName;       ///< The name of the working cluster list (current list for daughter algs)

    pandora::StringVector       m_reconstructionAlgNames;       ///< The name of the algorithms to run for each neutrino interaction
    std::string                 m_listManagementAlgName;        ///< The name of the algorithm to perform list management operations
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoParentAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoParentAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_PARENT_ALGORITHM_H
