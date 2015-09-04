/**
 *  @file   LArContent/include/LArUtility/NeutrinoListManagementAlgorithm.h
 * 
 *  @brief  Header file for the neutrino list management algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_LIST_MANAGEMENT_ALGORITHM_H
#define LAR_NEUTRINO_LIST_MANAGEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoListManagementAlgorithm class
 */
class NeutrinoListManagementAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Perform list management for pfos
     */
    void PfoListManagement() const;

    /**
     *  @brief  Perform list management for vertices
     */
    void VertexListManagement() const;

    /**
     *  @brief  Perform list management for clusters
     */
    void ClusterListManagement() const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_neutrinoPfoListName;          ///< The neutrino pfo list name
    std::string     m_trackPfoListName;             ///< The track pfo list name
    std::string     m_showerPfoListName;            ///< The shower pfo list name
    std::string     m_neutrinoVertexListName;       ///< The neutrino vertex list name
    std::string     m_clusterListName;              ///< The cluster list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoListManagementAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoListManagementAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_LIST_MANAGEMENT_ALGORITHM_H
