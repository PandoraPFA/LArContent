/**
 *  @file   LArContent/include/LArCheating/CheatingCosmicRayShowerMatchingAlg.h
 * 
 *  @brief  Header file for the cosmic ray shower matching cheater class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_COSMIC_RAY_SHOWER_MATCHING_ALG_H
#define LAR_CHEATING_COSMIC_RAY_SHOWER_MATCHING_ALG_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  CheatingCosmicRayShowerMatchingAlg class
 */
class CheatingCosmicRayShowerMatchingAlg : public pandora::Algorithm
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

    /**
     *  @brief  Perform cosmic ray shower matching for a specific cluster in a pfo
     * 
     *  @param  clusterListNames the list of cluster list names
     *  @param  pPfoCluster the pfo cluster of interest
     *  @param  pPfoList the address of the list containing the cosmic ray pfos
     */
    pandora::StatusCode CosmicRayShowerMatching(const pandora::StringVector &clusterListNames, const pandora::Cluster *const pPfoCluster,
        pandora::ParticleFlowObject *pPfo) const;

    std::string             m_inputPfoListName;           ///< The input pfo list name

    pandora::StringVector   m_inputClusterListNamesU;     ///< The input cluster list names for the u view
    pandora::StringVector   m_inputClusterListNamesV;     ///< The input cluster list names for the v view
    pandora::StringVector   m_inputClusterListNamesW;     ///< The input cluster list names for the w view
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CheatingCosmicRayShowerMatchingAlg::Factory::CreateAlgorithm() const
{
    return new CheatingCosmicRayShowerMatchingAlg();
}

} // namespace lar

#endif // #ifndef LAR_CHEATING_COSMIC_RAY_SHOWER_MATCHING_ALG_H
