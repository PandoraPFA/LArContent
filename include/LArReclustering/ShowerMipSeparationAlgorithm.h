/**
 *  @file   LArContent/include/LArReclustering/ShowerMipSeparationAlgorithm.h
 * 
 *  @brief  Header file for the shower-mip separation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SHOWER_MIP_SEPARATION_ALGORITHM_H
#define LAR_SHOWER_MIP_SEPARATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ShowerMipSeparationAlgorithm class
 */
class ShowerMipSeparationAlgorithm : public pandora::Algorithm
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

    std::string     m_clusteringAlgorithmName;          ///< The name of the clustering algorithm to run
    std::string     m_associationAlgorithmName;         ///< The name of the cluster association algorithm to run
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ShowerMipSeparationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ShowerMipSeparationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_SHOWER_MIP_SEPARATION_ALGORITHM_H
