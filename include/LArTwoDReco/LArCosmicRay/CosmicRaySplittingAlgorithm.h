/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h
 * 
 *  @brief  Header file for the cosmic ray splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_SPLITTING_ALGORITHM_H
#define LAR_COSMIC_RAY_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

namespace lar
{

/**
 *  @brief  CosmicRaySplittingAlgorithm class
 */
class CosmicRaySplittingAlgorithm : public pandora::Algorithm
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



    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    float          m_clusterMinLength;           ///<
    unsigned int   m_halfWindowLayers;           ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRaySplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRaySplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SPLITTING_ALGORITHM_H
