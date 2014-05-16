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



    pandora::StatusCode FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, pandora::CartesianVector &splitPosition,
        pandora::CartesianVector &splitDirection1, pandora::CartesianVector &splitDirection2) const; 




    pandora::StatusCode ConfirmSplitPosition(const TwoDSlidingFitResult &slidingFitResult, const pandora::CartesianVector &splitPosition,
        const pandora::CartesianVector &splitDirection1, const pandora::CartesianVector &splitDirection2) const; 


    unsigned int   m_halfWindowLayers;           ///<
    float          m_samplingPitch;              ///< 
    float          m_clusterMinLength;           ///<
    float          m_maxSplitCosTheta;           ///<
    float          m_minMergeCosTheta;           ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRaySplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRaySplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SPLITTING_ALGORITHM_H
