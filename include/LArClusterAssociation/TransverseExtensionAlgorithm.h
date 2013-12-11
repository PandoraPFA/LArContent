/**
 *  @file   LArContent/include/LArClusterAssociation/TransverseExtensionAlgorithm.h
 * 
 *  @brief  Header file for the hello world algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRANSVERSE_EXTENSION_ALGORITHM_H
#define LAR_TRANSVERSE_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

namespace lar
{

/**
 *  @brief  TransverseExtensionAlgorithm class
 */
class TransverseExtensionAlgorithm : public pandora::Algorithm
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


    void SortInputClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &longVector, pandora::ClusterVector &shortVector) const;



    
    bool IsAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const pandora::Cluster *const pCluster) const;


    bool IsEndAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const pandora::Cluster *const pCluster) const;



    bool IsMidAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const pandora::Cluster *const pCluster) const;


  



    float m_minClusterLength;

    float m_maxTransverseDisplacement;

    float m_maxLongitudinalDisplacement;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TransverseExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new TransverseExtensionAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRANSVERSE_EXTENSION_ALGORITHM_H
