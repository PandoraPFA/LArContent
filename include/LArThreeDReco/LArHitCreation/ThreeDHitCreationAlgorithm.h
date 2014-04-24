/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional hit creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
#define LAR_THREE_D_HIT_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDHitCreationAlgorithm::Algorithm class
 */
class ThreeDHitCreationAlgorithm : public pandora::Algorithm
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
     *  @brief  Get the primary u, v and w clusters in a provided pfo
     * 
     *  @param  pPfo address of the pfo
     *  @param  pClusterU to receive the address of cluster u
     *  @param  pClusterV to receive the address of cluster v
     *  @param  pClusterW to receive the address of cluster w
     */
    void GetClusters(const pandora::ParticleFlowObject *const pPfo, pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV, pandora::Cluster *&pClusterW) const;

    /**
     *  @brief  Create three dimensional hits, using an input list of two dimensional hits and two associated sliding fit results
     * 
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  fitResult1 the first sliding fit result
     *  @param  fitResult2 the second sliding fit result
     *  @param  newThreeDHits to receive the new three dimensional hits
     *  @param  omittedTwoDHits to receive the two dimensional hits for which no three dimensional hits could be created
     */
    void CreateThreeDHits(const pandora::CaloHitList &inputTwoDHits, const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
        pandora::CaloHitList &newThreeDHits, pandora::CaloHitList &omittedTwoDHits) const;

    /**
     *  @brief  Create a new three dimensional cluster, using a list of provided three dimensional hits
     * 
     *  @param  caloHitList the list of three dimensional hits
     *  @param  pCluster to receive the address of the new cluster
     */
    void CreateThreeDCluster(const pandora::CaloHitList &caloHitList, pandora::Cluster *&pCluster) const;

    /**
     *  @brief  Create new three dimensional hits for those hits previously omitted, using extrapolation of the three dimensional cluster
     * 
     *  @param  omittedHits the list of previously omitted two dimensional hits
     *  @param  pCluster address of the three dimensional cluster
     *  @param  extrapolatedHits to receive the list of extrapolated hits
     */
    void CreateExtrapolatedHits(const pandora::CaloHitList &omittedHits, pandora::Cluster *pCluster, pandora::CaloHitList &extrapolatedHits) const;

    std::string     m_inputPfoListName;         ///< The name of the input pfo list
    std::string     m_outputCaloHitListName;    ///< The name of the output calo hit list
    std::string     m_outputClusterListName;    ///< The name of the output cluster list

    unsigned int    m_slidingFitWindow;         ///< The layer window for the sliding linear fits
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDHitCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDHitCreationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
