/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDStraightTracksAlgorithm.h
 * 
 *  @brief  Header file for the three dimension straight tracksalgorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_STRAIGHT_TRACKS_ALGORITHM_H
#define LAR_THREE_D_STRAIGHT_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

#include "ThreeDBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDStraightTracksAlgorithm class
 */
class ThreeDStraightTracksAlgorithm : public ThreeDBaseAlgorithm
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
    void SelectInputClusters();

    /**
     *  @brief  Select input clusters for 2D->3D operations
     * 
     *  @param  pClusterList address of the input cluster list
     *  @param  clusterVector to receive clusters to be used in 2D->3D operations
     */
    void SelectInputClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    void ModifyInputClusters();
    void InitializeTensor();
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

    /**
     *  @brief  Calculate overlap result for special case with clusters at constant x
     * 
     *  @param  slidingFitResultU sliding fit result for u cluster
     *  @param  slidingFitResultV sliding fit result for v cluster
     *  @param  slidingFitResultW sliding fit result for w cluster
     */
    void CalculateConstantXOverlapResult(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
        const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW);

    bool ExamineTensor();
    void UpdateTensor();
    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    OverlapTensor<float>    m_overlapTensor;        ///< The overlap tensor
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDStraightTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDStraightTracksAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_STRAIGHT_TRACKS_ALGORITHM_H
