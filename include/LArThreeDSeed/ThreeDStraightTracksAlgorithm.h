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
     *  @brief  
     * 
     * 
     * 
     */
    void SelectInputClusters() const;

    void ModifyInputClusters();
    void InitializeTensor();
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);
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
