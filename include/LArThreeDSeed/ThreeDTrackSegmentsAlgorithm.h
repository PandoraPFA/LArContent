/**
 *  @file   LArContent/include/ThreeDSeed/ThreeDTrackSegmentsAlgorithm.h
 * 
 *  @brief  Header file for the three dimension track segments algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRACK_SEGMENTS_ALGORITHM_H
#define LAR_THREE_D_TRACK_SEGMENTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "ThreeDBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDTrackSegmentsAlgorithm class
 */
class ThreeDTrackSegmentsAlgorithm : public ThreeDBaseAlgorithm
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

inline pandora::Algorithm *ThreeDTrackSegmentsAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDTrackSegmentsAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRACK_SEGMENTS_ALGORITHM_H
