/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDShowersAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional showers algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_SHOWERS_ALGORITHM_H
#define LAR_THREE_D_SHOWERS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArHelpers/LArClusterHelper.h"

#include "ThreeDBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDShowersAlgorithm class
 */
class ThreeDShowersAlgorithm : public ThreeDBaseAlgorithm<float>
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
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);
    bool ExamineTensor();

    typedef std::map<unsigned int, pandora::CartesianVector> ShowerEdgeMap;

    /**
     *  @brief  GetIncludedHitFraction
     * 
     *  @param  pCluster
     *  @param  minX
     *  @param  maxX
     *  @param  xPitch
     *  @param  edgeMap1
     *  @param  edgeMap2
     */
    float GetIncludedHitFraction(const pandora::Cluster *const pCluster, const float minX, const float maxX, const float xPitch,
        const ShowerEdgeMap &edgeMap1, const ShowerEdgeMap &edgeMap2) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDShowersAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDShowersAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_SHOWERS_ALGORITHM_H
