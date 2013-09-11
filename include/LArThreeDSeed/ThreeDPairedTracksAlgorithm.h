/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDPairedTracksAlgorithm.h
 * 
 *  @brief  Header file for the 3D paired tracks algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_PAIRED_TRACKS_ALGORITHM_H
#define LAR_THREE_D_PAIRED_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArOverlapTensor.h"

#include "PandoraMonitoringApi.h"

namespace lar
{

/**
 *  @brief  ThreeDPairedTracksAlgorithm class
 */
class ThreeDPairedTracksAlgorithm : public pandora::Algorithm
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
    void CalculateOverlapResult(const pandora::ClusterVector &clusterVectorU, const pandora::ClusterVector &clusterVectorV, const pandora::ClusterVector &clusterVectorW);
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

    /**
     *  @brief  Identify next particle and try to build it; return true if a particle if created
     * 
     *  @return boolean
     */
    bool BuildNextParticle() const;

    OverlapTensor<TrackOverlapResult> m_overlapTensor;          ///< The overlap tensor

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string                 m_outputPfoListName;            ///< The output pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDPairedTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDPairedTracksAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_PAIRED_TRACKS_ALGORITHM_H
