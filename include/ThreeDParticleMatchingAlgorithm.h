/**
 *  @file   ThreeDParticleMatchingAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional particle creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_PARTICLE_MATCHING_ALGORITHM_H
#define LAR_THREE_D_PARTICLE_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  ThreeDParticleMatchingAlgorithm class
 */
class ThreeDParticleMatchingAlgorithm : public pandora::Algorithm
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


    void ProcessSingleView(const pandora::ClusterList* const pClusterList, LArPointingClusterVertexList& pointingClusterVertexList );

    void GetListOfCleanClusters(const pandora::ClusterList* const pClusterList, pandora::ClusterVector &clusterVector);  

    std::string         m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string         m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string         m_inputClusterListNameW;        ///< The name of the view W cluster list
   
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDParticleMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDParticleMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_PARTICLE_MATCHING_ALGORITHM_H
