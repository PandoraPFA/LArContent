/**
 *  @file   LArContent/include/LArThreeDReco/LArClusterSplitting/ThreeDTrackSplittingAlgorithm.h
 * 
 *  @brief  Header file for the cosmic ray identification algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRACK_SPLITTING_ALGORITHM_H
#define LAR_THREE_D_TRACK_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar
{

/**
 *  @brief  ThreeDTrackSplittingAlgorithm class
 */
class ThreeDTrackSplittingAlgorithm : public pandora::Algorithm
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



    void GetListOfCleanPointingClusters(const pandora::ClusterList *const pClusterList, LArPointingClusterMap &pointingClusterMap,
        LArPointingClusterVertexList &cleanVertexList) const;


    void CollectClusters(const LArPointingClusterMap& inputList, pandora::ClusterList& outputList) const;

 
    void CollectClusters(const LArPointingClusterVertexList& inputList, pandora::ClusterList& outputList) const;


    void CollectMarkers(const LArPointingClusterVertexList& inputList, pandora::CartesianPointList& outputList) const;



    float m_minClusterLength;

    std::string m_clusterListNameU;
    std::string m_clusterListNameV;
    std::string m_clusterListNameW;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTrackSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDTrackSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRACK_SPLITTING_ALGORITHM_H
