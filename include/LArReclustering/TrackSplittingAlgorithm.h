/**
 *  @file   LArContent/include/LArReclustering/TrackSplittingAlgorithm.h
 * 
 *  @brief  Header file for the track splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRACK_SPLITTING_ALGORITHM_H
#define LAR_TRACK_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  TrackSplittingAlgorithm class
 */
class TrackSplittingAlgorithm : public pandora::Algorithm
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

    std::string             m_inputSeedClusterListName;         ///< The seed cluster list name
    std::string             m_showerSeedClusterListName;        ///< The shower seed cluster list name
    std::string             m_trackSeedClusterListName;         ///< The track seed cluster list name
    std::string             m_nonSeedClusterListName;           ///< The non seed cluster list name

    std::string             m_clusteringAlgorithmName;          ///< The name of the clustering algorithm to run
    pandora::StringVector   m_clusterAssociationAlgorithms;     ///< The ordered list of cluster association algorithms to be used
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TrackSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new TrackSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRACK_SPLITTING_ALGORITHM_H
