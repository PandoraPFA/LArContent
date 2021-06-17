/**
 *  @file   larpandoracontent/LArTwoDReco/StreamSelectionAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_STREAM_SELECTION_ALGORITHM_H
#define LAR_STREAM_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  StreamSelectionAlgorithm class
 */
class StreamSelectionAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    StreamSelectionAlgorithm();

    virtual ~StreamSelectionAlgorithm() = default;

protected:
    typedef std::map<std::string, pandora::ClusterList> ClusterListMap;

    /**
     *  @brief  Allocate a cluster to the appropriate streams.
     *
     *  @param  pCluster The cluster to allocate to a stream
     *
     *  @return The StatusCode
     */
    virtual pandora::StatusCode AllocateToStreams(const pandora::Cluster *const pCluster) = 0;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_listType;            ///< The type of the input lists (currently only Cluster is supported)
    pandora::StringVector m_listNames; ///< The name of the output lists
    ClusterListMap m_clusterListMap;   ///< The map from cluster list names to cluster lists

private:
    pandora::StatusCode Run();
};

} // namespace lar_content

#endif // LAR_STREAM_SELECTION_ALGORITHM_H
