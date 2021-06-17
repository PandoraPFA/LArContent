/**
 *  @file   larpandoradlcontent/LArTwoDReco/StreamingAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_STREAMING_ALGORITHM_H
#define LAR_STREAMING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  StreamingAlgorithm class
 */
class StreamingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    StreamingAlgorithm();

    virtual ~StreamingAlgorithm();

private:
    typedef std::map<std::string, pandora::StringVector> StreamAlgorithmMap;
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_outputListName;            ///< The name of the output list
    std::string m_listType;                  ///< The type of the input lists (currently only Cluster is supported)
    pandora::StringVector m_inputListNames;  ///< The names of the input lists
    pandora::StringVector m_outputListNames; ///< Names of the output lists if not combining into a single list at the end
    StreamAlgorithmMap m_streamAlgorithmMap; ///< A map from individual streams to the algorithms that stream should run
};

} // namespace lar_content

#endif // LAR_STREAMING_ALGORITHM_H
