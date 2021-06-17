/**
 *  @file   larpandoradlcontent/LArTwoDReco/DlTrackShowerStreamSelectionAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_TRACK_SHOWER_STREAM_SELECTION_ALGORITHM_H
#define LAR_DL_TRACK_SHOWER_STREAM_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/StreamSelectionAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DlTrackShowerStreamSelectionAlgorithm class
 */
class DlTrackShowerStreamSelectionAlgorithm : public lar_content::StreamSelectionAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlTrackShowerStreamSelectionAlgorithm() = default;

    virtual ~DlTrackShowerStreamSelectionAlgorithm() = default;

protected:
    /**
     *  @brief  Allocate a cluster to the appropriate streams.
     *
     *  @param  pCluster The cluster to allocate to a stream
     *
     *  @return The StatusCode
     */
    virtual pandora::StatusCode AllocateToStreams(const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_trackListName;  ///< The name of the track list
    std::string m_showerListName; ///< The name of the shower list
};

} // namespace lar_dl_content

#endif // LAR_DL_TRACK_SHOWER_STREAM_SELECTION_ALGORITHM_H
