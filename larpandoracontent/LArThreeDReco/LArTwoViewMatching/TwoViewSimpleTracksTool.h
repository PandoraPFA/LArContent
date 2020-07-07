/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/TwoViewSimpleTracksTool.h
 *
 *  @brief  Header file for the simple showers tool class.
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_SIMPLE_TRACKS_TOOL_H
#define TWO_VIEW_SIMPLE_TRACKS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArShowerMatching/ThreeViewShowersAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewSimpleTracksTool class
 */
class TwoViewSimpleTracksTool : public ShowerTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewSimpleTracksTool();

    bool Run(ThreeViewShowersAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    /**
     *  @brief  Find best shower match as a simple way to (try to) resolve ambiguities in the tensor
     *
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindBestShower(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Whether a provided (iterator to a) tensor element passes the selection cuts for particle creation
     *
     *  @param  eIter the iterator to the tensor element
     */
    bool PassesElementCuts(TensorType::ElementList::const_iterator eIter) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (in each view) for particle creation
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_SIMPLE_TRACKS_TOOL_H
