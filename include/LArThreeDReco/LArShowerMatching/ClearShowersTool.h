/**
 *  @file   LArContent/include/LArThreeDReco/LArShowerMatching/ClearShowersTool.h
 * 
 *  @brief  Header file for the clear showers tool class.
 * 
 *  $Log: $
 */
#ifndef CLEAR_SHOWERS_TOOL_H
#define CLEAR_SHOWERS_TOOL_H 1

#include "LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClearShowersTool class
 */
class ClearShowersTool : public ShowerTensorTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    ClearShowersTool();

    /**
     *  @brief  Whether a large shower-like element shares clusters with any other long elements
     * 
     *  @param  iIter specifies the large element under consideration
     *  @param  iteratorList list of iterators to other large elements
     * 
     *  @return boolean
     */
    static bool HasLargeDirectConnections(IteratorList::const_iterator iIter, const IteratorList &iteratorList);

    /**
     *  @brief  Whether a large shower-like element is significantly larger that other elements with which it shares a cluster
     * 
     *  @param  iIter specifies the large element under consideration
     *  @param  elementList the full list of connected tensor elements
     *  @param  minMatchedSamplingPointRatio the min ratio between 1st and 2nd highest msps for simple ambiguity resolution
     *  @param  minMatchedSamplingPointRatio the min ratio between 1st and 2nd highest x-overlap spans for simple ambiguity resolution
     *  @param  usedClusters the list of clusters already marked as to be added to a pfo
     */
    static bool IsLargerThanDirectConnections(IteratorList::const_iterator iIter, const TensorType::ElementList &elementList,
        const unsigned int minMatchedSamplingPointRatio, const float minXOverlapSpanRatio, const pandora::ClusterList &usedClusters);

    bool Run(ThreeDShowersAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find clear shower matches, hidden by simple ambiguities in the tensor
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  protoParticleVector to receive the list of proto particles
     */
    void FindClearShowers(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Select a list of large shower-like elements from a set of connected tensor elements
     * 
     *  @param  elementList the full list of connected tensor elements
     *  @param  usedClusters the list of clusters already marked as to be added to a pfo
     *  @param  iteratorList to receive a list of iterators to large shower-like elements
     */
    void SelectLargeShowerElements(const TensorType::ElementList &elementList, const pandora::ClusterList &usedClusters,
        IteratorList &iteratorList) const;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for particle creation
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for particle creation
    float           m_minXOverlapFraction;              ///< The min x overlap fraction (in each view) for particle creation
    unsigned int    m_minMatchedSamplingPointRatio;     ///< The min ratio between 1st and 2nd highest msps for simple ambiguity resolution
    float           m_minXOverlapSpanRatio;             ///< The min ratio between 1st and 2nd highest x-overlap spans for simple ambiguity resolution
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearShowersTool::Factory::CreateAlgorithmTool() const
{
    return new ClearShowersTool();
}

} // namespace lar_content

#endif // #ifndef CLEAR_SHOWERS_TOOL_H
