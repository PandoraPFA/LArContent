/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/OvershootTracksTool.h
 * 
 *  @brief  Header file for the overshoot tracks tool class.
 * 
 *  $Log: $
 */
#ifndef OVERSHOOT_TRACKS_TOOL_H
#define OVERSHOOT_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  OvershootTracksTool class
 */
class OvershootTracksTool : public TensorManipulationTool
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

private:
    /**
     *  @brief  Split particle class
     */
    class SplitParticle
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  elementA the tensor element A
         *  @param  elementB the tensor element B
         */
        SplitParticle(const TensorType::Element &elementA, const TensorType::Element &elementB);

        pandora::Cluster           *m_pCommonCluster;       ///< Address of the common cluster
        pandora::Cluster           *m_pClusterA1;           ///< Address of cluster in element A, view 1
        pandora::Cluster           *m_pClusterA2;           ///< Address of cluster in element A, view 2
        pandora::Cluster           *m_pClusterB1;           ///< Address of cluster in element B, view 1
        pandora::Cluster           *m_pClusterB2;           ///< Address of cluster in element B, view 2
        pandora::CartesianVector    m_splitPosition;        ///< The candidate split position
    };

    typedef std::vector<SplitParticle> SplitParticleList;

    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Find overshoot tracks, by receiving a list of details about possible required cluster splits
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  splitParticleList to receive the split particle list
     */
    void FindOvershootTracks(const TensorType &overlapTensor, SplitParticleList &splitParticleList) const;

    /**
     *  @brief  Apply the changes cached in a split particle list and update the tensor accordingly
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  splitParticleList the split particle list
     * 
     *  @return whether changes to the tensor have been made
     */
    bool ApplyChanges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const SplitParticleList &splitParticleList) const;

    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Select elements representing possible components of a two particles, merged due to an overshoot in the clustering
     * 
     *  @param  eIter iterator to a candidate element
     *  @param  elementList the provided element list
     *  @param  usedClusters the list of used clusters
     *  @param  iteratorList to receive a list of iterators to relevant elements
     */
    void SelectOvershootElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
        const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Get split particle objects, identifying required splits for clusters merged due to overshoots in the clustering.
     * 
     *  @param  iteratorList list of iterators to relevant elements
     *  @param  splitParticleList to be populated with split particles for 
     */
    void GetSplitParticles(const IteratorList &iteratorList, SplitParticleList &splitParticleList) const;

    /**
     *  @brief  Whether a provided (iterator to a) tensor element passes the selection cuts for overshoot identification
     * 
     *  @param  eIter the iterator to the tensor element
     *  @param  usedClusters the list of used clusters
     */
    bool PassesElementCuts(TensorType::ElementList::const_iterator eIter, const pandora::ClusterList &usedClusters) const;

    /**
     *  @brief  Whether a pair of vertices pass longitudinal projection cuts
     * 
     *  @param  vertexA vertex from cluster in tensor element a
     *  @param  vertexB vertex from cluster in tensor element b
     */
    bool PassesVertexCuts(const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB) const;

    /**
     *  @brief  Set split position for a provided split particle
     * 
     *  @param  vertexA1 vertex for tensor element a in view 1
     *  @param  vertexA2 vertex for tensor element a in view 2
     *  @param  vertexB1 vertex for tensor element b in view 1
     *  @param  vertexB2 vertex for tensor element b in view 2
     *  @param  splitParticle the split particle
     */
    void SetSplitPosition(const LArPointingCluster::Vertex &vertexA1, const LArPointingCluster::Vertex &vertexA2,
        const LArPointingCluster::Vertex &vertexB1, const LArPointingCluster::Vertex &vertexB2, SplitParticle &splitParticle) const;

    typedef std::map<pandora::Cluster*, pandora::CartesianPointList> SplitPositionMap;

    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for use as a key tensor element
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for use as a key tensor element
    float           m_minLongitudinalImpactParameter;   ///< The min longitudinal impact parameter for connecting accompanying clusters
    float           m_maxVertexXSeparation;             ///< The max separation between accompanying clusters vertex x positions to make split
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *OvershootTracksTool::Factory::CreateAlgorithmTool() const
{
    return new OvershootTracksTool();
}

} // namespace lar

#endif // #ifndef OVERSHOOT_TRACKS_TOOL_H
