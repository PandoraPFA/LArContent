/**
 *  @file   LArContent/include/LArThreeDReco/LArTransverseTrackMatching/ThreeDKinkBaseTool.h
 * 
 *  @brief  Header file for the three d kink base tool
 * 
 *  $Log: $
 */
#ifndef THREE_D_KINK_BASE_TOOL_H
#define THREE_D_KINK_BASE_TOOL_H 1

#include "LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ThreeDKinkBaseTool class
 */
class ThreeDKinkBaseTool : public TransverseTensorTool
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  nCommonClusters the number of common clusters to select
     */
    ThreeDKinkBaseTool(const unsigned int nCommonClusters);

    /**
     *  @brief  Destructor
     */
    virtual ~ThreeDKinkBaseTool();

    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);

protected:
    /**
     *  @brief  Modification class
     */
    class Modification
    {
    public:
        SplitPositionMap            m_splitPositionMap;     ///< The split position map
        ClusterMergeMap             m_clusterMergeMap;      ///< The cluster merge map
        pandora::ClusterList        m_affectedClusters;     ///< The list of affected clusters
    };

    typedef std::vector<Modification> ModificationList;

    /**
     *  @brief  Whether a provided (iterator to a) tensor element passes the selection cuts for overshoot identification
     * 
     *  @param  eIter the iterator to the tensor element
     *  @param  usedClusters the list of used clusters
     */
    virtual bool PassesElementCuts(TensorType::ElementList::const_iterator eIter, const pandora::ClusterList &usedClusters) const;

    /**
     *  @brief  Get modification objects for a specific elements of the tensor, identifying required splits and merges for clusters
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  iteratorList list of iterators to relevant tensor elements
     *  @param  modificationList to be populated with modifications
     */
    virtual void GetIteratorListModifications(ThreeDTransverseTracksAlgorithm *pAlgorithm, const IteratorList &iteratorList,
        ModificationList &modificationList) const = 0;

    /**
     *  @brief  Get a sampling point in x that is common to sliding linear fit objects in three views
     * 
     *  @param  splitPosition1 the split position in view 1
     *  @param  isForwardInX whether to work forwards (or backwards) in x
     *  @param  fitResult1 the sliding fit result in view 1
     *  @param  fitResult2 the sliding fit result in view 2
     *  @param  fitResult3 the sliding fit result in view 3
     * 
     *  @return the sampling point
     */
    float GetXSamplingPoint(const pandora::CartesianVector &splitPosition1, const bool isForwardInX, const TwoDSlidingFitResult &fitResult1,
        const TwoDSlidingFitResult &fitResult2, const TwoDSlidingFitResult &fitResult3) const;

    /**
     *  @brief  Whether pointing cluster labelled A extends to lowest x positions (as opposed to that labelled B)
     * 
     *  @param  pointingClusterA pointing cluster A
     *  @param  pointingClusterB pointing cluster B
     */
    static bool IsALowestInX(const LArPointingCluster &pointingClusterA, const LArPointingCluster &pointingClusterB);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool            m_majorityRulesMode;                ///< Whether to run in majority rules mode (always split overshoots, always merge undershoots)
    unsigned int    m_nCommonClusters;                  ///< The number of common clusters
    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for use as a key tensor element
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for use as a key tensor element
    float           m_minLongitudinalImpactParameter;   ///< The min longitudinal impact parameter for connecting accompanying clusters
    int             m_nLayersForKinkSearch;             ///< The number of sliding fit layers to step in the kink search
    float           m_additionalXStepForKinkSearch;     ///< An additional (safety) step to tack-on when choosing x sampling points

private:
    /**
     *  @brief  Get modification objects, identifying required splits and merges for clusters
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  modificationList to be populated with modifications
     */
    void GetModifications(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor, ModificationList &modificationList) const;

    /**
     *  @brief  Select elements representing possible components of interest due to overshoots or undershoots in clustering
     * 
     *  @param  eIter iterator to a candidate element
     *  @param  elementList the provided element list
     *  @param  usedClusters the list of used clusters
     *  @param  iteratorList to receive a list of iterators to relevant elements
     */
    void SelectTensorElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
        const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Apply the changes cached in a modification list and update the tensor accordingly
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  modificationList the modification list
     * 
     *  @return whether changes to the tensor have been made
     */
    bool ApplyChanges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ModificationList &modificationList) const;
};

} // namespace lar_content

#endif // #ifndef THREE_D_KINK_BASE_TOOL_H
