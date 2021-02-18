/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewThreeDKinkBaseTool.h
 *
 *  @brief  Header file for the three d kink base tool
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_THREE_D_KINK_BASE_TOOL_H
#define TWO_VIEW_THREE_D_KINK_BASE_TOOL_H 1

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewThreeDKinkBaseTool class
 */
class TwoViewThreeDKinkBaseTool : public TransverseMatrixTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  nCommonClusters the number of common clusters to select
     */
    TwoViewThreeDKinkBaseTool(const unsigned int nCommonClusters);

    /**
     *  @brief  Destructor
     */
    virtual ~TwoViewThreeDKinkBaseTool();

    bool Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);

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
     *  @brief  Whether a provided (iterator to a) matrix element passes the selection cuts for overshoot identification
     *
     *  @param  eIter the iterator to the matrix element
     *  @param  usedClusters the list of used clusters
     */
    virtual bool PassesElementCuts(MatrixType::ElementList::const_iterator eIter, const pandora::ClusterSet &usedClusters) const;

    /**
     *  @brief  Get modification objects for a specific elements of the matrix, identifying required splits and merges for clusters
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  iteratorList list of iterators to relevant matrix elements
     *  @param  modificationList to be populated with modifications
     */
    virtual void GetIteratorListModifications(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList,
        ModificationList &modificationList) const = 0;

    /**
     *  @brief  Get a sampling point in x that is common to sliding linear fit objects in two views
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
        const TwoDSlidingFitResult &fitResult2/*, const TwoDSlidingFitResult &fitResult3*/) const;

    /**
     *  @brief  Whether pointing cluster labelled A extends to lowest x positions (as opposed to that labelled B)
     *
     *  @param  pointingClusterA pointing cluster A
     *  @param  pointingClusterB pointing cluster B
     */
    static bool IsALowestInX(const LArPointingCluster &pointingClusterA, const LArPointingCluster &pointingClusterB);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_nCommonClusters;                  ///< The number of common clusters
    bool            m_majorityRulesMode;                ///< Whether to run in majority rules mode (always split overshoots, always merge undershoots)
    float           m_minXOverlapFraction;           ///< The min x overlap fraction value for particle creation		//FROM CALO MATCHING - clear tracks tool!
    float           m_minMatchingScore;              ///< The min global matching score for particle creation			//FROM CALO MATCHING - clear tracks tool!
    float           m_minLocallyMatchedFraction;     ///< The min locally matched fraction for particle creation		//FROM CALO MATCHING - clear tracks tool!
    //float           m_minMatchedFraction;               ///< The min matched sampling point fraction for use as a key matrix element
    //unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for use as a key matrix element
    float           m_minLongitudinalImpactParameter;   ///< The min longitudinal impact parameter for connecting accompanying clusters
    int             m_nLayersForKinkSearch;             ///< The number of sliding fit layers to step in the kink search
    float           m_additionalXStepForKinkSearch;     ///< An additional (safety) step to tack-on when choosing x sampling points

private:
    /**
     *  @brief  Get modification objects, identifying required splits and merges for clusters
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapMatrix the overlap matrix
     *  @param  modificationList to be populated with modifications
     */
    void GetModifications(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const MatrixType &overlapMatrix, ModificationList &modificationList) const;

    /**
     *  @brief  Select elements representing possible components of interest due to overshoots or undershoots in clustering
     *
     *  @param  eIter iterator to a candidate element
     *  @param  elementList the provided element list
     *  @param  usedClusters the list of used clusters
     *  @param  iteratorList to receive a list of iterators to relevant elements
     */
    void SelectMatrixElements(MatrixType::ElementList::const_iterator eIter, const MatrixType::ElementList &elementList,
        const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Apply the changes cached in a modification list and update the matrix accordingly
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  modificationList the modification list
     *
     *  @return whether changes to the matrix have been made
     */
    bool ApplyChanges(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const ModificationList &modificationList) const;
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_THREE_D_KINK_BASE_TOOL_H
