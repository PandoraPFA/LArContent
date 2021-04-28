/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewThreeDKinkTool.h
 *
 *  @brief  Header file for the two view three d kink base tool
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_THREE_D_KINK_TOOL_H
#define TWO_VIEW_THREE_D_KINK_TOOL_H 1

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewThreeDKinkTool class
 */
class TwoViewThreeDKinkTool : public TransverseMatrixTool
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  nCommonClusters the number of common clusters to select
     */
    TwoViewThreeDKinkTool();

    /**
     *  @brief  Destructor
     */
    ~TwoViewThreeDKinkTool();

    bool Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);

private:
    /**
     *  @brief  Modification class
     */
    class Modification
    {
    public:
        SplitPositionMap m_splitPositionMap;     ///< The split position map
        ClusterMergeMap m_clusterMergeMap;       ///< The cluster merge map
        pandora::ClusterList m_affectedClusters; ///< The list of affected clusters
    };

    typedef std::vector<Modification> ModificationList;

    /**
     *  @brief  Whether a provided (iterator to a) matrix element passes the
     * selection cuts for overshoot identification
     *
     *  @param  eIter the iterator to the matrix element
     *  @param  usedClusters the list of used clusters
     */
    bool PassesElementCuts(MatrixType::ElementList::const_iterator eIter, const pandora::ClusterSet &usedClusters) const;

    /**
     *  @brief  Get a sampling point in x that is common to sliding linear fit
     * objects in two views
     *
     *  @param  splitPosition1 the split position in view 1
     *  @param  isForwardInX whether to work forwards (or backwards) in x
     *  @param  fitResult1 the sliding fit result in view 1
     *  @param  fitResult2 the sliding fit result in view 2
     *
     *  @return the sampling point
     */
    float GetXSamplingPoint(const pandora::CartesianVector &splitPosition1, const bool isForwardInX, const TwoDSlidingFitResult &fitResult1,
        const TwoDSlidingFitResult &fitResult2) const;

    /**
     *  @brief  Whether pointing cluster labelled A extends to lowest x positions
     * (as opposed to that labelled B)
     *
     *  @param  pointingClusterA pointing cluster A
     *  @param  pointingClusterB pointing cluster B
     */
    static bool IsALowestInX(const LArPointingCluster &pointingClusterA, const LArPointingCluster &pointingClusterB);

    /**
     *  @brief  Get modification objects, identifying required splits and merges
     * for clusters
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapMatrix the overlap matrix
     *  @param  modificationList to be populated with modifications
     */
    void GetModifications(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const MatrixType &overlapMatrix, ModificationList &modificationList) const;

    /**
     *  @brief  Apply the changes cached in a modification list and update the
     * matrix accordingly
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  modificationList the modification list
     *
     *  @return whether changes to the matrix have been made
     */
    bool ApplyChanges(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const ModificationList &modificationList) const;

    /**
     *  @brief  Select elements representing possible components of interest due
     * to overshoots or undershoots in clustering
     *
     *  @param  eIter iterator to a candidate element
     *  @param  elementList the provided element list
     *  @param  usedClusters the list of used clusters
     *  @param  iteratorList to receive a list of iterators to relevant elements
     */
    void SelectMatrixElements(MatrixType::ElementList::const_iterator eIter, const MatrixType::ElementList &elementList,
        const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Particle class
     */
    class Particle
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  elementA the matrix element A
         *  @param  elementB the matrix element B
         */
        Particle(const MatrixType::Element &elementA, const MatrixType::Element &elementB);

        const pandora::Cluster *m_pClusterA;      ///< Address of non-shared cluster in element A
        const pandora::Cluster *m_pClusterB;      ///< Address of non-shared cluster in element B
        const pandora::Cluster *m_pCommonCluster; ///< Address of the common cluster
    };

    void GetIteratorListModifications(
        TwoViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const;

    /**
     *  @brief  Whether the provided particle is consistent with being a kink,
     * when examined in three dimensions at the provided split position
     *
     *  @param  pAlgorithm the calling algorithm
     *  @param  particle the particle
     *  @param  splitPosition the candidate split position
     *  @param  isALowestInX whether cluster associated with matrix element a
     * extends to lowest x positions
     *
     *  @return boolean
     */
    bool IsThreeDKink(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const Particle &particle,
        const pandora::CartesianVector &splitPosition, const bool isALowestInX) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minXOverlapFraction;            ///< The min x overlap fraction value for particle creation
    float m_minMatchingScore;               ///< The min global matching score for particle creation
    float m_minLocallyMatchedFraction;      ///< The min locally matched fraction for particle creation
    float m_minLongitudinalImpactParameter; ///< The min longitudinal impact parameter for connecting accompanying clusters
    int m_nLayersForKinkSearch;             ///< The number of sliding fit layers to step in the kink search
    float m_additionalXStepForKinkSearch;   ///< An additional (safety) step to tack-on when choosing x sampling
    bool m_splitMode;                       ///< Whether to run in cluster splitting mode, as opposed to cluster merging mode
    float m_maxTransverseImpactParameter;   ///< The maximum transverse impact parameter for connecting broken clusters
    float m_minImpactParameterCosTheta;     ///< The minimum cos theta (angle between vertex directions) for connecting broken clusters
    float m_cosThetaCutForKinkSearch;       ///< The cos theta cut used for the kink search in three dimensions


};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_THREE_D_KINK_TOOL_H
