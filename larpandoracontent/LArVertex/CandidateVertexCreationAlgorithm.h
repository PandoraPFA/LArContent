/**
 *  @file   larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.h
 * 
 *  @brief  Header file for the candidate vertex creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
#define LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CandidateVertexCreationAlgorithm::Algorithm class
 */
class CandidateVertexCreationAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Default constructor
     */
    CandidateVertexCreationAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Select a subset of input clusters (contained in the input list names) for processing in this algorithm
     *
     *  @param  clusterVectorU to receive the selected clusters in the u view
     *  @param  clusterVectorV to receive the selected clusters in the v view
     *  @param  clusterVectorW to receive the selected clusters in the w view
     */
    void SelectClusters(pandora::ClusterVector &clusterVectorU, pandora::ClusterVector &clusterVectorV, pandora::ClusterVector &clusterVectorW);

    /**
     *  @brief  Create candidate vertex positions by comparing pairs of cluster end positions
     *
     *  @param  clusterVector1 the clusters in view 1
     *  @param  clusterVector1 the clusters in view 2
     */
    void ClusterEndPointComparison(const pandora::ClusterVector &clusterVector1, const pandora::ClusterVector &clusterVector2) const;

    /**
     *  @brief  Create a candidate vertex position, using an end-point position from one cluster and the fit for a second cluster
     *
     *  @param  position1 an end-point position for the first cluster
     *  @param  hitType1 the hit type of the first cluster
     *  @param  fitResult2 the two dimensional sliding fit result for the second cluster
     */
    void CreateVertex(const pandora::CartesianVector &position1, const pandora::HitType hitType1, const TwoDSlidingFitResult &fitResult2) const;

    /**
     *  @brief  Add a new sliding fit result, for the specified cluster, to the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get a sliding fit result from the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    const TwoDSlidingFitResult &GetCachedSlidingFitResult(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Clear relevant algorithm member variables between events
     */
    void TidyUp();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector       m_inputClusterListNames;        ///< The list of cluster list names
    std::string                 m_outputVertexListName;         ///< The name under which to save the output vertex list
    bool                        m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list

    unsigned int                m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap     m_slidingFitResultMap;          ///< The sliding fit result map

    unsigned int                m_minClusterCaloHits;           ///< The min number of hits in base cluster selection method
    float                       m_minClusterLengthSquared;      ///< The min length (squared) in base cluster selection method
    float                       m_maxClusterXDiscrepancy;       ///< The max cluster end-point discrepancy
    float                       m_chiSquaredCut;                ///< The chi squared cut (accept only values below the cut value)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CandidateVertexCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CandidateVertexCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
