/**
 *  @file   LArContent/include/LArVertex/CandidateVertexCreationAlgorithm.h
 * 
 *  @brief  Header file for the candidate vertex creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
#define LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H 1

#include "LArObjects/LArTwoDSlidingFitResult.h"
#include "LArObjects/LArThreeDSlidingFitResult.h"

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
     *  @param  clusterListU to receive the selected clusters in the u view
     *  @param  clusterListV to receive the selected clusters in the v view
     *  @param  clusterListW to receive the selected clusters in the w view
     */
    void SelectClusters(pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW);

    /**
     *  @brief  Create candidate vertex positions by comparing pairs of cluster end positions
     *
     *  @param  clusterList1 the list of clusters in view 1
     *  @param  clusterList1 the list of clusters in view 2
     */
    void ClusterEndPointComparison(const pandora::ClusterList &clusterList1, const pandora::ClusterList &clusterList2) const;

    /**
     *  @brief  Create a candidate vertex position, using an end-point position from one cluster and the fit for a second cluster
     *
     *  @param  position1 an end-point position for the first cluster
     *  @param  hitType1 the hit type of the first cluster
     *  @param  fitResult2 the two dimensional sliding fit result for the second cluster
     */
    void CreateVertex(const pandora::CartesianVector &position1, const pandora::HitType hitType1, const TwoDSlidingFitResult &fitResult2) const;

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     *
     *  @param  clusterListName the cluster list name
     *  @param  selectedClusterList to receive the selected cluster list
     */
    void SelectClusters(const std::string &clusterListName, pandora::ClusterList &selectedClusterList);

    /**
     *  @brief  Add a new sliding fit result, for the specified cluster, to the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
     
    /**
     *  @brief  Create candidate vertex positions by comparing pairs of cluster end positions
     *
     *  @param  clusterList1 the list of clusters in view 1
     *  @param  clusterList1 the list of clusters in view 2
     */
    void Find2DClusterCrossings(const pandora::ClusterList &clusterList, std::vector<pandora::CartesianVector> &crossingsVector) const;
    
    /**
     *  @brief  Create a candidate vertex position, using an end-point position from one cluster and the fit for a second cluster
     *
     *  @param  position1 an end-point position for the first cluster
     *  @param  hitType1 the hit type of the first cluster
     *  @param  fitResult2 the two dimensional sliding fit result for the second cluster
     */
    void CreateEnergySpikeVertices(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV, const pandora::ClusterList &clusterListW);
    

    /**
     *  @brief  Create a candidate vertex position, using an end-point position from one cluster and the fit for a second cluster
     *
     *  @param  position1 an end-point position for the first cluster
     *  @param  hitType1 the hit type of the first cluster
     *  @param  fitResult2 the two dimensional sliding fit result for the second cluster
     */
    void CreateMatchedVertices(std::vector<pandora::CartesianVector> &crossingsVector1, std::vector<pandora::CartesianVector> &crossingsVector2,  pandora::HitType hitType1, pandora::HitType hitType2) const;
     
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
    
    static bool SortSpacePointsByZ(pandora::CartesianVector &vector1, pandora::CartesianVector &vector2);
    
    static bool SortEnergyVectorByRL(pandora::CartesianVector &position1, pandora::CartesianVector &position2);
    
    void CreateEnergyAlongRLVector(const TwoDSlidingFitResult &slidingFitResult, const pandora::CaloHitList &caloHitList, std::vector<pandora::CartesianVector> &energyAlongRLvector);
    
    void DrawEnergyVector(std::vector<pandora::CartesianVector> &energyAlongRLvector, const pandora::Cluster* pCluster);
    
    void FilterEnergyVector(const std::vector<pandora::CartesianVector> &unfilteredEnergyVector, std::vector<pandora::CartesianVector> &filteredEnergyVector);
    
    void FindMatchingHitsInDifferentView(const pandora::ClusterList &clusterList, pandora::CartesianVector &energySpikeVector, std::vector<pandora::CartesianVector> &matchedHits);
    
    void FindEnergySpike(std::vector<pandora::CartesianVector> &energyAlongRLvector, std::vector<float> &spikeRLvector, bool &foundSplit) const;
    
    void ConvertRLtoCaloHit(const float &spikeRL, const TwoDSlidingFitResult &slidingFitResult, const pandora::CaloHitList &caloHitList, pandora::CartesianVector &hitPosition);
    
    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string                 m_outputVertexListName;         ///< The name under which to save the output vertex list
    bool                        m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list

    unsigned int                m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap     m_slidingFitResultMap;          ///< The sliding fit result map

    unsigned int                m_minClusterCaloHits;           ///< The min number of hits in base cluster selection method
    float                       m_minClusterLengthSquared;      ///< The min length (squared) in base cluster selection method
    float                       m_maxClusterXDiscrepancy;       ///< The max cluster end-point discrepancy
    float                       m_chiSquaredCut;                ///< The chi squared cut (accept only values below the cut value)
    
    bool                        m_enableCrossingCandidates;
    bool                        m_enableEnergyCandidates;
    
    bool                        m_strictMatching;
    bool                        m_energyPlot;
    
    float                       m_maxScatterRms;
    float                       m_maxScatterCosTheta;
    float                       m_maxSlidingCosTheta;
    
    float                       m_energyDifferenceThreshold;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CandidateVertexCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CandidateVertexCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
