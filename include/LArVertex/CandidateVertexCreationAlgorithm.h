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
     *  @brief  Create candidate vertex positions by comparing pairs of cluster end positions
     *
     *  @param  clusterList1 the list of clusters in view 1
     *  @param  clusterList1 the list of clusters in view 2
     */
    void CreateClusterEndPointComparisonVertices(const pandora::ClusterList &clusterList1, const pandora::ClusterList &clusterList2) const;
    
    /**
     *  @brief  Extrapolates 2D clusters, finds where they cross, and matches these crossing points between views to create vertex candidates
     *
     *  @param  clusterListU the clusters in the u view
     *  @param  clusterListV the clusters in the v view
     *  @param  clusterListW the clusters in the w view
     */
    void CreateCrossingVertices(pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW);
    
    /**
     *  @brief  Looks at long 2D clusters to see if a proton and muon cluster have been merged. It looks at a jump in hit energy and matches these jumps between views to crate energy-based vertex candidates.
     *
     *  @param  clusterListU the clusters in the u view
     *  @param  clusterListV the clusters in the v view
     *  @param  clusterListW the clusters in the w view
     */
    void CreateEnergyVertices(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV, const pandora::ClusterList &clusterListW);

    /**
     *  @brief  Select a subset of input clusters (contained in the input list names) for processing in this algorithm
     *
     *  @param  clusterListU the clusters in the u view
     *  @param  clusterListV the clusters in the v view
     *  @param  clusterListW the clusters in the w view
     */
    void SelectClusters(pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW);
    
    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     *
     *  @param  clusterListName the cluster list name
     *  @param  selectedClusterList to receive the selected cluster list
     */
    void SelectClusters(const std::string &clusterListName, pandora::ClusterList &selectedClusterList);

    /**
     *  @brief  Create a candidate vertex position, using an end-point position from one cluster and the fit for a second cluster
     *
     *  @param  position1 an end-point position for the first cluster
     *  @param  hitType1 the hit type of the first cluster
     *  @param  fitResult2 the two dimensional sliding fit result for the second cluster
     */
    void CreateVertex(const pandora::CartesianVector &position1, const pandora::HitType hitType1, const TwoDSlidingFitResult &fitResult2) const;

    /**
     *  @brief  Extrapolates 2D clusters and finds where the clusters cross in 2D
     *
     *  @param  clusterList the list of input clusters
     *  @param  crossingsVector stores the crossing points in 2D
     */
    void Find2DClusterCrossings(const pandora::ClusterList &clusterList, std::vector<pandora::CartesianVector> &crossingsVector);
    
    /**
     *  @brief  The method for extrapolating the clusters in 2D using 2D sliding linear fits
     *
     *  @param  spacePointVector the list containing all cluster 2D hit positions as well as extrapolated positions
     *  @param  pCluster address of the cluster to be extrapolated
     */
    void GetExtrapolatedClusterSpacepoints(std::vector<pandora::CartesianVector> &spacePointVector, const pandora::Cluster *const pCluster);
    
    /**
     *  @brief  The method for finding where the extrapolated clusters cross
     *
     *  @param  spacePointVector1 the list containing all extrapolated cluster 2D hit positions of cluster 1
     *  @param  spacePointVector2 the list containing all extrapolated cluster 2D hit positions of cluster 2
     *  @param  crossingsVector stores the crossing points in 2D
     */
    void FindCrossingsFromSpacepoints(std::vector<pandora::CartesianVector> &spacePointVector1, std::vector<pandora::CartesianVector> &spacePointVector2, std::vector<pandora::CartesianVector> &crossingsVector);
    
    /**
     *  @brief  Create a candidate vertex position, using 2D crossing points in 2 views, by attempting to match these crossing points between the two views
     *
     *  @param  crossingsVector1 the crossing points in 2D in view 1
     *  @param  crossingsVector2 the crossing points in 2D in view 2
     *  @param  hitType1 the view in which the points in crossingsVector1 exist
     *  @param  hitType2 the view in which the points in crossingsVector2 exist
     */
    void CreateMatchedVertices(std::vector<pandora::CartesianVector> &crossingsVector1, std::vector<pandora::CartesianVector> &crossingsVector2, pandora::HitType hitType1, pandora::HitType hitType2) const;
    
    /**
     *  @brief  Creates a Cartesian vector where the X coordinate is the longitudinal distance along a 2D sliding fit of a hit and the Z cooridnate is the hit energy
     *
     *  @param  pCluster address of the cluster for which the vector is to be made
     *  @param  energyAlongRLvector the output vector by reference
     */
    void CreateEnergyAlongRLVector(const pandora::Cluster *const pCluster, std::vector<pandora::CartesianVector> &energyAlongRLvector);
    
    /**
     *  @brief  Draws a plot where the X coordinate is the longitudinal distance along a 2D sliding fit of a hit and the Z cooridnate is the hit energy
     *
     *  @param  energyAlongRLvector the RL vector to be drawn
     *  @param  pCluster address of the cluster from which the vector was made
     * 
     */
    void DrawEnergyVector(std::vector<pandora::CartesianVector> &energyAlongRLvector, const pandora::Cluster* pCluster);
    
    /**
     *  @brief  Filters out hits with outlying hit energies from an RL, hit energy Cartesian vector
     *
     *  @param  unfilteredEnergyVector the unfiltered input RL, hit energy vector
     *  @param  filteredEnergyVector the filtered output vector by reference
     * 
     */
    void FilterEnergyVector(const std::vector<pandora::CartesianVector> &unfilteredEnergyVector, std::vector<pandora::CartesianVector> &filteredEnergyVector);
    
    /**
     *  @brief  Bins an input RL, hit energy vector into bins 1cm wide, with X coordinate the average RL cooridnate of the bin, and as Z coordinate the average hit energy of the bin
     *
     *  @param  energyAlongRLvector the unbinned input RL, hit energy vector
     *  @param  binnedEnergyAlongRLvector the binned output vector by reference
     * 
     */
    void BinEnergyRLVector(const std::vector<pandora::CartesianVector> &energyAlongRLvector, std::vector<pandora::CartesianVector> &binnedEnergyAlongRLvector);
    
    /**
     *  @brief  Attempts to find the longitudinal coordinate along a 2D sliding lionear fit of an energy spike
     *
     *  @param  energyAlongRLvector the input RL, hit energy vector
     *  @param  spikeRLvector the vector containing all the RL coordinates of the energy jumps
     * 
     */
    void FindEnergySpike(const std::vector<pandora::CartesianVector> &energyAlongRLvector, std::vector<float> &spikeRLvector);
    
    /**
     *  @brief  Attempts to find the bin in which there is an energy jump in a binned RL, hit energy vector
     *
     *  @param  energyAlongRLvector the input binned RL, hit energy vector
     *  @param  spikeRLvector the vector containing all the average RL coordinates of the bins in which there are energy jumps
     * 
     */
    void FindBinWithSpike(const std::vector<pandora::CartesianVector> &binnedEnergyAlongRLvector, std::vector<float> &binRLvector);
    
    /**
     *  @brief  Converts an average bin RL coordinate to the RL coordinate of the hit that best matches the energy jump
     *
     *  @param  binRLvector the vector containing all the average RL coordinates of the bins in which there are energy jumps
     *  @param  spikeRLvector the vector containing all the RL coordinates of the energy jumps
     *  @param  energyAlongRLvector the original energy along RL vector
     * 
     */
    void ConvertBinRLToSpikeRL(const std::vector<float> &binRLvector, std::vector<float> &spikeRLvector, const std::vector<pandora::CartesianVector> &energyAlongRLvector);
    
    /**
     *  @brief  Converts an RL coordinate along a 2D sliding linear fit to a 2D hit position in the corresponding cluster
     *
     *  @param  filteredEnergyAlongRLvector the RL, hit energy vector that was used to find the energy jump
     *  @param  caloHitList the caloHitList of the cluster used to form the RL, hit energy vector
     *  @param  hitPosition the output hit position in 2D by reference
     * 
     */
    void ConvertRLtoCaloHit(const float &spikeRL, std::vector<pandora::CartesianVector> &filteredEnergyAlongRLvector, const pandora::CaloHitList &caloHitList, pandora::CartesianVector &hitPosition);
    
    /**
     *  @brief  When a hit corresponding to an energy jump has been found, this method find corresponding hits in the other two views by matching the hits very narrowly by their X coordinates (drift times)
     *
     *  @param  clusterList all the clusters in the target view
     *  @param  energySpikeVector the vector containing the 2D hit positions of the energy spikes
     *  @param  matchedHits the output matched hits by reference
     * 
     */
    void FindMatchingHitsInDifferentView(const pandora::ClusterList &clusterList, pandora::CartesianVector &energySpikeVector, std::vector<pandora::CartesianVector> &matchedHits);
    
    /**
     *  @brief  The final step that creates energy vertex candidates using the output from the previous methods
     *
     *  @param  energySpikeRLvector vector containing all the RL positions of the energy spikes
     *  @param  filteredEnergyAlongRLvector the filtered RL, hit energy vector used to find the energy spikes
     *  @param  caloHitList the calo hits corresponding to the cluster used to find the energy spikes
     *  @param  clusterList1 all the clusters in view 1
     *  @param  clusterList2 all the clusters in view 2
     *  @param  clusterList3 all the clusters in view 3
     * 
     */
    void CreateVerticesFromSpikes(const std::vector<float>energySpikeRLvector, std::vector<pandora::CartesianVector> filteredEnergyAlongRLvector, pandora::CaloHitList &caloHitList, const pandora::ClusterList &clusterList1, const pandora::ClusterList &clusterList2, const pandora::ClusterList &clusterList3);
    
    /**
     *  @brief  Creates a 2D sliding fit of a cluster and stores it for later use
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
    
    /**
     *  @brief  Sorts 2D positions in two Cartesian vectors by their Z coordinates
     * 
     *  @param  vector1 the first vector to be sorted
     *  @param  vector2 the second vector to be sorted
     */
    static bool SortSpacePointsByZ(pandora::CartesianVector &vector1, pandora::CartesianVector &vector2);
    
    /**
     *  @brief  Sorts 2D positions in two Cartesian vectors by their X coordinates
     * 
     *  @param  vector1 the first vector to be sorted
     *  @param  vector2 the second vector to be sorted
     */
    static bool SortEnergyVectorByRL(pandora::CartesianVector &vector1, pandora::CartesianVector &vector2);
    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string                 m_topologyVertexListName;       ///< The name under which to save the output topology-based vertex list
    std::string                 m_energyVertexListName;         ///< The name under which to save the output energy-based vertex list

    bool                        m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list

    unsigned int                m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap     m_slidingFitResultMap;          ///< The sliding fit result map

    unsigned int                m_minClusterCaloHits;           ///< The min number of hits in base cluster selection method
    float                       m_minClusterLengthSquared;      ///< The min length (squared) in base cluster selection method
    float                       m_maxClusterXDiscrepancy;       ///< The max cluster end-point discrepancy
    float                       m_chiSquaredCut;                ///< The chi squared cut (accept only values below the cut value)
    
    bool                        m_enableCrossingCandidates;     ///< Whether to create crossing vertex candidates
    bool                        m_enableEnergyCandidates;       ///< Whether to create energy-based candidates
    bool                        m_energyPlot;                   ///< Whether to draw 
    
    unsigned int                m_minCrossingClusterSize;       ///< The minimum number of hits a cluster needs to have to be considered in the crossing vertex procedure
    
    float                       m_extrapolationLength;          ///< How far in cm to extrapolate a cluster in order to look for 2D crossing points
    float                       m_extrapolationStepSize;        ///< The distance between extrapolation points in cm
    float                       m_minClusterCrossingApproach;   ///< The minimum impact parameter between two extrapolated clusters to be marked as crossing
    float                       m_postCrossingSkipDistance;     ///< How far to skip along an exrapolated cluster after finding a crossing point (for speed)
    
    unsigned int                m_minEnergyVertexClusterSize;   ///< The minimum number of hits a cluster needs to have to be considered in the crossing vertex procedure
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CandidateVertexCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CandidateVertexCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
