/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/DLThreeDClusterSplittingAlgorithm.h
 *
 *  @brief  Header file for the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_THREE_D_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_DL_THREE_D_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_dl_content
{
/**
 *  @brief  DLThreeDClusterSplittingAlgorithm class
 */
class DLThreeDClusterSplittingAlgorithm : public pandora::Algorithm
{
/**
 *  @brief  MCContaminant class
 */
class MCContaminant
{
 public:
    /**
     *  @brief  Constructor
     */
    MCContaminant(const pandora::MCParticle *const pMCParticle, const pandora::CartesianVector &contStartPosition, 
                  const pandora::CartesianVector &contEndPosition, const pandora::CartesianVector &startPosition, 
                  const pandora::CartesianVector &endPosition, const pandora::CartesianVector &startDirection,
                  const pandora::CartesianVector &endDirection);

    const pandora::MCParticle *m_pMCParticle;     ///< The characterised MCParticle
    pandora::CartesianVector m_contStartPosition; ///< The start position within the cluster calculated from hits the MCParticle contributes to
    pandora::CartesianVector m_contEndPosition;   ///< The end position within the cluster calculated from hits the MCParticle contributes to
    pandora::CartesianVector m_startPosition;     ///< The start position within the cluster calculated from hits the MCParticle 'owns'
    pandora::CartesianVector m_endPosition;       ///< The end position within the cluster calculated from hits the MCParticle 'owns'
    pandora::CartesianVector m_startDirection;    ///< The direction at the start position
    pandora::CartesianVector m_endDirection;      ///< The direction at the end position
};
/**
 *  @brief  ClusterHit class
 */
class ClusterHit
{
public:
    /**
     *  @brief  Constructor
     */
    ClusterHit(const pandora::CaloHit *const pCaloHit, const float l, const float t);

    /**
     *  @brief  ClusterHit == operator
     *
     *  @param  rhs the ClusterHit to compare
     */
    bool operator==(const ClusterHit &rhs) const;

    const pandora::CaloHit *m_pHit; ///< The characterised CaloHit
    float m_l;                      ///< The longitudinal position of the CaloHit from the cluster fit
    float m_t;                      ///< The transverse position of the CaloHit from the cluster fit
};
/**
 *  @brief  ClusterHit class
 */
class Feature
{
public:
    /**
     *  @brief  Constructor
     */
    Feature(const float mean, const float std, const pandora::FloatVector &sequence);

    /**
     *  @brief  Normalise the feature vector
     */
    void Normalise();

    /**
     *  @brief  Smooth the feature vector
     */
    void Smooth();

    /**
     *  @brief  Obtain the feature vector
     *
     *  @return FloatVector a const reference to the feature vector
     */
    const pandora::FloatVector* GetSequence() const;

    float m_mean;                    ///< The feature mean
    float m_std;                     ///< The feature standard deviation
    pandora::FloatVector m_sequence; ///< The feature vector
};
public:
    /**
     *  @brief  Default constructor
     */
    DLThreeDClusterSplittingAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~DLThreeDClusterSplittingAlgorithm();

    pandora::StatusCode Run();

private:
    typedef std::vector<MCContaminant> MCContaminantVector;
    typedef std::vector<ClusterHit> ClusterPath;
    typedef lar_content::KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;
    typedef lar_content::KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;
    typedef std::unordered_map<const pandora::MCParticle *, pandora::CaloHitList> MCParticleToHitListMap;
    typedef std::map<std::string, Feature> Features;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Fill algorithm lists
     *
     *  @return StatusCode whether all required lists have been filled
     */
    pandora::StatusCode GetLists();

    /**
     *  @brief  Identify true or predicted split positions 
     *
     *  @param  pCluster the input cluster
     *  @param  clusterFit the TwoDSlidingFitResult of the cluster 
     *  @param  kdTree kdTree of view hits
     *  @param  splitPositions to store found split positions
     */
    void ProcessCluster(const pandora::Cluster *const pCluster, const lar_content::TwoDSlidingFitResult &clusterFit, 
         HitKDTree2D &kdTree, pandora::CartesianPointVector &splitPositions);

    /**
     *  @brief  Identify the 'sequence' hits of the cluster i.e those used to fill the feature vectors
     *
     *  @param  pCluster the input cluster
     *  @param  clusterFit the TwoDSlidingFitResult of the cluster 
     *  @param  clusterPath the output hit sequence 
     */
    void FindPath(const pandora::Cluster *const pCluster, const lar_content::TwoDSlidingFitResult &clusterFit, ClusterPath &clusterPath) const;

    /**
     *  @brief  Initialise the map of feature vectors
     *
     *  @param  features the input feature map
     */
    void InitialiseFeatures(Features &features) const;

    /**
     *  @brief  From the clusterPath, fill the feature vectors
     *
     *  @param  clusterPath the input hit sequence 
     *  @param  clusterFit the input TwoDSlidingFitResult of the cluster 
     *  @param  features the input feature map
     *  @param  kdTree kdTree of view hits
     */
    void FillFeatures(const ClusterPath &clusterPath, const lar_content::TwoDSlidingFitResult &clusterFit, Features &features, HitKDTree2D &kdTree) const;

    /**
     *  @brief  Initialise the KalmanFilter2D, and seed its fit
     *
     *  @param  clusterPath the input hit sequence 
     */
    lar_content::KalmanFilter2D InitialiseKalmanFilter(const ClusterPath &clusterPath) const;

    /**
     *  @brief  Project the secondary vertices into the considered view and rule out those that do not lie within the cluster
     *
     *  @param  clusterFit the TwoDSlidingFitResult of the cluster 
     *  @param  viewSecVtx the container to store the filtered 2D secondary vertices
     */
    void GetFilteredViewSecVertices(const lar_content::TwoDSlidingFitResult &clusterFit, pandora::CartesianPointVector &viewSecVtx) const;

    /**
     *  @brief  Get the separation between CaloHit and the closest secondary vertex
     *
     *  @param  pCaloHit the input CaloHit
     *  @param  viewSecVtx the filtered 2D secondary vertices
     *
     *  @return the separation between CaloHit and the closest secondary vertex
     */
    float GetDistanceToSecVertex(const pandora::CaloHit *const pCaloHit, const pandora::CartesianPointVector &viewSecVtx) const;

    /**
     *  @brief  Get the separation between CaloHit and the closest non-clusterPath hit
     *
     *  @param  kdTree kdTree of view hits
     *  @param  pCaloHit the input CaloHit
     *  @param  clusterPathHits the hits of the cluster path
     *
     *  @return the separation between CaloHit and the closest non-clusterPath hit
     */
    float GetDistanceToEventHit(HitKDTree2D &kdTree, const pandora::CaloHit *const pCaloHit, 
        const pandora::CaloHitList &clusterPathHits) const;

    /**
     *  @brief  Get the separation between a ClusterHit and the next in the ClusterPath
     *
     *  @param  currentHit the considered ClusterHit
     *  @param  nextHit the subsequent ClusterHit
     *
     *  @return the separation between a ClusterHit and the next in the ClusterPath
     */
    float GetDistanceToClusterHit(const ClusterHit &currentHit, const ClusterHit &nextHit) const;

    /**
     *  @brief  Determine whether two hits are in the same TPC volume
     *
     *  @param  pPrevHit the first CaloHit
     *  @param  pCurrentHit the second CaloHit
     *
     *  @return whether two hits are in the same TPC volume
     */
    float IsInSameVolume(const pandora::CaloHit *const pPrevHit, const pandora::CaloHit *const pCurrentHit) const;

    /**
     *  @brief  Identify predicted split positions in terms of their sequence position
     *
     *  @param  features the input feature map
     *  @param  splitIndices the container to store found split indices
     */
    void GetPredictedSplitIndices(const Features &features, pandora::IntVector &splitIndices);

    /**
     *  @brief  Run the window & split position models to identify predicted split indices and their scores
     *
     *  @param  features the input feature map
     *  @param  windowStart the sequence index of window starts
     *  @param  splitIndices the container to store found split indices
     *  @param  splitScores the container to store associated split scores
     */
    void GetPredictedSplitIndices(const Features &features, const pandora::IntVector &windowStart, 
        pandora::IntVector &splitIndices, pandora::FloatVector &splitScores);

    /**
     *  @brief  Reduce runs of subsequence split indices to a single index
     *
     *  @param  splitIndices the container of found split indices
     *  @param  splitScores the container of associated split scores
     */
    void FilterModelOutput(const pandora::FloatVector &splitScores, pandora::IntVector &splitIndices) const;

    /**
     *  @brief  Filter predicted split positions by model unrelated cuts
     *
     *  @param  clusterPath the hit sequence 
     *  @param  hitType the type of 2D view
     *  @param  features the feature map
     *  @param  splitIndices the idenitfied split sequence indices
     */
    void FilterSplitIndices(const ClusterPath &clusterPath, const pandora::HitType hitType, const Features &features, 
        pandora::IntVector &splitIndices) const;

    /**
     *  @brief  Consider the split positions of all views, and accept those that form a consistent 3D position
     *
     *  @param  splitPosU split positions in the U view
     *  @param  splitPosV split positions in the V view
     *  @param  splitPosW split positions in the W view
     *  @param  usedU accepted split positions in the U view
     *  @param  usedV accepted split positions in the V view
     *  @param  usedW accepted split positions in the W view
     */
    void ThreeViewMatching(const pandora::CartesianPointVector &splitPosU, const pandora::CartesianPointVector &splitPosV, 
        const pandora::CartesianPointVector &splitPosW, pandora::IntVector &usedU, pandora::IntVector &usedV, pandora::IntVector &usedW);

    /**
     *  @brief  Consider the split positions of two views, and accept those that map onto a hit in the third view
     *
     *  @param  splitPos1 split positions in view 1
     *  @param  splitPos2 split positions in view 2
     *  @param  hitType1 the view type of view 1
     *  @param  hitType2 the view type of view 2
     *  @param  kdTree kdTree of the hits in view 3
     *  @param  used1 accepted split positions in view 1
     *  @param  used2 accepted split positions in view 2
     */
    void TwoViewMatching(const pandora::CartesianPointVector &splitPos1, const pandora::CartesianPointVector &splitPos2, 
        const pandora::HitType hitType1, const pandora::HitType hitType2, HitKDTree2D &kdTree, pandora::IntVector &used1, pandora::IntVector &used2);

    /**
     *  @brief  Split the hits of the cluster into groups determined by the found split indices 
     *
     *  @param  clusterFit the TwoDSlidingFitResult of the considered cluster 
     *  @param  splitPositions identified split positions
     *
     *  @return the groups of CaloHits from which to form the new clusters
     */
    std::vector<pandora::CaloHitList> DivideCaloHits(const lar_content::TwoDSlidingFitResult &clusterFit, 
        const pandora::CartesianPointVector &splitPositions) const;

    /**
     *  @brief  Split the cluster creating new clusters
     *
     *  @param  pCluster the cluster to split
     *  @param  splitClusterHits the groups of CaloHits from which to form the new cluster
     */
    void SplitCluster(const pandora::Cluster *const pCluster, const std::vector<pandora::CaloHitList> &splitClusterHits) const;

//------------------------------------------------------------------------------------------------------------------------------------------
// Training Functions
//------------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Identify the MCParticles that dominate and contribute to the energy of the cluster hits
     *
     *  @param  pCluster the cluster to consider
     *  @param  mainHitListMap the map to fill with dominant contributions
     *  @param  contHitListMap the map to fill with all 'significant' contributions
     */
    void FillHitListMaps(const pandora::Cluster *const pCluster, MCParticleToHitListMap &mainHitListMap, MCParticleToHitListMap &contHitListMap);

    /**
     *  @brief  Find the MCContaminants of the considered cluster
     *
     *  @param  mainHitListMap the map of dominant contributions
     *  @param  contHitListMap the map of 'significant' contributions
     *  @param  mcContaminants the container to store found MCContaminants
     */
    void FindContaminants(const MCParticleToHitListMap &mainHitListMap, const MCParticleToHitListMap &contHitListMap,
        MCContaminantVector &mcContaminants);

    /**
     *  @brief  Create from an input contaminant MCParticle a MCContaminant, 
     *          a parameterisation of a MCParticle that contaminants the cluster
     *
     *  @param  pMCContaminant the input MCParticle
     *  @param  mainHitListMap the map of dominant contributions
     *  @param  contHitListMap the map of 'significant' contributions
     *  @param  mcContaminants the container to store MCContaminants
     */
    void BuildMCContaminant(const pandora::MCParticle *const pMCContaminant, const pandora::CaloHitList &mainHitList, 
        const pandora::CaloHitList &contHitList, MCContaminantVector &mcContaminants);

    /**
     *  @brief  Identify true split positions within the considered cluster
     *
     *  @param  clusterFit the TwoDSlidingFitResult of the considered cluster 
     *  @param  mcContaminants the container to of MCContaminants
     *  @param  trueSplitPositions the container to store found true split positions
     */
    void FindTrueSplitPositions(const lar_content::TwoDSlidingFitResult &clusterFit, const MCContaminantVector &mcContaminants,
        pandora::CartesianPointVector &trueSplitPositions);

    /**
     *  @brief  Fill the training tree
     *
     *  @param  splitPositions the found true split positions
     *  @param  mainHitListMap the map of dominant contributions
     *  @param  clusterFit the TwoDSlidingFitResult of the considered cluster 
     *  @param  features the input feature map
     */
    void FillTree(const pandora::CartesianPointVector &splitPositions, const MCParticleToHitListMap &mainHitListMap, 
        const lar_content::TwoDSlidingFitResult &clusterFit, const Features &features);

    std::string m_caloHitListNameU;   ///< The name of the U view CaloHitList
    std::string m_caloHitListNameV;   ///< The name of the V view CaloHitList
    std::string m_caloHitListNameW;   ///< The name of the W view CaloHitList
    std::string m_clusterListNameU;   ///< The name of the U view ClusterList
    std::string m_clusterListNameV;   ///< The name of the V view ClusterList
    std::string m_clusterListNameW;   ///< The name of the W view ClusterList
    std::string m_nuVertexListName;   ///< The name of the neutrino VertexList
    std::string m_secVertexListName;  ///< The name of the secondary VertexList
    std::string m_mcParticleListName; ///< The name of the MCParticleList
    const pandora::CaloHitList *m_pCaloHitListU;      ///< The U view CaloHitList
    const pandora::CaloHitList *m_pCaloHitListV;      ///< The V view CaloHitList
    const pandora::CaloHitList *m_pCaloHitListW;      ///< The W view CaloHitList
    const pandora::ClusterList *m_pClusterListU;      ///< The U view ClusterList
    const pandora::ClusterList *m_pClusterListV;      ///< The V view ClusterList
    const pandora::ClusterList *m_pClusterListW;      ///< The W view ClusterList
    const pandora::VertexList *m_pNuVertexList;       ///< The neutrino VertexList
    const pandora::VertexList *m_pSecVertexList;      ///< The secondary VertexList
    const pandora::MCParticleList *m_pMCParticleList; ///< The MCParticleList
    bool m_trainingMode;                   ///< Whether to run the algorithm in training mode
    std::string m_fileName;                ///< The name of the output training file
    std::string m_treeName;                ///< The name of the output training tree
    unsigned int m_minClusterHits;         ///< The min number of hits of a considered cluster
    int m_slidingWindow;                   ///< The window length used in TwoDSlidingFits
    float m_lBinSize;                      ///< The longitudinal bin size used to obtain the clusterPath
    float m_kalmanDelta;                   ///< the time step for the kalman filter
    float m_kalmanProcessVarCoeff;         ///< the process variance coefficient for the Kalman filter 
    float m_kalmanMeasurementVarCoeff;     ///< the measurement variance coefficient for the Kalman filter
    float m_searchRegion1D;                ///< The width/length of the region in which to search for event hits
    unsigned int m_windowLength;           ///< The length of the window considered by the window/split point models
    float m_bkgThreshold;                  ///< The maximum background 'probability' of a contaminated window
    float m_isContaminatedThreshold;       ///< The threshold 'probability' of a contaminated window
    float m_isSplitThreshold;              ///< The threshold 'probability' of a predicted split position
    float m_minNuVertexSep;                ///< The minimum separation of a split position from the nu vertex
    float m_gapSearchRegion;               ///< The points either side of the split position in which to search for a gap
    float m_threeViewMatchMaxXSpan;        ///< the max x-span of the split positions in a three-view match
    float m_threeViewMatchMaxChi2;         ///< the max chi2 of real-proj split positions obtained from a triplet match
    float m_twoViewMatchMaxXSpan;          ///< the max x-span of the split positions in a two-view match
    float m_twoViewMatchSearchRegion;      ///< the width of the search region used to find hits in a projected view
    float m_contFraction;                  ///< The min hit charge fraction of a 'contributing' MCParticle
    float m_mainFraction;                  ///< The min hit charge fraction of a 'main' MCParticle
    unsigned int m_minTargetMCHits;        ///< The min number of MC hits of a 'significant' MC contaminant
    float m_endpointBuffer;                ///< The min distance to the min hit charge fraction of a 'contributing' MCParticle
    float m_collinearMaxSep;               ///< The max separation of the endpoints of collinear MCParticles
    float m_collinearMaxOpeningAngle;      ///< The max opening angle between collinear MCParticles
    unsigned int m_reducedMinTargetMCHits; ///< The min number of MC hits of a 'noted' MC contaminant
    std::string m_windowModelName;           ///< The name of the window model weight file
    LArDLHelper::TorchModel m_windowModel;   ///< The window model
    std::string m_splitPosModelName;         ///< The name of the split position model weight file
    LArDLHelper::TorchModel m_splitPosModel; ///< The split position model
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline DLThreeDClusterSplittingAlgorithm::MCContaminant::MCContaminant(const pandora::MCParticle *const pMCParticle, 
    const pandora::CartesianVector &contStartPosition, const pandora::CartesianVector &contEndPosition, 
    const pandora::CartesianVector &startPosition, const pandora::CartesianVector &endPosition, 
    const pandora::CartesianVector &startDirection, const pandora::CartesianVector &endDirection) :
        m_pMCParticle(pMCParticle),
        m_contStartPosition(contStartPosition),
        m_contEndPosition(contEndPosition),
        m_startPosition(startPosition),
        m_endPosition(endPosition),
        m_startDirection(startDirection),
        m_endDirection(endDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DLThreeDClusterSplittingAlgorithm::ClusterHit::ClusterHit(const pandora::CaloHit *const pCaloHit, const float l, const float t) :
    m_pHit(pCaloHit),
    m_l(l),
    m_t(t)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool DLThreeDClusterSplittingAlgorithm::ClusterHit::operator==(const ClusterHit &rhs) const
{
    return this->m_l == rhs.m_l;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DLThreeDClusterSplittingAlgorithm::Feature::Feature(const float mean, const float std, const pandora::FloatVector &sequence) :
    m_mean(mean),
    m_std(std),
    m_sequence(sequence)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DLThreeDClusterSplittingAlgorithm::Feature::Normalise()
{
    for (float &val : m_sequence)
        val = (val - m_mean) / m_std;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DLThreeDClusterSplittingAlgorithm::Feature::Smooth()
{
    pandora::FloatVector temp(m_sequence);
    
    for (unsigned int iEntry = 0; iEntry < m_sequence.size(); ++iEntry)
    {
        float total(temp.at(iEntry));
        int nEntries(1);

        for (unsigned int iWindow = 1; iWindow <= 4; ++iWindow)
        {
            const int belowIndex(iEntry - iWindow);
            const unsigned int aboveIndex(iEntry + iWindow);
            
            if (belowIndex >= 0)
            {
                total += temp.at(belowIndex);
                ++nEntries;
            }
                
            if (aboveIndex < temp.size())
            {
                total += temp.at(aboveIndex);
                ++nEntries;
            }
        }

        m_sequence[iEntry] = (total / nEntries);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::FloatVector* DLThreeDClusterSplittingAlgorithm::Feature::GetSequence() const
{
    return &m_sequence;
}

} // namespace lar_content

#endif // #ifndef LAR_DL_THREE_D_CLUSTER_SPLITTING_ALGORITHM_H
