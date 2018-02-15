/**
 *  @file   larpandoracontent/LArVertex/SvmVertexSelectionAlgorithm.h
 *
 *  @brief  Header file for the svm vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SVM_VERTEX_SELECTION_ALGORITHM_H
#define LAR_SVM_VERTEX_SELECTION_ALGORITHM_H 1

#include "Api/PandoraContentApi.h"

#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

#include <random>

namespace lar_content
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SvmVertexSelectionAlgorithm class
 */
class SvmVertexSelectionAlgorithm : public VertexSelectionBaseAlgorithm
{
public:
    /**
     *  @brief Vertex feature info class
     */
    class VertexFeatureInfo
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  beamDeweighting the beam deweighting feature
         *  @param  rPhiFeature the r/phi feature
         *  @param  energyKick the energy kick feature
         *  @param  localAsymmetry the local asymmetry feature
         *  @param  globalAsymmetry the global asymmetry feature
         *  @param  showerAsymmetry the shower asymmetry feature
         */
        VertexFeatureInfo(const float beamDeweighting, const float rPhiFeature, const float energyKick, const float localAsymmetry,
                          const float globalAsymmetry, const float showerAsymmetry);

        float    m_beamDeweighting;    ///< The beam deweighting feature
        float    m_rPhiFeature;        ///< The r/phi feature
        float    m_energyKick;         ///< The energy kick feature
        float    m_localAsymmetry;     ///< The local asymmetry feature
        float    m_globalAsymmetry;    ///< The global asymmetry feature
        float    m_showerAsymmetry;    ///< The shower asymmetry feature
    };

    typedef std::map<const pandora::Vertex *const, VertexFeatureInfo> VertexFeatureInfoMap;

    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief Event feature info class
     */
    class EventFeatureInfo
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  eventShoweryness the event showeryness feature
         *  @param  eventEnergy the energy of the event
         *  @param  eventVolume the volume of the event
         *  @param  longitudinality the longitudinality of the event
         *  @param  nHits the number of hits in the event
         *  @param  nClusters the number of clusters in the event
         *  @param  nCandidates the total number of vertex candidates
         */
        EventFeatureInfo(const float eventShoweryness, const float eventEnergy, const float eventVolume, const float longitudinality,
            const unsigned int nHits, const unsigned int nClusters, const unsigned int nCandidates);

        float           m_eventShoweryness;        ///< The event showeryness feature
        float           m_eventEnergy;             ///< The event energy
        float           m_eventVolume;             ///< The volume of the event
        float           m_longitudinality;         ///< The longitudinality of the event
        unsigned int    m_nHits;                   ///< The number of hits in the event
        unsigned int    m_nClusters;               ///< The number of clusters in the event
        unsigned int    m_nCandidates;             ///< The total number of vertex candidates
    };

    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Default constructor
     */
    SvmVertexSelectionAlgorithm();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    typedef std::pair<pandora::CartesianVector, pandora::CartesianVector> ClusterEndPoints;
    typedef std::map<const pandora::Cluster *const, ClusterEndPoints> ClusterEndPointsMap;
    typedef std::vector<DoubleVector> FeatureListVector;
    typedef std::vector<pandora::VertexVector> VectorOfVertexVectors;

    /**
     *  @brief  Get the vertex score list
     *
     *  @param  vertexVector the vector of vertices
     *  @param  beamConstants the beam constants
     *  @param  kdTreeU the hit kd tree for the U view
     *  @param  kdTreeV the hit kd tree for the V view
     *  @param  kdTreeW the hit kd tree for the W view
     *  @param  vertexScoreList the vertex score list to fill
     */
    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;

    /**
     *  @brief  Calculate the shower cluster map for a cluster list
     *
     *  @param  inputClusterList the input cluster list
     *  @param  showerClusterList the shower cluster list to populate
     */
    void CalculateShowerClusterList(const pandora::ClusterList &inputClusterList, ShowerClusterList &showerClusterList) const;

    /**
     *  @brief  Add the endpoints of any shower-like clusters to the map
     *
     *  @param  clusterList the list of clusters
     *  @param  showerLikeClusters the list of shower-like clusters to populate
     *  @param  clusterEndPointsMap the map of shower-like cluster endpoints to populate
     */
    void GetShowerLikeClusterEndPoints(const pandora::ClusterList &clusterList, pandora::ClusterList &showerLikeClusters,
        ClusterEndPointsMap &clusterEndPointsMap) const;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;

    /**
     * @brief   Populate kd tree with information about hits in a provided list of clusters
     *
     * @param   clusterList the list of clusters
     * @param   kdTree to receive the populated kd tree
     * @param   hitToClusterMap to receive the populated hit to cluster map
     */
    void PopulateKdTree(const pandora::ClusterList &clusterList, HitKDTree2D &kdTree, HitToClusterMap &hitToClusterMap) const;

    /**
     *  @brief  Try to add an available cluster to a given shower cluster, using shower clustering approximation
     *
     *  @param  clusterEndPointsMap the map of shower-like cluster endpoints
     *  @param  availableShowerLikeClusters the list of shower-like clusters still available
     *  @param  pCluster the cluster in the shower cluster from which to consider distances
     *  @param  showerCluster the shower cluster
     *
     *  @return boolean
     */
    bool AddClusterToShower(const ClusterEndPointsMap &clusterEndPointsMap, pandora::ClusterList &availableShowerLikeClusters,
        const pandora::Cluster *const pCluster, pandora::ClusterList &showerCluster) const;

    /**
     *  @brief  Try to add an available cluster to a given shower cluster, using cluster hit positions cached in kd tree
     *
     *  @param  kdTree the kd tree, used purely for efficiency in events with large hit multiplicity
     *  @param  hitToClusterMap the hit to cluster map, used to interpret kd tree findings
     *  @param  availableShowerLikeClusters the list of shower-like clusters still available
     *  @param  pCluster the cluster in the shower cluster from which to consider distances
     *  @param  showerCluster the shower cluster
     *
     *  @return boolean
     */
    bool AddClusterToShower(HitKDTree2D &kdTree, const HitToClusterMap &hitToClusterMap, pandora::ClusterList &availableShowerLikeClusters,
        const pandora::Cluster *const pCluster, pandora::ClusterList &showerCluster) const;

    /**
     *  @brief  Calculate the event parameters
     *
     *  @param  clusterListU the U-view cluster list
     *  @param  clusterListV the V-view cluster list
     *  @param  clusterListW the W-view cluster list
     *  @param  vertexVector the vector of vertex candidates
     *
     *  @return the EventFeatureInfo object
     */
    EventFeatureInfo CalculateEventFeatures(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV,
        const pandora::ClusterList &clusterListW, const pandora::VertexVector &vertexVector) const;

    /**
     *  @brief  Increment the showery hit parameters for a cluster list
     *
     *  @param  clusterList the cluster list
     *  @param  nShoweryHits the number of showery hits
     *  @param  nHits the number of hits
     *  @param  eventEnergy the event energy
     */
    void IncrementShoweryParameters(const pandora::ClusterList &clusterList, unsigned int &nShoweryHits, unsigned int &nHits,
       float &eventEnergy) const;

    /**
     *  @brief  Find whether a cluster is shower-like
     *
     *  @param  pCluster the cluster
     *
     *  @return boolean
     */
    bool IsClusterShowerLike(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get the event shape features
     *
     *  @param  clusterList the cluster list
     *  @param  eventVolume the event volume
     *  @param  longitudinality the event longitudinality
     */
    void GetEventShapeFeatures(const pandora::ClusterList &clusterList, float &eventVolume, float &longitudinality) const;

    /**
     *  @brief  Update the min/max coordinate spans
     *
     *  @param  minPositionCoord the min position coordinate
     *  @param  maxPositionCoord the max position coordinate
     *  @param  minCoord the current min coordinate
     *  @param  maxCoord the current max coordinate
     */
    void UpdateSpanCoordinate(const float minPositionCoord, const float maxPositionCoord, pandora::InputFloat &minCoord,
        pandora::InputFloat &maxCoord) const;

    /**
     *  @brief  Get the coordinate span
     *
     *  @param  minCoord the min coordinate
     *  @param  maxCoord the max coordinate
     *
     *  @return the coordinate span
     */
    float GetCoordinateSpan(const pandora::InputFloat &minCoord, const pandora::InputFloat &maxCoord) const;

    /**
     *  @brief  Add the event features to a vector in the correct order
     *
     *  @param  eventFeatureInfo the event feature info
     *  @param  featureVector the vector of doubles to append
     */
    void AddEventFeaturesToVector(const EventFeatureInfo &eventFeatureInfo, DoubleVector &featureVector) const;

    /**
     *  @brief  Populate the vertex feature info map for a given vertex
     *
     *  @param  beamConstants the beam constants
     *  @param  clusterListMap the cluster list map
     *  @param  slidingFitDataListMap the sliding fit data list map
     *  @param  showerClusterListMap the shower cluster list map
     *  @param  kdTreeMap the kd tree map
     *  @param  pVertex the vertex
     *  @param  vertexFeatureInfoMap the map to populate
     */
    void PopulateVertexFeatureInfoMap(const BeamConstants &beamConstants, const ClusterListMap &clusterListMap,
        const SlidingFitDataListMap &slidingFitDataListMap, const ShowerClusterListMap &showerClusterListMap, const KDTreeMap &kdTreeMap,
        const pandora::Vertex *const pVertex, VertexFeatureInfoMap &vertexFeatureInfoMap) const;

    /**
     *  @brief  Populate the initial vertex score list for a given vertex
     *
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  pVertex the vertex
     *  @param  initialScoreList the score list to populate
     */
    void PopulateInitialScoreList(VertexFeatureInfoMap &vertexFeatureInfoMap, const pandora::Vertex *const pVertex,
        VertexScoreList &initialScoreList) const;

    /**
     *  @brief  Get the list of top-N separated vertices
     *
     *  @param  initialScoreList the initial score list
     *  @param  bestRegionVertices the list of best region vertices populate
     */
    void GetBestRegionVertices(VertexScoreList &initialScoreList, pandora::VertexVector &bestRegionVertices) const;

    /**
     *  @brief  Produce the region and vertex training sets
     *
     *  @param  vertexVector the vector of all vertices
     *  @param  bestRegionVertices the best region vertices
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  eventFeatureList the list of event features
     *  @param  kdTreeMap
     */
    void ProduceTrainingSets(const pandora::VertexVector &vertexVector, const pandora::VertexVector &bestRegionVertices,
        VertexFeatureInfoMap &vertexFeatureInfoMap, const DoubleVector &eventFeatureList,const KDTreeMap &kdTreeMap) const;

    /**
     *  @brief  Calculate the r/phi scores for the vertices in a vector, possibly erasing those that fail the fast score test
     *
     *  @param  vertexVector the vector of vertices
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  kdTreeMap the kd tree map
     */
    void CalculateRPhiScores(pandora::VertexVector &vertexVector, VertexFeatureInfoMap &vertexFeatureInfoMap, const KDTreeMap &kdTreeMap) const;

    /**
     *  @brief  Get the interaction type string
     *
     *  @return the interaction type string
     */
    std::string GetInteractionType() const;

    /**
     *  @brief  Produce a set of training examples for a binary classifier
     *
     *  @param  vertexVector the vector of vertices to use
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  coinFlip a distribution for producing coin flips
     *  @param  generator the random number generator
     *  @param  interactionType the interaction type string
     *  @param  trainingOutputFile the training set output file
     *  @param  eventFeatureList the event feature list
     *  @param  maxRadius the maximum allowed radius for the 'best' vertex
     *  @param  useRPhi whether to include the r/phi feature
     *
     *  @return address of the vertex used as the 'best' vertex in the classifier
     */
    const pandora::Vertex * ProduceTrainingExamples(const pandora::VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
        std::bernoulli_distribution &coinFlip, std::mt19937 &generator, const std::string &interactionType, const std::string &trainingOutputFile,
        const DoubleVector &eventFeatureList, const float maxRadius, const bool useRPhi) const;

    /**
     *  @brief  Use the MC information to get the best vertex from a list
     *
     *  @param  vertexVector the vector of vertices
     *  @param  pBestVertex address of the best vertex
     *  @param  bestVertexDr dR of the best vertex
     */
    void GetBestVertex(const pandora::VertexVector &vertexVector, const pandora::Vertex *&pBestVertex, float &bestVertexDr) const;

    /**
     *  @brief  Add the vertex features to a vector in the correct order
     *
     *  @param  vertexFeatureInfo the vertex feature info
     *  @param  featureVector the vector of floats to append
     *  @param  useRPhi whether to include the r/phi feature
     */
    void AddVertexFeaturesToVector(const VertexFeatureInfo &vertexFeatureInfo, DoubleVector &featureVector, const bool useRPhi) const;

    /**
     *  @brief  Used a binary classifier to compare a set of vertices and pick the best one
     *
     *  @param  vertexVector the vector of vertices
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  eventFeatureList the event feature list
     *  @param  supportVectorMachine the support vector machine classifier
     *  @param  useRPhi whether to include the r/phi feature
     *
     *  @return address of the best vertex
     */
    const pandora::Vertex * CompareVertices(const pandora::VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
        const DoubleVector &eventFeatureList, const SupportVectorMachine &supportVectorMachine, const bool useRPhi) const;

    /**
     *  @brief  Populate the final vertex score list using the r/phi score to find the best vertex in the vicinity
     *
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  pFavouriteVertex address of the favourite vertex
     *  @param  vertexVector the vector of all vertex candidates
     *  @param  finalVertexScoreList the final vertex score list to populate
     */
    void PopulateFinalVertexScoreList(const VertexFeatureInfoMap &vertexFeatureInfoMap, const pandora::Vertex *const pFavouriteVertex,
        const pandora::VertexVector &vertexVector, VertexScoreList &finalVertexScoreList) const;

    VertexFeatureTool::FeatureToolVector    m_featureToolVector;            ///< The feature tool vector
    std::string                             m_filePathEnvironmentVariable;  ///< The environment variable providing a list of paths to svm files
    std::string                             m_svmFileName;                  ///< The Svm file name
    std::string                             m_regionSvmName;                ///< The name of the region Svm to find
    std::string                             m_vertexSvmName;                ///< The name of the vertex Svm to find
    SupportVectorMachine                    m_svMachineRegion;              ///< The region support vector machine
    SupportVectorMachine                    m_svMachineVertex;              ///< The vertex support vector machine

    bool                  m_trainingSetMode;                      ///< Whether to train
    bool                  m_allowClassifyDuringTraining;          ///< Whether classification is allowed during training
    float                 m_mcVertexXCorrection;                  ///< The correction to the x-coordinate of the MC vertex position
    std::string           m_trainingOutputFileRegion;             ///< The training output file for the region Svm
    std::string           m_trainingOutputFileVertex;             ///< The training output file for the vertex Svm
    std::string           m_mcParticleListName;                   ///< The MC particle list for creating training examples
    std::string           m_caloHitListName;                      ///< The 2D CaloHit list name

    pandora::StringVector m_inputClusterListNames;                ///< The list of cluster list names
    unsigned int          m_minClusterCaloHits;                   ///< The min number of hits parameter in the energy score
    unsigned int          m_slidingFitWindow;                     ///< The layer window for the sliding linear fits
    float                 m_minShowerSpineLength;                 ///< The minimum length at which all are considered to be tracks

    float                 m_beamDeweightingConstant;              ///< The beam deweighting constant for the initial region score list
    float                 m_localAsymmetryConstant;               ///< The local asymmetry constant for the initial region score list
    float                 m_globalAsymmetryConstant;              ///< The global asymmetry constant for the initial region score list
    float                 m_showerAsymmetryConstant;              ///< The shower asymmetry constant for the initial region score list
    float                 m_energyKickConstant;                   ///< The energy kick constant for the initial region score list

    float                 m_showerClusteringDistance;             ///< The shower clustering distance
    unsigned int          m_minShowerClusterHits;                 ///< The minimum number of shower cluster hits
    bool                  m_useShowerClusteringApproximation;     ///< Whether to use the shower clustering distance approximation
    float                 m_regionRadius;                         ///< The radius for a vertex region
    float                 m_rPhiFineTuningRadius;                 ///< The maximum distance the r/phi tune can move a vertex
    float                 m_maxTrueVertexRadius;                  ///< The maximum distance at which a vertex candidate can be considered the 'true' vertex
    bool                  m_useRPhiFeatureForRegion;              ///< Whether to use the r/phi feature for the region vertex
    bool                  m_dropFailedRPhiFastScoreCandidates;    ///< Whether to drop candidates that fail the r/phi fast score test
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline SvmVertexSelectionAlgorithm::VertexFeatureInfo::VertexFeatureInfo(const float beamDeweighting, const float rPhiFeature, const float energyKick,
    const float localAsymmetry, const float globalAsymmetry, const float showerAsymmetry) :
    m_beamDeweighting(beamDeweighting),
    m_rPhiFeature(rPhiFeature),
    m_energyKick(energyKick),
    m_localAsymmetry(localAsymmetry),
    m_globalAsymmetry(globalAsymmetry),
    m_showerAsymmetry(showerAsymmetry)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SvmVertexSelectionAlgorithm::EventFeatureInfo::EventFeatureInfo(const float eventShoweryness, const float eventEnergy,
    const float eventVolume, const float longitudinality, const unsigned int nHits, const unsigned int nClusters, const unsigned int nCandidates) :
    m_eventShoweryness(eventShoweryness),
    m_eventEnergy(eventEnergy),
    m_eventVolume(eventVolume),
    m_longitudinality(longitudinality),
    m_nHits(nHits),
    m_nClusters(nClusters),
    m_nCandidates(nCandidates)
{
}

} // namespace lar_content

#endif // #ifndef LAR_SVM_VERTEX_SELECTION_ALGORITHM_H
