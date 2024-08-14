/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/MvaLowEClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the mva lowe cluster merging base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MVA_LOW_E_CLUSTER_MERGING_ALGORITHM_H
#define LAR_MVA_LOW_E_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"
#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"

#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
 *  @brief  MvaLowEClusterMergingAlgorithm class
 */
template <typename T>
class MvaLowEClusterMergingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MvaLowEClusterMergingAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~MvaLowEClusterMergingAlgorithm();

protected:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
      *  @brief  Get the MC particle for a given cluster, caching to a map.
      *
      *  @param  cluster         The current cluster to lookup
      *  @param  clusterToMCMap  Map from Cluster to MC to cache results.
      */
    const pandora::MCParticle* GetMCForCluster(const pandora::Cluster *const cluster, std::map<const pandora::Cluster*,
        const pandora::MCParticle*> &clusterToMCMap) const;

    bool IsValidToUse(const pandora::Cluster *const cluster, std::map<const pandora::Cluster*, bool> &clusterIsUsed) const;

    double Angle(const pandora::CartesianVector vector1, const pandora::CartesianVector vector2) const;

    const pandora::CaloHitList EdgeHitFinder(const pandora::Cluster *const cluster, pandora::CaloHitList &clusterEdgeHits) const;

    bool ClusterTool(std::vector<std::string> featureOrder, LArMvaHelper::MvaFeatureMap featureMap) const;

    pandora::StatusCode EdgeHitComparer(const pandora::ClusterList *const pClusterList, const std::string &listName) const;

    pandora::StringVector  m_inputClusterListNames; ///< The names of the input cluster lists.
    std::string            m_mcParticleListName;    ///< Input MC particle list name.
    bool                   m_trainingSetMode;       ///< Whether to train
    bool                   m_enableProbability;     ///< Whether to use probabilities instead of binary classification
    float                  m_minProbabilityCut;     ///< The minimum probability to label a cluster as track-like
    std::string            m_treeName;              ///< Input tree name for ROOT.
    std::string            m_fileName;              ///< Input file name for ROOT.
    int                    m_event;                 ///< Event Number Counter
    float                  m_maxClusterFraction;    ///< The maximum fraction a cluster can be contaminated by to be considered clean.
    float                  m_minNCaloHits;          ///< The minimum number of hits for a cluster to be deemed true for IsAvailableToUse.
    bool                   m_writeTree;             ///< Whether a tree should be output with recorded parameters.
    float                  m_upperHitThreshold;     ///< Max number of hits for cluster to be considered
    std::string            m_trainingOutputFile;    ///< The training output file
    std::string            m_filePathEnvironmentVariable; ///< The environment variable providing a list of paths to mva files
    std::string            m_mvaFileName;           ///< The mva input file
    std::string            m_mvaName;               ///< The name of the mva to find
    T                      m_mva;                   ///< The mva
    float                  m_countHitsThreshold;    ///< A cut on whether cluster merges will occur depending on total event hits
    std::string            m_vertexListName;        ///< Input Vertex List name for vertex based calculation
    float                  m_contactThreshold;      ///< Distance value for hits to be considered in contact
    float                  m_proximityThreshold;    ///< Distance value for hits to be considered in proximity   
    float                  m_divisions;             ///< Number of sectors to search in with Edge Hit Finder Function
    float                  m_sectorTolerance;       ///< Tolerance in radians for dot product between sector and centroid to CaloHit vector



};

typedef MvaLowEClusterMergingAlgorithm<AdaBoostDecisionTree> BdtLowEClusterMergingAlgorithm;

} // namespace lar_content

#endif // #ifndef LAR_MVA_LOW_E_CLUSTER_MERGING_ALGORITHM_H
