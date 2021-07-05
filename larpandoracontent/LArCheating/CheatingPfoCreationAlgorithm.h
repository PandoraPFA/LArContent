/**
 *  @file   larpandoracontent/LArCheating/CheatingPfoCreationAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_PFO_CREATION_ALGORITHM_H
#define LAR_CHEATING_PFO_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingPfoCreationAlgorithm class
 */
class CheatingPfoCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingPfoCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::MCParticle *, pandora::ClusterList> MCParticleToClusterListMap;

    /**
     *  @brief  Get a map relating mc particles to a list of daughter clusters
     *
     *  @param  pClusterList address of a cluster list
     *  @param  mcPrimaryMap the mapping between mc particles and their parents
     *  @param  mcParticleToClusterListMap to receive the populated mc particle to cluster list map
     */
    void GetMCParticleToClusterListMap(const pandora::ClusterList *const pClusterList,
        const LArMCParticleHelper::MCRelationMap &mcPrimaryMap, MCParticleToClusterListMap &mcParticleToClusterListMap) const;

    /**
     *  @brief  Create pfos corresponding to the details in a provided mc particle to cluster list map
     *
     *  @param  mcParticleToClusterListMap the mc particle to cluster list map
     */
    void CreatePfos(const MCParticleToClusterListMap &mcParticleToClusterListMap) const;

    /**
     *  @brief  Get the number of hit types containing more than a specified number of hits
     *
     *  @param  clusterList the cluster list, consider all hits in clusters in this list
     *  @param  nHitsThreshold the threshold number of hits of a specified hit type
     *
     *  @return the number of good hit types
     */
    unsigned int GetNHitTypesAboveThreshold(const pandora::ClusterList &clusterList, const unsigned int nHitsThreshold) const;

    bool IsTrack(const pandora::MCParticle *const pMCParticle) const;
    bool IsShower(const pandora::MCParticle *const pMCParticle) const;

    typedef std::map<pandora::HitType, unsigned int> HitTypeMap;
    typedef std::set<int> ParticleIdList;

    pandora::StringVector m_inputClusterListNames; ///< The names of the input cluster lists
    std::string m_outputPfoListName;               ///< The output pfo list name
    std::string m_outputVertexListName;            ///< The output vertex list name

    bool m_collapseToPrimaryMCParticles; ///< Whether to collapse mc particle hierarchies to primary particles
    std::string m_mcParticleListName;    ///< The mc particle list name

    bool m_useOnlyAvailableClusters;    ///< Whether to consider unavailable clusters when identifying cheated pfos
    bool m_addVertices;                 ///< Whether to add the start vertex to the cheated pfo
    bool m_replaceCurrentVertexList;    ///< Whether to replace current vertex list
    unsigned int m_minGoodHitTypes;     ///< The min number of good hit types in the clusters collected for a given mc particle
    unsigned int m_nHitsForGoodHitType; ///< The min number of hits of a particular hit type in order to declare the hit type is good
    ParticleIdList m_particleIdList;    ///< The list of particle ids to consider for pfo creation; will consider all ids if empty
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_PFO_CREATION_ALGORITHM_H
