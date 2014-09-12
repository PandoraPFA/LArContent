/**
 *  @file   LArContent/include/LArCheating/CheatingPfoCreationAlgorithm.h
 * 
 *  @brief  Header file for the cheating cluster creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_PFO_CREATION_ALGORITHM_H
#define LAR_CHEATING_PFO_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingPfoCreationAlgorithm class
 */
class CheatingPfoCreationAlgorithm : public pandora::Algorithm
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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<int, pandora::ClusterList> IdToClusterListMap;
    typedef std::map<int, pandora::MCParticle*> IdToMCParticleMap;

    /**
     *  @brief  Get a map relating three dimensional mc particle ids to a list of daughter clusters
     * 
     *  @param  pClusterList address of a cluster list
     *  @param  idOffset the offset between two dimensional cluster mc particle ids and the three dimensional mc particle id
     *  @param  idToClusterListMap to receive the populated id to cluster list map
     */
    void GetIdToClusterListMap(const pandora::ClusterList *const pClusterList, const int idOffset, IdToClusterListMap &idToClusterListMap) const;

    /**
     *  @brief  Get a map relating three dimensional mc particle ids to the three dimensional mc particles
     * 
     *  @param  idToMCParticleMap to receive the populated id to mc particle map
     */
    void GetIdToMCParticleMap(IdToMCParticleMap &idToMCParticleMap) const;

    /**
     *  @brief  Create pfos corresponding to the details in a provided id to cluster list map
     * 
     *  @param  idToClusterListMap the id to cluster list map
     *  @param  idToMCParticleMap the id to mc particle map
     */
    void CreatePfos(const IdToClusterListMap &idToClusterListMap, const IdToMCParticleMap &idToMCParticleMap) const;

    /**
     *  @brief  Get the number of hit types containing more than a specified number of hits
     * 
     *  @param  clusterList the cluster list, consider all hits in clusters in this list
     *  @param  nHitsThreshold the threshold number of hits of a specified hit type
     * 
     *  @return the number of good hit types
     */
    unsigned int GetNHitTypesAboveThreshold(const pandora::ClusterList &clusterList, const unsigned int nHitsThreshold) const;

    typedef std::map<pandora::HitType, unsigned int> HitTypeMap;
    typedef std::set<int> ParticleIdList;

    std::string     m_inputClusterListNameU;    ///< The name of the view U cluster list
    std::string     m_inputClusterListNameV;    ///< The name of the view V cluster list
    std::string     m_inputClusterListNameW;    ///< The name of the view W cluster list
    std::string     m_outputPfoListName;        ///< The output pfo list name

    int             m_idOffsetU;                ///< The mc particle parent address offset for u clusters
    int             m_idOffsetV;                ///< The mc particle parent address offset for v clusters
    int             m_idOffsetW;                ///< The mc particle parent address offset for w clusters
    std::string     m_mcParticle3DListName;     ///< The name of the three d mc particle list name

    bool            m_useOnlyAvailableClusters; ///< Whether to consider unavailable clusters when identifying cheated pfos
    unsigned int    m_minGoodHitTypes;          ///< The min number of good hit types in the clusters collected for a given mc particle
    unsigned int    m_nHitsForGoodHitType;      ///< The min number of hits of a particular hit type in order to declare the hit type is good
    ParticleIdList  m_particleIdList;           ///< The list of particle ids to consider for pfo creation; will consider all ids if empty
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CheatingPfoCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CheatingPfoCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_PFO_CREATION_ALGORITHM_H
