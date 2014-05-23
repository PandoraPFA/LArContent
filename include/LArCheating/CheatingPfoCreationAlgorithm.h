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

namespace lar
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

    std::string     m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string     m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string     m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string     m_outputPfoListName;            ///< The output pfo list name

    int             m_idOffsetU;                    ///< The mc particle parent address offset for u clusters
    int             m_idOffsetV;                    ///< The mc particle parent address offset for v clusters
    int             m_idOffsetW;                    ///< The mc particle parent address offset for w clusters
    std::string     m_mcParticle3DListName;         ///< The name of the three d mc particle list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CheatingPfoCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CheatingPfoCreationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_CHEATING_PFO_CREATION_ALGORITHM_H
