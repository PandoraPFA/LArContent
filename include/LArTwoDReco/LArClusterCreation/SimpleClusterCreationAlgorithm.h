/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SIMPLE_CLUSTER_CREATION_ALGORITHM_H
#define LAR_SIMPLE_CLUSTER_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SimpleClusterCreationAlgorithm class
 */
class SimpleClusterCreationAlgorithm : public pandora::Algorithm
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
    SimpleClusterCreationAlgorithm();

private:
    pandora::StatusCode Run();

    typedef std::map<const pandora::CaloHit*, pandora::CaloHitList> HitAssociationMap;

    /**
     *  @brief Select calo hits for clustering
     *
     *  @param inputList The input list of calo hits
     *  @param outputList The output list of selected calo hits
     */
    void SelectCaloHits(const pandora::CaloHitList* inputList, pandora::CaloHitList &outputList) const;

    /**
     *  @brief Create map of associations between calo hits
     *
     *  @param caloHitList The input list of calo hits
     *  @param hitAssociationMap The map of associations between calo hits
     */
    void BuildAssociationMap(const pandora::CaloHitList &caloHitList, HitAssociationMap &hitAssociationMap) const;

    /**
     *  @brief Create clusters from selected calo hits and their associations
     *
     *  @param caloHitList The input list of calo hits
     *  @param hitAssociationMap The map of associations between calo hits
     */
    void CreateClusters(const pandora::CaloHitList &caloHitList, const HitAssociationMap &hitAssociationMap) const;

    /**
     *  @brief For a given seed calo hits, collect up all the associated calo hits
     *
     *  @param pSeedCaloHit the seed calo hits
     *  @param pCurrentCaloHit a possible associated calo hit
     *  @param hitAssociationMap the map of associations between hits
     *  @param vetoList the list of used calo hits
     *  @param mergeList the list of hits associated with the seed hit
     */
    void CollectAssociatedHits(pandora::CaloHit *pSeedCaloHit, pandora::CaloHit *pCurrentCaloHit,
        const HitAssociationMap &hitAssociationMap, const pandora::CaloHitList &vetoList,
        pandora::CaloHitList &mergeList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float                m_clusteringWindowSquared;      ///< Maximum distance (squared) for two hits to be joined
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SimpleClusterCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new SimpleClusterCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_SIMPLE_CLUSTER_CREATION_ALGORITHM_H
