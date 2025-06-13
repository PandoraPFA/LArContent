/**
 *  @file   larpandoracontent/LArUtility/RecoTree.h
 *
 *  @brief  Header file for the reco tree class.
 *
 *  $Log: $
 */
#ifndef LAR_RECO_TREE_H
#define LAR_RECO_TREE_H 1

#include "Objects/OrderedCaloHitList.h"

#include "larpandoracontent/LArUtility/KalmanFilter.h"

#include <memory>
#include <vector>

namespace lar_content
{

/**
 *  @brief  The RecoTree class builds sets of candidate provisional clusters
 */
class RecoTree
{
public:
    /**
     *  @brief  Default constructor
     *
     *  @param  orderedCaloHits the ordered calo hit list
     *  @param  ambiguousHits the set of ambiguous hits
     */
    RecoTree(const pandora::OrderedCaloHitList &orderedCaloHits, const pandora::CaloHitSet &ambiguousHits);

    /**
     *  @brief  Populate the RecoTree from the collection of ordered calo hits
     */
    void Populate();

private:
    class Node
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pCaloHit the seed calo hit associated with this node
         *  @param  tree the reco tree to which this node belongs
         */
        Node(const pandora::CaloHit *const pSeedHit, RecoTree &tree);

        /**
         *  @brief  Populate the node based on the seed hit and collection of ordered calo hits
         */
        void Populate();

    private:
        /**
         *  @brief  Get the proximity of a calo hit to the candidate cluster
         *
         *  @param  pCaloHit the calo hit to check
         *  @return the proximity value
         */
        float GetProximity(const pandora::CaloHit *const pCaloHit) const;

        /**
         *  @brief  Get the Mahalanobis distance of a calo hit to the candidate cluster
         *
         *  @param  pCaloHit the calo hit to check
         *  @return the Mahalanobis distance value
         */
        float GetMahalanobisDistance(const pandora::CaloHit *const pCaloHit);

        const pandora::CaloHit *const m_pSeedHit; ///< The seed calo hit associated with this node
        RecoTree &m_tree; ///< The reco tree to which this node belongs
        pandora::CaloHitVector m_candidateCluster; ///< The candidate cluster associated with this node
        KalmanFilter2D m_kalmanFilter; ///< The Kalman filter used to grow the candidate cluster
    };
    friend class Node;

    typedef std::vector<std::unique_ptr<Node>> NodeVector;

    const pandora::OrderedCaloHitList &m_orderedCaloHits; ///< The ordered calo hit list
    const pandora::CaloHitSet &m_ambiguousHits; ///< The set of ambiguous hits
    pandora::CaloHitSet m_usedHits; ///< The set of used hits in the reco tree
    NodeVector m_rootNodes; ///< The vector of root nodes in the reco tree
};

} // namespace lar_content

#endif // #ifndef LAR_RECO_TREE_H
