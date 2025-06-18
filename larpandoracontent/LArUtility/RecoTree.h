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
#include "Pandora/AlgorithmHeaders.h"

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
    class Node;
    typedef std::vector<std::unique_ptr<Node>> NodeVector;

    /**
     *  @brief  Default constructor
     *
     *  @param  orderedCaloHits the ordered calo hit list
     *  @param  ambiguousHits the set of ambiguous hits
     *  @param  pitch the pitch of the relevant channel (used for proximity checks)
     */
    RecoTree(const pandora::OrderedCaloHitList &orderedCaloHits, const pandora::CaloHitSet &ambiguousHits, const float pitch, const pandora::Pandora &pandora);

    /**
     *  @brief  Populate the RecoTree from the collection of ordered calo hits
     */
    void Populate();

    /**
     *  @brief  Add the ambiguous hits to the appropriate node
     */
    void ClusterAmbiguousHits();

    /**
     *  @brief  Get the root nodes of the reco tree
     *
     *  @return the vector of root nodes
     */
    const NodeVector &GetRootNodes() const;

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

        /**
         *  @brief  Get the hits associated with this node
         *
         *  @return the hits associated with this node
         */
        const pandora::CaloHitVector &GetHits() const;

        /**
         *  @brief  Add a single (ambiguous) hit to the candidate cluster associated with this node
         *
         *  @param  pCaloHit the calo hit to add
         *  @param  addAtEnd whether to add the hit at the end of the candidate cluster (true) or at the front (false)
         */
        void AddHit(const pandora::CaloHit *const pCaloHit, const bool addAtEnd);

        /**
         *  @brief  Get the closest approach of a calo hit to the candidate cluster
         *
         *  @param  pCaloHit the calo hit to check
         *  @param  pClosestHit the closest hit in the candidate cluster (output parameter)
         *
         *  @return the closest approach value
         */
        float GetClosestApproach(const pandora::CaloHit *const pCaloHit, const pandora::CaloHit *&pClosestHit) const;

        /**
         *  @brief  Get the proximity of a calo hit to the candidate cluster. This is an ordered check, so it's the proximity
         *          to the back of the cluster.
         *
         *  @param  pCaloHit the calo hit to check
         *  @param  centralProximity the proximity when considering hit centres
         *  @param  boundaryProximity the proximity when considering hit boundaries
         *
         *  @return the proximity value (the function attempts to pick the "most appropriate" of central or boundary depending
         *          on the configuration of the hits)
         */
        float GetProximity(const pandora::CaloHit *const pCaloHit, float &centralProximity, float &boundaryProximity) const;

        /**
         *  @brief  Get the Mahalanobis distance of a calo hit to the candidate cluster
         *
         *  @param  pCaloHit the calo hit to check
         *  @return the Mahalanobis distance value
         */
        float GetMahalanobisDistance(const pandora::CaloHit *const pCaloHit);

    private:
        const pandora::CaloHit *const m_pSeedHit; ///< The seed calo hit associated with this node
        RecoTree &m_tree; ///< The reco tree to which this node belongs
        pandora::CaloHitVector m_candidateCluster; ///< The candidate cluster associated with this node
        KalmanFilter2D m_kalmanFilter; ///< The Kalman filter used to grow the candidate cluster
    };
    friend class Node;

private:
    template <class T>
    double WalkThroughCluster(T iter, const T endIter, const Eigen::VectorXd &t, KalmanFilter2D &kalmanFilter);

    const pandora::OrderedCaloHitList &m_orderedCaloHits; ///< The ordered calo hit list
    const pandora::CaloHitSet &m_ambiguousHits; ///< The set of ambiguous hits
    const float m_pitch; ///< The pitch of the relevant channel (used for proximity checks)
    pandora::CaloHitSet m_usedHits; ///< The set of used hits in the reco tree
    NodeVector m_rootNodes; ///< The vector of root nodes in the reco tree
    const pandora::Pandora &m_pandora; ///< The Pandora instance associated with this reco tree 
};

} // namespace lar_content

#endif // #ifndef LAR_RECO_TREE_H
