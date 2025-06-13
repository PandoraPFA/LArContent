/**
 *  @file   larpandoracontent/LArUtility/RecoTree.cc
 *
 *  @brief  Implementation of the RecoTree class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArUtility/RecoTree.h"

using namespace pandora;

namespace lar_content
{

RecoTree::RecoTree(const pandora::OrderedCaloHitList &orderedCaloHits, const CaloHitSet &ambiguousHits) :
    m_orderedCaloHits(orderedCaloHits),
    m_ambiguousHits(ambiguousHits),
    m_usedHits(),
    m_rootNodes()
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void RecoTree::Populate()
{
    // Clear any existing root nodes
    m_rootNodes.clear();

    // Iterate through the ordered calo hits and build the tree
    for (OrderedCaloHitList::const_iterator iter = m_orderedCaloHits.begin(); iter != m_orderedCaloHits.end(); ++iter)
    {
        const CaloHitList &currentLayerHits(*iter->second);
        if (currentLayerHits.empty())
            continue;

        // Create a new for any unused, unambiguous hit
        for (const CaloHit *const pCaloHit : currentLayerHits)
        {
            if (m_ambiguousHits.count(pCaloHit) == 0 && m_usedHits.count(pCaloHit) == 0)
            {
                m_usedHits.insert(pCaloHit);
                // ToDo, figure out the appropriate interface for a Node constructor
                m_rootNodes.emplace_back(std::make_unique<Node>(pCaloHit, *this));
                // Basically we should create a seed based on the closest hit in the current or downstream layer
                // and thereafter collect up the best unambiguous hits, stopping if we hit an ambiguous one that
                // could otherwise match
                // At that point we should probably bail and fall back to the RecoTree class to continue, but there is
                // scope to make this truly tree like and construct branches based on nearest neighbours of the ambiguous hits
                // Ambiguous hits could form nodes in their own right
            }
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

RecoTree::Node::Node(const pandora::CaloHit *const pSeedHit, RecoTree &tree) :
    m_pSeedHit(pSeedHit),
    m_tree(tree),
    m_candidateCluster({pSeedHit}),
    m_kalmanFilter{KalmanFilter2D(1, 0.125 * 0.125, 0.25 * 0.25,  Eigen::VectorXd(2), 1000.f)}
{
    const CartesianVector &pos{pSeedHit->GetPositionVector()};
    Eigen::VectorXd seed(2);
    seed << pos.GetX(), pos.GetZ();
    m_kalmanFilter.Predict();
    m_kalmanFilter.Update(seed);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void RecoTree::Node::Populate()
{
    CaloHitList *currentLayerHitsPtr{nullptr}, *nextLayerHitsPtr{nullptr};
    m_tree.m_orderedCaloHits.GetCaloHitsInPseudoLayer(m_pSeedHit->GetPseudoLayer(), currentLayerHitsPtr);
    m_tree.m_orderedCaloHits.GetCaloHitsInPseudoLayer(m_pSeedHit->GetPseudoLayer() + 1, nextLayerHitsPtr);
    const CaloHitList currentLayerHits(currentLayerHitsPtr ? *currentLayerHitsPtr : CaloHitList());
    const CaloHitList nextLayerHits(nextLayerHitsPtr ? *nextLayerHitsPtr : CaloHitList());

    const CaloHit *pBestHit{nullptr};
    // This can either be proximity or Mahalanobis distance based, depending on the size of the cluster
    float bestMetric{std::numeric_limits<float>::max()};
    for (const CaloHit *const pCurrentHit : currentLayerHits)
    {
        if (pCurrentHit == m_pSeedHit)
            continue;

        // In the current layer we should skip used hits and we should stop if we hit an ambiguous hit
        if (m_tree.m_usedHits.count(pCurrentHit) > 0)
            continue;
        if (m_tree.m_ambiguousHits.count(pCurrentHit) > 0)
            break;

        // Otherwise, check proximity/consistency and see if we should add this hit to the candidate cluster
        if (m_candidateCluster.size() < 2)
        {
            const float proximity{this->GetProximity(pCurrentHit)};
            if (proximity < 0.5f && proximity < bestMetric)
            {
                pBestHit = pCurrentHit;
                bestMetric = proximity;
            }
        }
        else
        {
            const float proximity{this->GetProximity(pCurrentHit)};
            if (proximity < 0.5f)
            {
                const float mahalanobisDistance{this->GetMahalanobisDistance(pCurrentHit)};
                if (mahalanobisDistance < bestMetric)
                {
                    // If we have a candidate hit, check if this one is better
                    if (!pBestHit || mahalanobisDistance < bestMetric)
                    {
                        pBestHit = pCurrentHit;
                        bestMetric = mahalanobisDistance;
                    }
                }
            }
        }
    }

    // Here, having processed the current layer and the next layer, we should check to see if a best hit was found and update the candidate
    // cluster if so. Otherwise, we stop building this candidate cluster.
    if (pBestHit)
    {
        m_candidateCluster.emplace_back(pBestHit);
        m_tree.m_usedHits.insert(pBestHit);

        // Update the Kalman filter with this hit
        const CartesianVector &pos{pBestHit->GetPositionVector()};
        Eigen::VectorXd x(2);
        x << pos.GetX(), pos.GetZ();
        m_kalmanFilter.Predict();
        m_kalmanFilter.Update(x);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float RecoTree::Node::GetProximity(const pandora::CaloHit *const pCaloHit) const
{
    const CartesianVector &position1{pCaloHit->GetPositionVector()};
    const CartesianVector &position2{m_candidateCluster.back()->GetPositionVector()};

    const float width1{0.5f * pCaloHit->GetCellSize1()};
    const float width2{0.5f * m_candidateCluster.back()->GetCellSize1()};
    if (pCaloHit->GetPseudoLayer() > m_candidateCluster.back()->GetPseudoLayer())
    {
        // Different pseudo layers, so the hits could have any ordering in x
        const float xlo1{position1.GetX() - width1}, xhi1{position1.GetX() + width1};
        const float xlo2{position2.GetX() - width2}, xhi2{position2.GetX() + width2};
        // If the hits have no overlap, check the distance between the closests edges, otherwise, use the centres
        if (xhi1 <= xlo2)
            return std::fabs(xhi1 - xlo2);
        else if (xhi2 <= xlo1)
            return std::fabs(xhi2 - xlo1);
        else
            return std::fabs(position1.GetZ() - position2.GetZ());
    }
    else
    {
        // Same pseudo layer, the centre of pCaloHit must be to the right of the previously added candidate hit
        const float xlo1{position1.GetX() - width1};
        const float xhi2{position2.GetX() + width2};
        // Note, it is possible to have hits on adjacent channels appear in the same pseudo layer, so we should also check for hit centre
        // proximity
        if (position1.GetZ() > position2.GetZ())
            return std::fabs(position1.GetZ() - position2.GetZ());
        else
            return std::fabs(xlo1 - xhi2);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float RecoTree::Node::GetMahalanobisDistance(const pandora::CaloHit *const pCaloHit)
{
    m_kalmanFilter.Predict();
    const CartesianVector &pos{pCaloHit->GetPositionVector()};
    Eigen::VectorXd x(2);
    x << pos.GetX(), pos.GetZ();
    return m_kalmanFilter.GetMahalanobisDistance(x);
}

} // namespace lar_content
