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

RecoTree::RecoTree(const OrderedCaloHitList &orderedCaloHits, const CaloHitSet &ambiguousHits, const float pitch) :
    m_orderedCaloHits(orderedCaloHits),
    m_ambiguousHits(ambiguousHits),
    m_pitch(pitch),
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
                m_rootNodes.emplace_back(std::make_unique<Node>(pCaloHit, *this));
                const std::unique_ptr<Node> &pNode{m_rootNodes.back()};
                pNode->Populate();
            }
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const RecoTree::NodeVector &RecoTree::GetRootNodes() const
{
    return m_rootNodes;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

RecoTree::Node::Node(const CaloHit *const pSeedHit, RecoTree &tree) :
    m_pSeedHit(pSeedHit),
    m_tree(tree),
    m_candidateCluster({pSeedHit}),
    m_kalmanFilter{KalmanFilter2D(1, 0.0625f * m_tree.m_pitch * m_tree.m_pitch, 0.25f * m_tree.m_pitch * m_tree.m_pitch,  Eigen::VectorXd(2), 10000.f)}
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
    const CaloHit *latestHit{m_pSeedHit};
    while (true)
    {
        const CaloHit *pBestHit{nullptr};
        CaloHitList *currentLayerHitsPtr{nullptr}, *nextLayerHitsPtr{nullptr};
        m_tree.m_orderedCaloHits.GetCaloHitsInPseudoLayer(latestHit->GetPseudoLayer(), currentLayerHitsPtr);
        m_tree.m_orderedCaloHits.GetCaloHitsInPseudoLayer(latestHit->GetPseudoLayer() + 1, nextLayerHitsPtr);
        const CaloHitList currentLayerHits(currentLayerHitsPtr ? *currentLayerHitsPtr : CaloHitList());
        const CaloHitList nextLayerHits(nextLayerHitsPtr ? *nextLayerHitsPtr : CaloHitList());

        // This can either be proximity or Mahalanobis distance based, depending on the size of the cluster
        float bestMetric{std::numeric_limits<float>::max()};
        for (const CaloHitList &hitList : {currentLayerHits, nextLayerHits})
        {
            for (const CaloHit *const pCurrentHit : hitList)
            {
                // In the current layer we should skip used and ambiguous hits
                if (m_tree.m_usedHits.count(pCurrentHit) > 0)
                    continue;
                if (m_tree.m_ambiguousHits.count(pCurrentHit) > 0)
                    continue;

                // Otherwise, check proximity/consistency and see if we should add this hit to the candidate cluster
                if (m_candidateCluster.size() < 2)
                {
                    float centralProximity{0.f}, boundaryProximity{0.f};
                    const float proximity{this->GetProximity(pCurrentHit, centralProximity, boundaryProximity)};
                    if (proximity < (1.07f * m_tree.m_pitch) && proximity < bestMetric)
                    {
                        if (proximity < boundaryProximity)
                        {
                            pBestHit = pCurrentHit;
                            bestMetric = proximity;
                        }
                        else if (boundaryProximity < 0.1f)
                        {
                            // If the returned proximity is the boundary proximity, we should be more strict to avoid
                            // gaps between hits being used to seed the cluster
                            pBestHit = pCurrentHit;
                            bestMetric = boundaryProximity;
                        }
                    }
                }
                else
                {
                    float centralProximity{0.f}, boundaryProximity{0.f};
                    const float proximity{this->GetProximity(pCurrentHit, centralProximity, boundaryProximity)};
                    if (proximity < (1.07f * m_tree.m_pitch))
                    {
                        const float mahalanobisDistance{this->GetMahalanobisDistance(pCurrentHit)};
                        if (mahalanobisDistance < (1.1f * m_tree.m_pitch) && mahalanobisDistance < bestMetric)
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
        }
        // Having processed the current layer and the next layer, we should check to see if a best hit was found and update the candidate
        // cluster if so. Otherwise, we stop building this candidate cluster.
        if (pBestHit)
        {
            latestHit = pBestHit;
            m_candidateCluster.emplace_back(pBestHit);
            m_tree.m_usedHits.insert(pBestHit);

            // Update the Kalman filter with this hit
            const CartesianVector &pos{pBestHit->GetPositionVector()};
            Eigen::VectorXd x(2);
            x << pos.GetX(), pos.GetZ();
            m_kalmanFilter.Predict();
            m_kalmanFilter.Update(x);
        }
        else
        {
            break;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CaloHitVector &RecoTree::Node::GetHits() const
{
    return m_candidateCluster;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float RecoTree::Node::GetProximity(const CaloHit *const pCaloHit, float &centralProximity, float &boundaryProximity) const
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
        {
            centralProximity = std::fabs(position1.GetZ() - position2.GetZ());
            boundaryProximity = std::fabs(xhi1 - xlo2);
            return boundaryProximity;
        }
        else if (xhi2 <= xlo1)
        {
            centralProximity = std::fabs(position1.GetZ() - position2.GetZ());
            boundaryProximity = std::fabs(xhi2 - xlo1);
            return boundaryProximity;
        }
        else
        {
            centralProximity = std::fabs(position1.GetZ() - position2.GetZ());
            boundaryProximity = std::numeric_limits<float>::max();
            return centralProximity;
        }
    }
    else
    {
        // Same pseudo layer, the centre of pCaloHit must be to the right of the previously added candidate hit
        const float xlo1{position1.GetX() - width1};
        const float xhi2{position2.GetX() + width2};
        // Note, it is possible to have hits on adjacent channels appear in the same pseudo layer, so we should also check for hit centre
        // proximity
        if (position1.GetZ() > position2.GetZ())
        {
            centralProximity = std::fabs(position1.GetZ() - position2.GetZ());
            boundaryProximity = std::numeric_limits<float>::max();
            return centralProximity;
        }
        else
        {
            centralProximity = std::fabs(position1.GetZ() - position2.GetZ());
            boundaryProximity = std::fabs(xlo1 - xhi2);
            return boundaryProximity;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float RecoTree::Node::GetMahalanobisDistance(const CaloHit *const pCaloHit)
{
    m_kalmanFilter.Predict();
    const CartesianVector &pos{pCaloHit->GetPositionVector()};
    Eigen::VectorXd x(2);
    x << pos.GetX(), pos.GetZ();
    return m_kalmanFilter.GetMahalanobisDistance(x);
}

} // namespace lar_content
