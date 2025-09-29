/**
 *  @file   larpandoracontent/LArUtility/RecoTree.cc
 *
 *  @brief  Implementation of the RecoTree class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArEigenHelper.h"

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

    this->ClusterAmbiguousHits();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void RecoTree::ClusterAmbiguousHits()
{
    // Iterate through the ambiguous hits and determine best matching node
    // First, identify the closest approach of each ambiguous hit to each root node
    std::vector<std::vector<float>> closestApproach(m_ambiguousHits.size(), std::vector<float>(m_rootNodes.size(), std::numeric_limits<float>::max()));
    size_t n{0};
    for (const auto &pNode : m_rootNodes)
    {
        const CaloHitVector &nodeHits{pNode->GetHits()};
        Eigen::MatrixXf hitMatrix(nodeHits.size(), 2);
        LArEigenHelper::Vectorize(nodeHits, hitMatrix);
        size_t h{0};
        for (const CaloHit *const pCaloHit : m_ambiguousHits)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            if (m_usedHits.count(pCaloHit) == 0)
            {
                Eigen::RowVectorXf row(2);
                row << pos.GetX(), pos.GetZ();
                Eigen::MatrixXf norms((hitMatrix.rowwise() - row).array().pow(2).rowwise().sum());
                Eigen::Index index;
                norms.col(0).minCoeff(&index);
                closestApproach[h][n] = norms(index, 0);
            }
            ++h;
        }
        ++n;
    }

    // For each ambiguous hit, collect all of the nodes within a given distance for more detailed consideration
    bool madeAllocation{true};
    while (madeAllocation)
    {
        madeAllocation = false;
        size_t h{0};
        for (const CaloHit *const pTargetHit : m_ambiguousHits)
        {
            if (m_usedHits.count(pTargetHit) > 0)
            {
                ++h;
                continue;
            }
            const CartesianVector &targetPos{pTargetHit->GetPositionVector()};
            Eigen::VectorXd t(2);
            t << targetPos.GetX(), targetPos.GetZ();
            double bestMahalanobisDistance{std::numeric_limits<double>::max()};
            float bestProximity{std::numeric_limits<float>::max()};
            Node *pBestNodeMahalanobis{nullptr}, *pBestNodeProximity{nullptr};
            (void)bestProximity;      // ToDo: Use this in the future for proximity based clustering
            (void)pBestNodeProximity; // ToDo: Use this in the future for proximity based clustering
            bool addAtEndMahalanobis{false}, addAtEndProximity{false};
            (void)addAtEndProximity; // ToDo: Use this in the future for proximity based clustering

            size_t c{0};
            for (const auto &pNode : m_rootNodes)
            {
                if (closestApproach[h][c] < m_closeApproachThreshold)
                {
                    const CaloHitVector &nodeHits{m_rootNodes[c]->GetHits()};
                    const float dFront{std::abs((pTargetHit->GetPositionVector() - nodeHits.front()->GetPositionVector()).GetMagnitudeSquared())};
                    const float dBack{std::abs((pTargetHit->GetPositionVector() - nodeHits.back()->GetPositionVector()).GetMagnitudeSquared())};
                    if (dFront > dBack)
                    {
                        // Walk from the front to the back of the cluster and see where the hit would best fit
                        KalmanFilter2D kalmanFilter(1, m_processVarianceCoeff * m_pitch * m_pitch,
                            m_measurementVarianceCoeff * m_pitch * m_pitch, Eigen::VectorXd(2), 10000.f);
                        if (nodeHits.size() > 1)
                        {
                            kalmanFilter.Predict();
                            auto iter{nodeHits.begin()};
                            const CartesianVector &pos{(*iter)->GetPositionVector()};
                            Eigen::VectorXd x(2);
                            x << pos.GetX(), pos.GetZ();
                            kalmanFilter.Update(x);
                            kalmanFilter.Predict();
                            ++iter;

                            const double mahalanobisDistance{this->WalkThroughCluster(iter, nodeHits.end(), t, kalmanFilter)};
                            if (mahalanobisDistance < bestMahalanobisDistance)
                            {
                                bestMahalanobisDistance = mahalanobisDistance;
                                pBestNodeMahalanobis = pNode.get();
                                addAtEndMahalanobis = true;
                            }
                        }
                        else
                        {
                            const CaloHit *pClosestHit{nullptr};
                            const float proximity{pNode->GetClosestApproach(pTargetHit, pClosestHit)};
                            if (proximity < bestProximity)
                            {
                                size_t minLayer{std::min(pClosestHit->GetPseudoLayer(), pTargetHit->GetPseudoLayer())};
                                size_t maxLayer{std::max(pClosestHit->GetPseudoLayer(), pTargetHit->GetPseudoLayer())};
                                CaloHitVector regionHits;
                                for (size_t i = minLayer; i <= maxLayer; ++i)
                                {
                                    CaloHitList *pHits{nullptr};
                                    m_orderedCaloHits.GetCaloHitsInPseudoLayer(i, pHits);
                                    if (pHits)
                                        regionHits.insert(regionHits.end(), pHits->begin(), pHits->end());
                                }
                                if (!LArClusterHelper::HasBlockedPath(regionHits, pClosestHit, pTargetHit))
                                {
                                    bestProximity = proximity;
                                    pBestNodeProximity = pNode.get();
                                    addAtEndProximity = true;
                                }
                            }
                        }
                    }
                    else
                    {
                        // Walk from the back to the front of the cluster and see where the hit would best fit
                        KalmanFilter2D kalmanFilter(1, m_processVarianceCoeff * m_pitch * m_pitch,
                            m_measurementVarianceCoeff * m_pitch * m_pitch, Eigen::VectorXd(2), 10000.f);
                        if (nodeHits.size() > 1)
                        {
                            kalmanFilter.Predict();
                            auto iter{nodeHits.rbegin()};
                            const CartesianVector &pos{(*iter)->GetPositionVector()};
                            Eigen::VectorXd x(2);
                            x << pos.GetX(), pos.GetZ();
                            kalmanFilter.Update(x);
                            kalmanFilter.Predict();
                            ++iter;

                            const double mahalanobisDistance{this->WalkThroughCluster(iter, nodeHits.rend(), t, kalmanFilter)};
                            if (mahalanobisDistance < bestMahalanobisDistance)
                            {
                                bestMahalanobisDistance = mahalanobisDistance;
                                pBestNodeMahalanobis = pNode.get();
                                addAtEndMahalanobis = false;
                            }
                        }
                        else
                        {
                            const CaloHit *pClosestHit{nullptr};
                            const float proximity{pNode->GetClosestApproach(pTargetHit, pClosestHit)};
                            if (proximity < bestProximity)
                            {
                                size_t minLayer{std::min(pClosestHit->GetPseudoLayer(), pTargetHit->GetPseudoLayer())};
                                size_t maxLayer{std::max(pClosestHit->GetPseudoLayer(), pTargetHit->GetPseudoLayer())};
                                CaloHitVector regionHits;
                                for (size_t i = minLayer; i <= maxLayer; ++i)
                                {
                                    CaloHitList *pHits{nullptr};
                                    m_orderedCaloHits.GetCaloHitsInPseudoLayer(i, pHits);
                                    if (pHits)
                                        regionHits.insert(regionHits.end(), pHits->begin(), pHits->end());
                                }
                                if (!LArClusterHelper::HasBlockedPath(regionHits, pClosestHit, pTargetHit))
                                {
                                    bestProximity = proximity;
                                    pBestNodeProximity = pNode.get();
                                    addAtEndProximity = false;
                                }
                            }
                        }
                    }
                }
                ++c;
            }
            ++h;

            // If we have both a proximity and a Mahalanobis distance based node, we should check which one is more appropriate
            // We should require the proximity based distance to be much better to favour it
            if (pBestNodeProximity && pBestNodeMahalanobis)
            {
                if (bestProximity < (m_proximityCoeff * m_pitch) && bestProximity < (m_mahalanobisRescaling * bestMahalanobisDistance))
                    pBestNodeMahalanobis = nullptr;
                else
                    pBestNodeProximity = nullptr;
            }
            // Might want to relax the distance and also potentially check for prediction inside hit for wide hits
            if (pBestNodeMahalanobis && bestMahalanobisDistance < (m_mahalanobisCoeff * m_pitch))
            {
                pBestNodeMahalanobis->AddHit(pTargetHit, addAtEndMahalanobis);
                m_usedHits.insert(pTargetHit);
                madeAllocation = true;
            }
            else if (pBestNodeProximity && bestProximity < (m_proximityCoeff * m_pitch))
            {
                pBestNodeProximity->AddHit(pTargetHit, addAtEndProximity);
                m_usedHits.insert(pTargetHit);
                madeAllocation = true;
            }
            else
            {
                // If we have no suitable node, we should create a new one
                m_rootNodes.emplace_back(std::make_unique<Node>(pTargetHit, *this));
                m_usedHits.insert(pTargetHit);
                madeAllocation = true;

                // Need to rebuild the closest approach matrix since we have added a new node
                closestApproach.clear();
                closestApproach = std::vector<std::vector<float>>(
                    m_ambiguousHits.size(), std::vector<float>(m_rootNodes.size(), std::numeric_limits<float>::max()));
                size_t h1{0};
                for (const auto &pNode : m_rootNodes)
                {
                    const CaloHitVector &nodeHits{pNode->GetHits()};
                    Eigen::MatrixXf hitMatrix(nodeHits.size(), 2);
                    LArEigenHelper::Vectorize(nodeHits, hitMatrix);
                    size_t c1{0};
                    for (const CaloHit *const pCaloHit : m_ambiguousHits)
                    {
                        const CartesianVector &pos{pCaloHit->GetPositionVector()};
                        if (m_usedHits.count(pCaloHit) == 0)
                        {
                            Eigen::RowVectorXf row(2);
                            row << pos.GetX(), pos.GetZ();
                            Eigen::MatrixXf norms((hitMatrix.rowwise() - row).array().pow(2).rowwise().sum());
                            Eigen::Index index;
                            norms.col(0).minCoeff(&index);
                            closestApproach[c1][h1] = norms(index, 0);
                        }
                        ++c1;
                    }
                    ++h1;
                }
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

template <class T>
double RecoTree::WalkThroughCluster(T iter, const T endIter, const Eigen::VectorXd &t, KalmanFilter2D &kalmanFilter)
{
    double bestMahalanobisDistance{std::numeric_limits<float>::max()};
    for (; iter != endIter; ++iter)
    {
        const CaloHit *const pCaloHit{*iter};
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        Eigen::VectorXd x(2);
        x << pos.GetX(), pos.GetZ();
        kalmanFilter.Update(x);
        kalmanFilter.Predict();
        const double mahalanobisDistance{kalmanFilter.GetMahalanobisDistance(t)};
        if (mahalanobisDistance < bestMahalanobisDistance)
            bestMahalanobisDistance = mahalanobisDistance;
    }

    return bestMahalanobisDistance;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void RecoTree::Configure(const float closeApproachThreshold, const float processVarianceCoeff, const float measurementVarianceCoeff,
    const float proximityCoeff, const float mahalanobisCoeff, const float mahalanobisRescaling, const float boundaryProximity)
{
    m_closeApproachThreshold = closeApproachThreshold;
    m_processVarianceCoeff = processVarianceCoeff;
    m_measurementVarianceCoeff = measurementVarianceCoeff;
    m_proximityCoeff = proximityCoeff;
    m_mahalanobisCoeff = mahalanobisCoeff;
    m_mahalanobisRescaling = mahalanobisRescaling;
    m_boundaryProximity = boundaryProximity;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

RecoTree::Node::Node(const CaloHit *const pSeedHit, RecoTree &tree) :
    m_pSeedHit(pSeedHit),
    m_tree(tree),
    m_candidateCluster({pSeedHit}),
    m_kalmanFilter{KalmanFilter2D(1, m_tree.m_processVarianceCoeff * m_tree.m_pitch * m_tree.m_pitch,
        m_tree.m_measurementVarianceCoeff * m_tree.m_pitch * m_tree.m_pitch, Eigen::VectorXd::Zero(2), 10000.f)}
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
                    if (proximity < (m_tree.m_proximityCoeff * m_tree.m_pitch) && proximity < bestMetric)
                    {
                        if (proximity < boundaryProximity)
                        {
                            pBestHit = pCurrentHit;
                            bestMetric = proximity;
                        }
                        else if (boundaryProximity < m_tree.m_boundaryProximity)
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
                    if (proximity < (m_tree.m_proximityCoeff * m_tree.m_pitch))
                    {
                        const float mahalanobisDistance{this->GetMahalanobisDistance(pCurrentHit)};
                        if (mahalanobisDistance < (m_tree.m_mahalanobisCoeff * m_tree.m_pitch) && mahalanobisDistance < bestMetric)
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

void RecoTree::Node::AddHit(const CaloHit *const pCaloHit, const bool addAtEnd)
{
    // Add the hit to the candidate cluster, either at the end or at the front
    if (addAtEnd)
    {
        m_candidateCluster.emplace_back(pCaloHit);
    }
    else
    {
        m_candidateCluster.insert(m_candidateCluster.begin(), pCaloHit);
    }
    // ATTN: We're choosing not to update the Kalman filter here under the assumption that this is run for ambiguous hit additions
    // and that therefore the internal Kalman state is less relevant
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float RecoTree::Node::GetClosestApproach(const CaloHit *const pCaloHit, const CaloHit *&pClosestHit) const
{
    Eigen::MatrixXf hitMatrix(m_candidateCluster.size(), 2);
    LArEigenHelper::Vectorize(m_candidateCluster, hitMatrix);
    const CartesianVector &pos{pCaloHit->GetPositionVector()};
    Eigen::RowVectorXf row(2);
    row << pos.GetX(), pos.GetZ();
    Eigen::MatrixXf norms((hitMatrix.rowwise() - row).array().pow(2).rowwise().sum());
    Eigen::Index index;
    norms.col(0).minCoeff(&index);
    pClosestHit = m_candidateCluster.at(index);

    return norms(index, 0);
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
