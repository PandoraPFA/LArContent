/**
 *  @file   larpandoracontent/LArHelpers/LArGraph.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArGraph.h"

#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

LArGraph::LArGraph(const bool fullyConnect, const int nSourceEdges) :
    m_fullyConnect{fullyConnect},
    m_nSourceEdges{nSourceEdges}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArGraph::~LArGraph()
{
    for (const Edge *pEdge : m_edges)
        delete pEdge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArGraph::EdgeVector &LArGraph::GetEdges() const
{
    return this->m_edges;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGraph::MakeGraph(const CaloHitList &caloHitList)
{
    Eigen::MatrixXf hitMatrix(caloHitList.size(), 2);
    this->Vectorize(caloHitList, hitMatrix);
    // Edges can be double-counted, so use map of maps to avoid this
    std::map<const CaloHit *, std::map<const CaloHit *, bool>> edgeMap;
    for (int r = 0; r < hitMatrix.rows(); ++r)
    {
        Eigen::RowVectorXf row(2);
        row << hitMatrix(r,0), hitMatrix(r,1);
        Eigen::MatrixXf norms((hitMatrix.rowwise() - row).array().pow(2).rowwise().sum());
        norms(r, 0) = std::numeric_limits<float>::max();
        Eigen::Index index1, index2;
        norms.col(0).minCoeff(&index1);
        auto iter0{caloHitList.begin()};
        std::advance(iter0, r);
        auto iter1{caloHitList.begin()};
        std::advance(iter1, index1);
        edgeMap[*iter0][*iter1] = true;
        edgeMap[*iter1][*iter0] = true;
        norms(index1, 0) = std::numeric_limits<float>::max();
        auto val2{norms.col(0).minCoeff(&index2)};
        auto val3{(hitMatrix.row(index1) - hitMatrix.row(index2)).array().pow(2).rowwise().sum()};
        if (val2 < val3(0))
        {
            auto iter2{caloHitList.begin()};
            std::advance(iter2, index2);
            edgeMap[*iter0][*iter2] = true;
            edgeMap[*iter2][*iter0] = true;
        }
    }
    std::map<const CaloHit *, std::map<const CaloHit *, bool>> usedEdges;
    for (const auto &[pCaloHit1, map] : edgeMap)
    {
        for (const auto &[pCaloHit2, dummy] : map)
        {
            if (usedEdges.find(pCaloHit1) != usedEdges.end())
            {
                std::map<const CaloHit *, bool> &innerMap{usedEdges[pCaloHit1]};
                if (innerMap.find(pCaloHit2) != innerMap.end())
                    continue;
            }
            usedEdges[pCaloHit1][pCaloHit2] = true;
            usedEdges[pCaloHit2][pCaloHit1] = true;
            this->m_edges.emplace_back(new Edge(pCaloHit1, pCaloHit2));
        }
    }
    if (m_fullyConnect)
    {
        HitEdgeMap hitToEdgesMap;
        for (const Edge *const pEdge : this->m_edges)
        {
            hitToEdgesMap[pEdge->m_v0].emplace_back(pEdge);
            hitToEdgesMap[pEdge->m_v1].emplace_back(pEdge);
        }
        HitConnectionsMap graphs;
        this->IdentifyDisconnectedRegions(hitToEdgesMap, graphs);
        while (graphs.size() > 1)
        {
            this->ConnectRegions(graphs, hitToEdgesMap);
            graphs.clear();
            this->IdentifyDisconnectedRegions(hitToEdgesMap, graphs);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGraph::Walk(const CaloHit *const pRootCaloHit, const HitEdgeMap &hitToEdgesMap, CaloHitList &graph, HitUseMap &connectedHits) const
{
    if (connectedHits.find(pRootCaloHit) != connectedHits.end())
        return;
    graph.emplace_back(pRootCaloHit);
    connectedHits[pRootCaloHit] = true;
    const EdgeVector &assocEdges{hitToEdgesMap.at(pRootCaloHit)};
    for (const Edge *const edge : assocEdges)
    {
        if (pRootCaloHit == edge->m_v0)
            this->Walk(edge->m_v1, hitToEdgesMap, graph, connectedHits);
        else
            this->Walk(edge->m_v0, hitToEdgesMap, graph, connectedHits);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGraph::IdentifyDisconnectedRegions(const HitEdgeMap &hitToEdgesMap, HitConnectionsMap &graphs) const
{
    HitUseMap connectedHits;
    for (const auto &[pCaloHit, assocEdges] : hitToEdgesMap)
    {
        if (connectedHits.find(pCaloHit) != connectedHits.end())
            continue;
        graphs[pCaloHit] = CaloHitList{pCaloHit};
        connectedHits[pCaloHit] = true;
        for (const Edge *const edge : assocEdges)
        {
            if (pCaloHit == edge->m_v0)
                this->Walk(edge->m_v1, hitToEdgesMap, graphs[pCaloHit], connectedHits);
            else
                this->Walk(edge->m_v0, hitToEdgesMap, graphs[pCaloHit], connectedHits);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGraph::ConnectRegions(const HitConnectionsMap &graphs, HitEdgeMap &hitToEdgesMap)
{
    std::map<int, std::vector<int>> connectedGraphMap;
    int i{0}, idx1{0}, idx2{0};
    for (const auto &[pCaloHit1, caloHitList1] : graphs)
    {
        ++i;
        Eigen::MatrixXf subGraphMatrix(caloHitList1.size(), 2);
        this->Vectorize(caloHitList1, subGraphMatrix);
        float closestApproach{std::numeric_limits<float>::max()};
        const CaloHit *pClosestHit1{nullptr};
        const CaloHit *pClosestHit2{nullptr};
        int j{0};
        for (const auto &[pCaloHit2, caloHitList2] : graphs)
        {
            ++j;
            auto begin{connectedGraphMap[i].begin()}, end{connectedGraphMap[i].end()};
            if (pCaloHit1 == pCaloHit2 || std::find(begin, end, j) != end)
                continue;
            for (const CaloHit *const pCaloHit : caloHitList2)
            {
                const CartesianVector &pos{pCaloHit->GetPositionVector()};
                Eigen::RowVectorXf row(2);
                row << pos.GetX(), pos.GetZ();
                Eigen::MatrixXf norms((subGraphMatrix.rowwise() - row).array().pow(2).rowwise().sum());
                Eigen::Index index;
                float val{norms.col(0).minCoeff(&index)};
                if (val < closestApproach)
                {
                    auto iter{caloHitList1.begin()};
                    std::advance(iter, index);
                    pClosestHit1 = *iter;
                    pClosestHit2 = pCaloHit;
                    closestApproach = val;
                    idx1 = i;
                    idx2 = j;
                }
            }
        }
        if (pClosestHit1 && pClosestHit2)
        {
            const Edge *pEdge{new Edge(pClosestHit1, pClosestHit2)};
            this->m_edges.emplace_back(pEdge);
            hitToEdgesMap[pEdge->m_v0].emplace_back(pEdge);
            hitToEdgesMap[pEdge->m_v1].emplace_back(pEdge);
            connectedGraphMap[idx1].emplace_back(idx2);
            connectedGraphMap[idx2].emplace_back(idx1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArGraph::Vectorize(const CaloHitList &caloHitList, Eigen::MatrixXf &hitMatrix) const
{
    int i{0};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        hitMatrix(i, 0) = pos.GetX();
        hitMatrix(i, 1) = pos.GetZ();
        ++i;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArGraph::Edge::Edge(const CaloHit *const pVertex0, const CaloHit *const pVertex1) :
    m_v0{pVertex0},
    m_v1{pVertex1}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGraph::Edge::LengthSquared() const
{
    const CartesianVector &pos0{this->m_v0->GetPositionVector()}, &pos1{this->m_v1->GetPositionVector()};
    const float dx{pos0.GetX() - pos1.GetX()}, dz{pos0.GetZ() - pos1.GetZ()};

    return dx*dx + dz*dz;
};

} // namespace lar_content
