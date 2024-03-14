/**
 *  @file   larpandoracontent/LArHelpers/LArGraphHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArGraphHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

void LArGraphHelper::MakeGraph(const CaloHitList &caloHitList, EdgeVector &edges)
{
    Eigen::MatrixXf hits(caloHitList.size(), 2);
    int i{0};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        hits(i, 0) = pos.GetX();
        hits(i, 1) = pos.GetZ();
        ++i;
    }
    // Edges can be double-counted, so use map of maps to avoid this
    std::map<const CaloHit *, std::map<const CaloHit *, bool>> edgeMap;
    for (int r = 0; r < hits.rows(); ++r)
    {
        Eigen::RowVectorXf row(2);
        row << hits(r,0), hits(r,1);
        Eigen::MatrixXf norms((hits.rowwise() - row).array().pow(2).rowwise().sum());
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
        auto val3{(hits.row(index1) - hits.row(index2)).array().pow(2).rowwise().sum()};
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
            edges.emplace_back(new Edge(pCaloHit1, pCaloHit2));
        }
    }
    std::map<const CaloHit *, EdgeVector> hitToEdgesMap;
    for (const Edge *const pEdge : edges)
    {
        hitsToEdgeMap[pEdge->m_v0].emplace_back(pEdge);
        hitsToEdgeMap[pEdge->m_v1].emplace_back(pEdge);
    }
    // Walk through the graph, starting from each unvisited hit, to determine disconnected regions
    // stored these hits in a map indexed on the first hit
    // Then find the closest approach between each pair of disconnected regions, noting which hits
    // should be connected
    // Identify the closest approach overall and connect the two regions
    // Update the close approaches to reflect the change in groups, i.e. if we just connected A and B
    // to form A', then previous close approaches to A and B are now close to A', but still with
    // the same hits
    // Find the next closest approach, connect, repeat until no disconnected regions
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArGraphHelper::Edge::Edge(const CaloHit *const pVertex0, const CaloHit *const pVertex1) :
    m_v0{pVertex0},
    m_v1{pVertex1}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArGraphHelper::Edge::LengthSquared() const
{
    const CartesianVector &pos0{this->m_v0->GetPositionVector()}, &pos1{this->m_v1->GetPositionVector()};
    const float dx{pos0.GetX() - pos1.GetX()}, dz{pos0.GetZ() - pos1.GetZ()};

    return dx*dx + dz*dz;
};

} // namespace lar_content
