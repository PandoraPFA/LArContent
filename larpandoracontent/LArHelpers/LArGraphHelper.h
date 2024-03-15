/**
 *  @file   larpandoracontent/LArHelpers/LArGraphHelper.h
 *
 *  @brief  Header file for the Delaunay triangulation helper class.
 *
 *  $Log: $
 */
#ifndef LAR_GRAPH_HELPER_H
#define LAR_GRAPH_HELPER_H 1

#include "Objects/Cluster.h"

#include <Eigen/Dense>

namespace lar_content
{

/**
 *  @brief  LArGraphHelper class
 */
class LArGraphHelper
{
public:
    class Edge
    {
    public:
        /**
         *  @brief Constructor
         *
         *  @param  vertex0 An endpoint of the edge
         *  @param  vertex1 An endpoint of the edge
         **/
        Edge(const pandora::CaloHit *const pVertex0, const pandora::CaloHit *const pVertex1);

        /**
         *  @brief Compute the square of the length of this edge
         *
         *  @return The square of thelength of this edge
         **/
        float LengthSquared() const;

        /**
         *  @brief Check if another edge shares the same calo hit end points
         *
         *  @param other The edge to compare to this edge
         *
         *  @return true if the edges share the same endpoints, irrespective of direction
         **/
        inline bool operator==(const Edge &other) const { return (this->m_v0 == other.m_v0 && this->m_v1 == other.m_v1) ||
            (this->m_v0 == other.m_v1 && this->m_v1 == other.m_v0); };

        const pandora::CaloHit *const m_v0; ///< An endpoint of the edge
        const pandora::CaloHit *const m_v1; ///< An endpoint of the edge
    };
    typedef std::vector<const Edge *> EdgeVector;

private:
    typedef std::map<const pandora::CaloHit *, EdgeVector> HitEdgeMap;
    typedef std::map<const pandora::CaloHit *, pandora::CaloHitList> HitConnectionsMap;
    typedef std::map<const pandora::CaloHit *, bool> HitUseMap;

public:
    /**
     *  @brief  Construct a graph from a list of calo hits.
     *
     *  @param  caloHitList the calo hit list containing the hits from which to construct a graph
     *  @param  edges the output vector of edges to populate
     */
    static void MakeGraph(const pandora::CaloHitList &caloHitList, EdgeVector &edges);

    /**
     *  @brief  Walk along node edges to define a connected graph.
     *
     *  @param  pCaloHit the root calo hit of the current subgraph
     *  @param  hitToEdgesMap the map of hits to the edges for which they are a vertex
     *  @param  graph the current graph to be augmented
     *  @param  connectedHits the map of hits in use
     */
    static void Walk(const pandora::CaloHit *const pCaloHit, const HitEdgeMap &hitToEdgesMap, pandora::CaloHitList &graph, HitUseMap &connectedHits);

private:
    /**
     *  @brief  Identify the different disconnected sub graphs in a set of eges.
     *
     *  @param  hitToEdgesMap the input maps from calo hits to edges
     *  @param  graphs the output map between a root calo hit and all connected hits
     */
    static void IdentifyDisconnectedRegions(const HitEdgeMap &hitToEdgesMap, HitConnectionsMap &graphs);

    /**
     *  @brief  Identify the different disconnected sub graphs in a set of eges.
     *
     *  @param  graphs the map between a root calo hit and all connected hits
     *  @param  hitToEdgesMap the (updatable) map from calo hits to edges
     *  @param  edges the (updatable) vector of edges describing the entire (set of) graph(s)
     */
    static void ConnectRegions(const HitConnectionsMap &graphs, HitEdgeMap &hitToEdgesMap, EdgeVector &edges);

    /**
     *  @brief  Convert a list of calo hits into an Eigen matrix.
     *
     *  @param  caloHitList the calo hit list containing the hits from which to construct a graph
     *  @param  hitMatrix the output Eigen matrix
     */
    static void Vectorize(const pandora::CaloHitList &caloHitList, Eigen::MatrixXf &hitMatrix);
};

} // namespace lar_content

#endif // #ifndef LAR_GRAPH_HELPER_H
