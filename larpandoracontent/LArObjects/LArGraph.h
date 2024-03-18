/**
 *  @file   larpandoracontent/LArHelpers/LArGraph.h
 *
 *  @brief  Header file for the LArGraph class.
 *
 *  $Log: $
 */
#ifndef LAR_GRAPH_H
#define LAR_GRAPH_H 1

#include "Objects/Cluster.h"

#include <Eigen/Dense>

namespace lar_content
{

/**
 *  @brief  LArGraph class
 */
class LArGraph
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
     *  @brief  Constructor
     *
     *  @param  fullyConnect Whether or not to connect disconnected regions
     *  @param  nSourceEdges The number of edges emerging from a source node
     */
    LArGraph(const bool fullyConnect = true, const int nSourceEdges = 2);

    /**
     *  @brief  Destructor
     */
    ~LArGraph();

    /**
     *  @brief  Retrieve the edges for this graph
     *
     *  @return The edges for this graph
     */
    const EdgeVector &GetEdges() const;

    /**
     *  @brief  Construct a graph from a list of calo hits.
     *
     *  @param  caloHitList the calo hit list containing the hits from which to construct a graph
     */
    void MakeGraph(const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Walk along node edges to define a connected graph.
     *
     *  @param  pCaloHit the root calo hit of the current subgraph
     *  @param  hitToEdgesMap the map of hits to the edges for which they are a vertex
     *  @param  graph the current graph to be augmented
     *  @param  connectedHits the map of hits in use
     */
    void Walk(const pandora::CaloHit *const pCaloHit, const HitEdgeMap &hitToEdgesMap, pandora::CaloHitList &graph, HitUseMap &connectedHits) const;

private:
    /**
     *  @brief  Identify the different disconnected sub graphs in a set of eges.
     *
     *  @param  hitToEdgesMap the input maps from calo hits to edges
     *  @param  graphs the output map between a root calo hit and all connected hits
     */
    void IdentifyDisconnectedRegions(const HitEdgeMap &hitToEdgesMap, HitConnectionsMap &graphs) const;

    /**
     *  @brief  Identify the different disconnected sub graphs in a set of eges.
     *
     *  @param  graphs the map between a root calo hit and all connected hits
     *  @param  hitToEdgesMap the (updatable) map from calo hits to edges
     */
    void ConnectRegions(const HitConnectionsMap &graphs, HitEdgeMap &hitToEdgesMap);

    /**
     *  @brief  Convert a list of calo hits into an Eigen matrix.
     *
     *  @param  caloHitList the calo hit list containing the hits from which to construct a graph
     *  @param  hitMatrix the output Eigen matrix
     */
    void Vectorize(const pandora::CaloHitList &caloHitList, Eigen::MatrixXf &hitMatrix) const;

    EdgeVector m_edges;     ///< The edges defining the graph
    bool m_fullyConnect;    ///< Whether or not to connect any disconnected regions
    int m_nSourceEdges;     ///< The number of edges to consider emerging from a source
};

} // namespace lar_content

#endif // #ifndef LAR_GRAPH_H
