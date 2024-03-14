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

    /**
     *  @brief  Construct a graph from a list of calo hits.
     *
     *  @param  caloHitList the calo hit list containing the hits from which to construct a graph
     *  @param  edges the output vector of edges to populate
     */
    static void MakeGraph(const pandora::CaloHitList &caloHitList, EdgeVector &edges);
};

} // namespace lar_content

#endif // #ifndef LAR_GRAPH_HELPER_H
