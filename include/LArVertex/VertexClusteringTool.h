/**
 *  @file   ExampleContent/include/VertexClusteringTools/VertexClusteringTool.h
 * 
 *  @brief  Header file for the example algorithm tool class.
 * 
 *  $Log: $
 */
#ifndef VERTEX_CLUSTERING_TOOL_H
#define VERTEX_CLUSTERING_TOOL_H 1

namespace lar_content
{

/**
 *  @brief  VertexClusteringTool class
 */
class VertexClusteringTool : public pandora::AlgorithmTool
{

public:

    /**
     *  @brief  Default constructor
     */
    VertexClusteringTool();

    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    class VertexCluster
    {
        public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        VertexCluster();

        /**
         *  @brief  Get the address of the vertex
         * 
         *  @return the address of the vertex
         */         
        const pandora::VertexList &GetVertexList() const;

        /**
         *  @brief  Get the score
         * 
         *  @return the score
         */
        pandora::CartesianVector GetCentroidPosition() const;

        /**
         *  @brief  Get the score
         * 
         *  @return the score
         */
        void AddVertex(const pandora::Vertex *const pVertex);

        /**
         *  @brief  Get the score
         * 
         *  @return the score
         */
         
         void AddVerticesFromCluster(const VertexCluster &targetVertexCluster);
         
         /**
         *  @brief  Get the score
         * 
         *  @return the score
         */
         
         void ClearVertexCluster();
         
         /**
         *  @brief  Get the score
         * 
         *  @return the score
         */
         
        void MergeVertexClusters(const VertexCluster &targetVertexCluster);

    private:
        pandora::VertexList m_vertexList;
    };

    typedef std::vector<VertexCluster*> VertexClusterList;

    /**
     *  @brief  Clusters vertices into vertex clusters
     * 
     *  @param  pVertexList the address of the vertex list to be divided into vertex clusters
     */
    std::vector<const pandora::VertexList*> ClusterVertices(const pandora::VertexList* pVertexList);
    
    /**
     *  @brief  Checks the distance from one vertex to the closest vertex in a vertex cluster
     * 
     *  @param  pVertex the address of the vertex
     *  @param  pVertexCluster the address of the vertex cluster
     */
    bool CheckVertexToClusterDistance(const pandora::Vertex *const pVertex, VertexCluster *const pVertexCluster) const;

    /**
     *  @brief  Boolean intended to sue with std::sort in order to sort a vertex list by the Z coordinates of its vertices
     * 
     *  @param  pLhs the address of the first vertex
     *  @param  pRhs the address of the second vertex
     */
    static bool SortVerticesByZ(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs);

private:

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_clusteringRadius;                 ///< The distance between individual vertex candidates that guarantees them being in the same cluster
    float           m_maxVertexToCentroidDistance;      ///< The radius around the centroid within which vertices will be added to the current vertex cluster
    unsigned int    m_minClusterSize;                   ///< Minimum size of cluster before switching to centroid as part of the distance clustering measure
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline VertexClusteringTool::VertexCluster::VertexCluster()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::VertexList &VertexClusteringTool::VertexCluster::GetVertexList() const
{
    return m_vertexList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VertexClusteringTool::VertexCluster::AddVertex(const pandora::Vertex *const pVertex)
{
    (void) m_vertexList.insert(pVertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VertexClusteringTool::VertexCluster::AddVerticesFromCluster(const VertexCluster &targetVertexCluster)
{
    m_vertexList.insert(targetVertexCluster.GetVertexList().begin(), targetVertexCluster.GetVertexList().end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VertexClusteringTool::VertexCluster::ClearVertexCluster()
{
    m_vertexList.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *VertexClusteringTool::Factory::CreateAlgorithmTool() const
{
    return new VertexClusteringTool();
}

} // namespace lar_content

#endif // #ifndef EXAMPLE_ALGORITHM_TOOL_H
