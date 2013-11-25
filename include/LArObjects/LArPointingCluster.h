/**
 *  @file   LArContent/include/LArObjects/LArPointingCluster.h
 * 
 *  @brief  Header file for the lar pointing cluster class.
 * 
 *  $Log: $
 */
#ifndef LAR_POINTING_CLUSTER_H
#define LAR_POINTING_CLUSTER_H 1

namespace pandora {class Cluster;};

namespace lar
{

/**
 *  @brief  LArPointingCluster class
 */
class LArPointingCluster
{
public:
    /**
     *  @brief  Vertex class
     */
    class Vertex
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster address of the cluster
         *  @param  useInnerVertex whether to use cluster inner or outer vertex
         *  @param  nLayersToSkip the number of layers to skip
         */
        Vertex(pandora::Cluster *const pCluster, const bool useInnerVertex, const unsigned int nLayersToSkip);

        /**
         *  @brief  Get the address of the cluster
         * 
         *  @return the address of the cluster
         */
        pandora::Cluster *GetCluster() const;

        /**
         *  @brief  Get the vertex position
         * 
         *  @return the vertex position
         */
        const pandora::CartesianVector &GetPosition() const;

        /**
         *  @brief  Get the vertex direction
         * 
         *  @return the vertex direction
         */
        const pandora::CartesianVector &GetDirection() const;

         /**
         *  @brief  Get rms from vertex fit
         * 
         *  @return the rms
         */
        float GetRms() const;

        /**
         *  @brief  Is this the inner vertex
         * 
         *  @return boolean
         */
        bool IsInnerVertex() const;

        /**
         *  @brief  Get the number of skipped layers
         * 
         *  @return the number of skipped layers
         */
        unsigned int GetNSkippedLayers() const;

    private:
        pandora::Cluster           *m_pCluster;             ///< The address of the cluster
        pandora::CartesianVector    m_position;             ///< The vertex position
        pandora::CartesianVector    m_direction;            ///< The vertex direction
        float                       m_rms;                  ///< rms from vertex fit
        bool                        m_isInner;              ///< The inner vertex
        unsigned int                m_nSkippedLayers;       ///< The number of skipped layers
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  pCluster address of the cluster
     */
    LArPointingCluster(pandora::Cluster *const pCluster);

    /**
     *  @brief  Get the address of the cluster
     * 
     *  @return the address of the cluster
     */
    pandora::Cluster *GetCluster() const;

    /**
     *  @brief  Get the inner vertex
     * 
     *  @return the inner vertex
     */
    const Vertex &GetInnerVertex() const;

    /**
     *  @brief  Get the outer vertex
     * 
     *  @return the outer vertex
     */
    const Vertex &GetOuterVertex() const;

    /**
     *  @brief  Get length squared of pointing cluster
     * 
     *  @return the length squared
     */
    float GetLengthSquared() const;

    /**
     *  @brief  Get length of pointing cluster
     * 
     *  @return the length
     */
    float GetLength() const;

private:
    pandora::Cluster               *m_pCluster;             ///< The address of the cluster
    Vertex                          m_innerVertex;          ///< The inner vertex
    Vertex                          m_outerVertex;          ///< The outer vertex
};

typedef std::vector<LArPointingCluster> LArPointingClusterList;
typedef std::vector<LArPointingCluster::Vertex> LArPointingClusterVertexList;
typedef std::map<const pandora::Cluster*, LArPointingCluster> LArPointingClusterMap;

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArPointingCluster::LArPointingCluster(pandora::Cluster *const pCluster) :
    m_pCluster(pCluster),
    m_innerVertex(pCluster, true, 0),
    m_outerVertex(pCluster, false, 0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *LArPointingCluster::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArPointingCluster::Vertex &LArPointingCluster::GetInnerVertex() const
{
    return m_innerVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArPointingCluster::Vertex &LArPointingCluster::GetOuterVertex() const
{
    return m_outerVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArPointingCluster::GetLengthSquared() const
{
    return (m_outerVertex.GetPosition() - m_innerVertex.GetPosition()).GetMagnitudeSquared();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArPointingCluster::GetLength() const
{
    return std::sqrt(this->GetLengthSquared());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *LArPointingCluster::Vertex::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArPointingCluster::Vertex::GetPosition() const
{
    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArPointingCluster::Vertex::GetDirection() const
{
    return m_direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArPointingCluster::Vertex::GetRms() const
{
    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArPointingCluster::Vertex::IsInnerVertex() const
{
    return m_isInner;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArPointingCluster::Vertex::GetNSkippedLayers() const
{
    return m_nSkippedLayers;
}

} // namespace lar

#endif // #ifndef LAR_POINTING_CLUSTER_H
