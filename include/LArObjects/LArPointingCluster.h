/**
 *  @file   LArContent/include/LArObjects/LArPointingCluster.h
 * 
 *  @brief  Header file for the lar pointing cluster class.
 * 
 *  $Log: $
 */
#ifndef LAR_POINTING_CLUSTER_H
#define LAR_POINTING_CLUSTER_H 1

#include "LArObjects/LArTwoDSlidingFitResult.h"
#include "LArObjects/LArThreeDSlidingFitResult.h"

namespace lar_content
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
         *  @brief  Default constructor
         */
        Vertex();

        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster address of the cluster
         *  @param  position the vertex position
         *  @param  direction the vertex direction
         *  @param  rms the rms from vertex fit
         *  @param  isInner whether this is a cluster inner or outer vertex
         */
        Vertex(const pandora::Cluster *const pCluster, const pandora::CartesianVector &position, const pandora::CartesianVector &direction,
            const float rms, const bool isInner);

        /**
         *  @brief  Copy constructor
         * 
         *  @param  rhs the vertex instance to copy
         */
        Vertex(const Vertex &rhs);

        /**
         *  @brief  Destructor
         * 
         */
        ~Vertex();

        /**
         *  @brief  Get the address of the cluster
         * 
         *  @return the address of the cluster
         */
        const pandora::Cluster *GetCluster() const;

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
         *  @brief  Whether the vertex has been initialized
         * 
         *  @return boolean
         */
        bool IsInitialized() const;

        /**
         *  @brief  Vertex assigment operator
         * 
         *  @param  rhs the vertex to assign
         */
        Vertex &operator=(const Vertex &rhs);

    private:
        const pandora::Cluster     *m_pCluster;             ///< The address of the cluster
        pandora::CartesianVector    m_position;             ///< The vertex position
        pandora::CartesianVector    m_direction;            ///< The vertex direction
        float                       m_rms;                  ///< Rms from vertex fit
        bool                        m_isInner;              ///< Whether this is the inner vertex
        bool                        m_isInitialized;        ///< Whether the vertex has been initialized
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  pCluster address of the cluster
     *  @param  fitHalfLayerWindow the fit layer half window
     *  @param  fitLayerPitch the fit layer pitch, units cm
     */
    LArPointingCluster(const pandora::Cluster *const pCluster, const unsigned int fitHalfLayerWindow = 10, const float fitLayerPitch = 0.3f);

    /**
     *  @brief  Constructor
     * 
     *  @param slidingFitResult the input sliding fit result
     */
    LArPointingCluster(const TwoDSlidingFitResult &slidingFitResult);   

    /**
     *  @brief  Constructor
     * 
     *  @param slidingFitResult the input sliding fit result
     */
    LArPointingCluster(const ThreeDSlidingFitResult &slidingFitResult);

    /**
     *  @brief  Get the address of the cluster
     * 
     *  @return the address of the cluster
     */
    const pandora::Cluster *GetCluster() const;

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
 
    /**
     *  @brief Build the pointing cluster object from the sliding fit result
     * 
     *  @param slidingFitResult the input sliding fit result
     */
    void BuildPointingCluster(const TwoDSlidingFitResult &slidingFitResult);  

    /**
     *  @brief Build the pointing cluster object from the sliding fit result
     * 
     *  @param slidingFitResult the input sliding fit result
     */
    void BuildPointingCluster(const ThreeDSlidingFitResult &slidingFitResult);

    const pandora::Cluster         *m_pCluster;             ///< The address of the cluster
    Vertex                          m_innerVertex;          ///< The inner vertex
    Vertex                          m_outerVertex;          ///< The outer vertex
};

typedef std::vector<LArPointingCluster> LArPointingClusterList;
typedef std::vector<LArPointingCluster::Vertex> LArPointingClusterVertexList;
typedef std::map<const pandora::Cluster*, LArPointingCluster> LArPointingClusterMap;

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *LArPointingCluster::GetCluster() const
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

inline const pandora::Cluster *LArPointingCluster::Vertex::GetCluster() const
{
    if (!m_isInitialized)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArPointingCluster::Vertex::GetPosition() const
{
    if (!m_isInitialized)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArPointingCluster::Vertex::GetDirection() const
{
    if (!m_isInitialized)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArPointingCluster::Vertex::GetRms() const
{
    if (!m_isInitialized)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArPointingCluster::Vertex::IsInnerVertex() const
{
    if (!m_isInitialized)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_isInner;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArPointingCluster::Vertex::IsInitialized() const
{
    return m_isInitialized;
}

} // namespace lar_content

#endif // #ifndef LAR_POINTING_CLUSTER_H
