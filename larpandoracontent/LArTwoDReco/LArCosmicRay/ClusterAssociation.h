/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/ClusterAssociation.h
 *
 *  @brief  Header file for the lar cluster association class
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_ASSOCIATION_H
#define LAR_CLUSTER_ASSOCIATION_H 1

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

namespace lar_content
{
/**
 *  @brief  ClusterAssociation class
 */
class ClusterAssociation
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterAssociation();

    /**
     *  @brief  Constructor
     *
     *  @param  pUpstreamCluster the upstream cluster of the two associated clusters
     *  @param  pDownstreamCluster the downstream cluster of the two associated clusters
     *  @param  upstreamMergePoint the upstream cluster point to be used in the merging process
     *  @param  upstreamMergeDirection the upstream cluster direction at the upstream merge point
     *  @param  downstreamMergePoint the downstream cluster point to be used in the merging process
     *  @param  downstreamMergeDirection the downstream cluster direction at the downstream merge point
     */
     ClusterAssociation(const pandora::CartesianVector &upstreamMergePoint, const pandora::CartesianVector &upstreamMergeDirection,
        const pandora::CartesianVector &downstreamMergePoint, const pandora::CartesianVector &downstreamMergeDirection);

    /**
     *  @brief  Returns the upstream cluster merge point
     *
     *  @return  CartesianVector the merge point of the upstream cluster
     */
    const pandora::CartesianVector GetUpstreamMergePoint() const;

    /**
     *  @brief  Returns the upstream cluster direction at the upstream merge point
     *
     *  @return  CartesianVector the direction at the merge point of the upstream cluster
     */        
    const pandora::CartesianVector GetUpstreamMergeDirection() const;

    /**
     *  @brief  Returns the downstream cluster merge point
     *
     *  @return  CartesianVector the merge point of the downstream cluster
     */
    const pandora::CartesianVector GetDownstreamMergePoint() const;

    /**
     *  @brief  Returns the downstream cluster direction at the downstream merge point
     *
     *  @return  CartesianVector the direction at the merge point of the downstream cluster
     */        
     const pandora::CartesianVector GetDownstreamMergeDirection() const;

    /**
     *  @brief  Returns the unit vector of the line connecting the upstream and downstream merge points (upstream -> downstream)
     *
     *  @return  CartesianVector the unit displacement vector from the upstream merge point to the downstream merge point
     */           
    const pandora::CartesianVector GetConnectingLineDirection() const;

    void SetUpstreamMergePoint(const pandora::CartesianVector &upstreamMergePoint);

    void SetDownstreamMergePoint(const pandora::CartesianVector &downstreamMergePoint);

    bool operator==(const ClusterAssociation &clusterAssociation) const;
    bool operator<(const ClusterAssociation &clusterAssociation) const;
        
protected:
    pandora::CartesianVector    m_upstreamMergePoint;          ///< The upstream cluster point to be used in the merging process
    pandora::CartesianVector    m_upstreamMergeDirection;      ///< The upstream cluster direction at the upstream merge point (points in the direction of the downstream cluster)
    pandora::CartesianVector    m_downstreamMergePoint;        ///< The downstream cluster point to be used in the merging process
    pandora::CartesianVector    m_downstreamMergeDirection;    ///< The downstream cluster direction at the downstream merge point (points in the direction of the upstream cluster)
    pandora::CartesianVector    m_connectingLineDirection;     ///< The unit vector of the line connecting the upstream and downstream merge points (upstream -> downstream)
};

class ClusterEndpointAssociation : public ClusterAssociation
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterEndpointAssociation();
    
    ClusterEndpointAssociation(const pandora::CartesianVector &upstreamMergePoint, const pandora::CartesianVector &upstreamMergeDirection,
        const pandora::CartesianVector &downstreamMergePoint, const pandora::CartesianVector &downstreamMergeDirection, const pandora::Cluster *pMainTrackCluster, const bool isEndUpstream);

    /**
     *  @brief  Returns the upstream cluster address
     *
     *  @return  Cluster the address of the upstream cluster
     */
    const pandora::Cluster *GetMainTrackCluster() const;

    void SetMainTrackCluster(const pandora::Cluster *const pMainTrackCluster);
       
    bool IsEndUpstream() const;        


private:
    const pandora::Cluster    *m_pMainTrackCluster;
    bool                       m_isEndUpstream;
};

//------------------------------------------------------------------------------------------------------------------------------------------
    
inline ClusterAssociation::ClusterAssociation() :
    m_upstreamMergePoint(pandora::CartesianVector(0.f, 0.f, 0.f)),
    m_upstreamMergeDirection(pandora::CartesianVector(0.f, 0.f, 0.f)),
    m_downstreamMergePoint(pandora::CartesianVector(0.f, 0.f, 0.f)),
    m_downstreamMergeDirection(pandora::CartesianVector(0.f, 0.f, 0.f)),
    m_connectingLineDirection(pandora::CartesianVector(0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterAssociation::ClusterAssociation(const pandora::CartesianVector &upstreamMergePoint, const pandora::CartesianVector &upstreamMergeDirection,
        const pandora::CartesianVector &downstreamMergePoint, const pandora::CartesianVector &downstreamMergeDirection) :
    m_upstreamMergePoint(upstreamMergePoint),
    m_upstreamMergeDirection(upstreamMergeDirection),
    m_downstreamMergePoint(downstreamMergePoint),
    m_downstreamMergeDirection(downstreamMergeDirection),
    m_connectingLineDirection(0.f, 0.f, 0.f)
{
    const pandora::CartesianVector connectingLineDirection(m_downstreamMergePoint.GetX() - m_upstreamMergePoint.GetX(), 0.f, m_downstreamMergePoint.GetZ() - m_upstreamMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ClusterAssociation::operator==(const ClusterAssociation &clusterAssociation) const
{
    // ISOBEL: ALSO CHECK DOWNSTREAM
    return (m_upstreamMergePoint == clusterAssociation.GetUpstreamMergePoint() &&  m_upstreamMergePoint == clusterAssociation.GetUpstreamMergePoint());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ClusterAssociation::operator<(const ClusterAssociation &clusterAssociation) const
{
    // ISOBEL: ALSO CHECK DOWNSTREAM
    return (LArClusterHelper::SortCoordinatesByPosition(m_upstreamMergePoint, clusterAssociation.GetUpstreamMergePoint()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector ClusterAssociation::GetUpstreamMergePoint() const
{
    return m_upstreamMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector ClusterAssociation::GetUpstreamMergeDirection() const
{
    return m_upstreamMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector ClusterAssociation::GetDownstreamMergePoint() const
{
    return m_downstreamMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector ClusterAssociation::GetDownstreamMergeDirection() const
{
    return m_downstreamMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector ClusterAssociation::GetConnectingLineDirection() const
{
    return m_connectingLineDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ClusterAssociation::SetUpstreamMergePoint(const pandora::CartesianVector &upstreamMergePoint)
{
    m_upstreamMergePoint = upstreamMergePoint;
    
    const pandora::CartesianVector connectingLineDirection(m_downstreamMergePoint.GetX() - m_upstreamMergePoint.GetX(), 0.f, m_downstreamMergePoint.GetZ() - m_upstreamMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ClusterAssociation::SetDownstreamMergePoint(const pandora::CartesianVector &downstreamMergePoint)
{
    m_downstreamMergePoint = downstreamMergePoint;
    
    const pandora::CartesianVector connectingLineDirection(m_downstreamMergePoint.GetX() - m_upstreamMergePoint.GetX(), 0.f, m_downstreamMergePoint.GetZ() - m_upstreamMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterEndpointAssociation::ClusterEndpointAssociation(const pandora::CartesianVector &upstreamMergePoint, const pandora::CartesianVector &upstreamMergeDirection,
    const pandora::CartesianVector &downstreamMergePoint, const pandora::CartesianVector &downstreamMergeDirection, const pandora::Cluster *pMainTrackCluster, const bool isEndUpstream) :
        ClusterAssociation(upstreamMergePoint, upstreamMergeDirection, downstreamMergePoint, downstreamMergeDirection),
        m_pMainTrackCluster(pMainTrackCluster),
        m_isEndUpstream(isEndUpstream)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterEndpointAssociation::ClusterEndpointAssociation() :
        ClusterAssociation(),
        m_pMainTrackCluster(nullptr),
        m_isEndUpstream(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *ClusterEndpointAssociation::GetMainTrackCluster() const
{
    return m_pMainTrackCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

inline bool ClusterEndpointAssociation::IsEndUpstream() const
{
    return m_isEndUpstream;
}

//------------------------------------------------------------------------------------------------------------------------------------------        

inline void ClusterEndpointAssociation::SetMainTrackCluster(const pandora::Cluster *const pMainTrackCluster)
{
    m_pMainTrackCluster = pMainTrackCluster;
}            
    
} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_ASSOCIATION_H
