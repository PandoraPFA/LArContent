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
     *  @param  upstreamMergePoint the upstream merge point
     *  @param  upstreamMergeDirection the cluster direction at the upstream merge point
     *  @param  downstreamMergePoint the downstream merge point
     *  @param  downstreamMergeDirection the cluster direction at the downstream merge point
     */
    ClusterAssociation(const pandora::CartesianVector &upstreamMergePoint, const pandora::CartesianVector &upstreamMergeDirection,
        const pandora::CartesianVector &downstreamMergePoint, const pandora::CartesianVector &downstreamMergeDirection);

    /**
     *  @brief  Returns the upstream cluster merge point
     *
     *  @return  the upstream merge point
     */
    const pandora::CartesianVector GetUpstreamMergePoint() const;

    /**
     *  @brief  Returns the cluster direction at the upstream merge point
     *
     *  @return  the cluster direction at the upstream merge point
     */
    const pandora::CartesianVector GetUpstreamMergeDirection() const;

    /**
     *  @brief  Returns the downstream cluster merge point
     *
     *  @return  the downstream merge point
     */
    const pandora::CartesianVector GetDownstreamMergePoint() const;

    /**
     *  @brief  Returns the cluster direction at the downstream merge point
     *
     *  @return  the cluster direction at the downstream merge point
     */
    const pandora::CartesianVector GetDownstreamMergeDirection() const;

    /**
     *  @brief  Returns the unit vector of the line connecting the upstream and downstream merge points (upstream -> downstream)
     *
     *  @return  the unit displacement vector from the upstream merge point to the downstream merge point
     */
    const pandora::CartesianVector GetConnectingLineDirection() const;

    /**
     *  @brief  Set the upstream merge point
     *
     *  @param  upstreamMergePoint the new upstream merge point
     */
    void SetUpstreamMergePoint(const pandora::CartesianVector &upstreamMergePoint);

    /**
     *  @brief  Set the downstream merge point
     *
     *  @param  downstreamMergePoint the new downstream merge point
     */
    void SetDownstreamMergePoint(const pandora::CartesianVector &downstreamMergePoint);

    bool operator==(const ClusterAssociation &clusterAssociation) const;
    bool operator<(const ClusterAssociation &clusterAssociation) const;

protected:
    /**
     *  @brief  Update the connecting line
     */
    void UpdateConnectingLine();

    pandora::CartesianVector m_upstreamMergePoint; ///< The upstream cluster point to be used in the merging process
    pandora::CartesianVector m_upstreamMergeDirection; ///< The upstream cluster direction at the upstream merge point (points in the direction of the downstream cluster)
    pandora::CartesianVector m_downstreamMergePoint; ///< The downstream cluster point to be used in the merging process
    pandora::CartesianVector m_downstreamMergeDirection; ///< The downstream cluster direction at the downstream merge point (points in the direction of the upstream cluster)
    pandora::CartesianVector m_connectingLineDirection; ///< The unit vector of the line connecting the upstream and downstream merge points (upstream -> downstream)
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ClusterPairAssociation class
 */
class ClusterPairAssociation : public ClusterAssociation
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterPairAssociation();

    /**
     *  @brief  Constructor
     *
     *  @param  upstreamMergePoint the upstream merge point
     *  @param  upstreamMergeDirection the cluster direction at the upstream merge point
     *  @param  downstreamMergePoint the downstream merge point
     *  @param  downstreamMergeDirection the cluster direction at the downstream merge point
     *  @param  pUpstreamCluster the address of the upstream cluster
     *  @param  pDownstreamCluster the address of the downstream cluster
     */
    ClusterPairAssociation(const pandora::CartesianVector &upstreamMergePoint, const pandora::CartesianVector &upstreamMergeDirection,
        const pandora::CartesianVector &downstreamMergePoint, const pandora::CartesianVector &downstreamMergeDirection,
        const pandora::Cluster *pUpstreamCluster, const pandora::Cluster *pDownstreamCluster);

    /**
     *  @brief  Returns the address of the upstream cluster
     *
     *  @return  the address of the upstream cluster
     */
    const pandora::Cluster *GetUpstreamCluster() const;

    /**
     *  @brief  Returns the address of the downstream cluster
     *
     *  @return  the address of the downstream cluster
     */
    const pandora::Cluster *GetDownstreamCluster() const;

private:
    const pandora::Cluster *m_pUpstreamCluster;   ///< The address of the upstream cluster
    const pandora::Cluster *m_pDownstreamCluster; ///< The address of the downstream cluster
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
    const pandora::CartesianVector connectingLineDirection(
        m_downstreamMergePoint.GetX() - m_upstreamMergePoint.GetX(), 0.f, m_downstreamMergePoint.GetZ() - m_upstreamMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ClusterAssociation::operator==(const ClusterAssociation &clusterAssociation) const
{
    return (m_upstreamMergePoint == clusterAssociation.GetUpstreamMergePoint() &&
        m_upstreamMergeDirection == clusterAssociation.GetUpstreamMergeDirection() &&
        m_downstreamMergePoint == clusterAssociation.GetDownstreamMergePoint() &&
        m_downstreamMergeDirection == clusterAssociation.GetDownstreamMergeDirection());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ClusterAssociation::operator<(const ClusterAssociation &clusterAssociation) const
{
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
    this->UpdateConnectingLine();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ClusterAssociation::SetDownstreamMergePoint(const pandora::CartesianVector &downstreamMergePoint)
{
    m_downstreamMergePoint = downstreamMergePoint;
    this->UpdateConnectingLine();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ClusterAssociation::UpdateConnectingLine()
{
    const pandora::CartesianVector connectingLineDirection(
        m_downstreamMergePoint.GetX() - m_upstreamMergePoint.GetX(), 0.f, m_downstreamMergePoint.GetZ() - m_upstreamMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterPairAssociation::ClusterPairAssociation(const pandora::CartesianVector &upstreamMergePoint,
    const pandora::CartesianVector &upstreamMergeDirection, const pandora::CartesianVector &downstreamMergePoint,
    const pandora::CartesianVector &downstreamMergeDirection, const pandora::Cluster *pUpstreamCluster, const pandora::Cluster *pDownstreamCluster) :
    ClusterAssociation(upstreamMergePoint, upstreamMergeDirection, downstreamMergePoint, downstreamMergeDirection),
    m_pUpstreamCluster(pUpstreamCluster),
    m_pDownstreamCluster(pDownstreamCluster)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterPairAssociation::ClusterPairAssociation() :
    ClusterAssociation(),
    m_pUpstreamCluster(nullptr),
    m_pDownstreamCluster(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *ClusterPairAssociation::GetUpstreamCluster() const
{
    return m_pUpstreamCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *ClusterPairAssociation::GetDownstreamCluster() const
{
    return m_pDownstreamCluster;
}

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_ASSOCIATION_H
