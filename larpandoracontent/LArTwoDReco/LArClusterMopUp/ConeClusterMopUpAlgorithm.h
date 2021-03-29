/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/ConeClusterMopUpAlgorithm.h
 *
 *  @brief  Header file for the cone cluster mop up algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CONE_CLUSTER_MOP_UP_ALGORITHM_H
#define LAR_CONE_CLUSTER_MOP_UP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/ClusterMopUpBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ConeClusterMopUpAlgorithm class
 */
class ConeClusterMopUpAlgorithm : public ClusterMopUpBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ConeClusterMopUpAlgorithm();

private:
    void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters) const;

    typedef std::pair<float, float> Coordinate;
    typedef std::vector<Coordinate> CoordinateList;

    /**
     *  @brief  Sort coordinates by increasing transverse displacement
     *
     *  @param  lhs the first coordinate for comparison
     *  @param  rhs the second coordinate for comparison
     *
     *  @return boolean
     */
    static bool SortCoordinates(const Coordinate &lhs, const Coordinate &rhs);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_slidingFitWindow; ///< The layer window for the sliding linear fits
    float m_showerEdgeMultiplier;    ///< Artificially tune width of shower envelope so as to make it more/less inclusive
    float m_coneAngleCentile;        ///< Cluster cone angle is defined using specified centile of distribution of hit half angles
    float m_maxConeLengthMultiplier; ///< Consider hits as bound if inside cone, with projected distance less than N times cone length
    float m_minBoundedFraction;      ///< The minimum cluster bounded fraction for merging
};

} // namespace lar_content

#endif // #ifndef LAR_CONE_CLUSTER_MOP_UP_ALGORITHM_H
