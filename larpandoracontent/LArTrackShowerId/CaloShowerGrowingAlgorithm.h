/**
 *  @file   larpandoracontent/LArTrackShowerId/CaloShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the calorimetric shower growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CALO_SHOWER_GROWING_ALGORITHM_H
#define LAR_CALO_SHOWER_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CaloShowerGrowingAlgorithm class
 */
class CaloShowerGrowingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CaloShowerGrowingAlgorithm();

private:
    class Bounds
    {
    public:
        /**
         *  Constructor. It doesn't really matter if the vertices correspond to top-left etc, as long as the coordinates are specified in a
         *  way that makes them cycle in either a counter-clockwise or clockwise direction about the chain tl -> bl -> br -> tr -> tl
         *
         *  @param  tl top left vertex
         *  @param  bl bottom left vertex
         *  @param  br bottom right bertex
         *  @param  tr top right vertex
         *  @param  pAlgorithm the parent algorithm
         */
        Bounds(const pandora::CartesianVector &tl, const pandora::CartesianVector &bl, const pandora::CartesianVector &br, const pandora::CartesianVector &tr);

        /*  @brief  Determines if a point is contained within a the bounding region
         *
         *  @return true if the point is within the bounds, false otherwise
         */
        bool Contains(const pandora::CartesianVector &point) const;

        const pandora::CartesianVector m_a;     ///< The top-left corner of the bounding region
        const pandora::CartesianVector m_b;     ///< The bottom-left corner of the bounding region
        const pandora::CartesianVector m_c;     ///< The bottom-right corner of the bounding region
        const pandora::CartesianVector m_d;     ///< The top-right corner of the bounding region

    private:
        /*  @brief  Determines if a point is contained within a triangle using the Barycentric technique, as described in
         *          Real-Time Collision Detection, C. Ericson (Morgan Kaufmann, 2005).
         *
         *  @param  ab the vector describing the edge ab of the triangle abc
         *  @param  ac the vector describing the edge ac of the triangle abc
         *  @param  ap the vector describing the displacement of a point p from vertex a of the triangle abc
         *
         *  @return true if the point p is within the bounds of the triangle abc, false otherwise
         */
        bool Contains(const pandora::CartesianVector &ab, const pandora::CartesianVector &ac, const pandora::CartesianVector &ap) const;

        // ATTN: The edges ab and ac describe one triangle of the bounding region, ac and ad describe the other
        const pandora::CartesianVector m_ab;    ///< Vector from tl to bl
        const pandora::CartesianVector m_ac;    ///< Vector from tl to br
        const pandora::CartesianVector m_ad;    ///< Vector from tl to tr
    };

    typedef std::map<const pandora::Cluster *, pandora::ClusterList> ClusterAssociationMap;

    pandora::StatusCode Run();

    /**
     *  @brief  Grow the showers in a cluster list
     *
     *  @param  clusterList the list of available seed candidates
     */
    void GrowShowers(const pandora::ClusterList &clusterList) const;

    /**
     *  @brief  Get the seed clusters
     *
     *  @param  clusterList the list of available seed candidates
     *  @param  seedClusterList the output list of seed clusters
     */
    void GetSeedClusters(const pandora::ClusterList &clusterList, pandora::ClusterList &seedClusterList) const;

    /**
     *  @brief  Get the bounding box for a seed cluster
     *
     *  @param  pSeed the seed whose bounds should be calculated
     *  
     *  @return the bounding region
     */
    Bounds GetSeedBounds(const pandora::Cluster *pSeed) const;

    /**
     *  @brief  Assess the plausibility of the group of clusters as representing a shower
     *
     *  @param  pSeed the seed around which the shower is to be formed
     *  @param  associatedClusterList the list of clusters associated with the seed
     *  @param  showerClusterList the output list of clusters definig a shower (empty if no plausible candidate found)
     *
     *  @return the chi2 for the best association
     */
    float AssessAssociation(const pandora::Cluster *pSeed, const pandora::ClusterList &associatedClusterlist, pandora::ClusterList &showerClusterList) const;

    /**
     *  @brief  Retrieve a reduced chi2 value for a fit of the given set of calo hits to photon and electron longitudinal energy profiles.
     *          Considers the forward and backward energy profiles for the hits against both photon and electron profiles and returns the
     *          reduced chi2 of the best match
     *
     *  @param  caloHitList the list of calo hits for which energy profiles are to be computed
     *  @return the reduced chi2 of the best profile match
     */
    float GetShowerProfileChi2(const pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Get the origin and unit direction of a PCA fit to the specified list of calo hits
     *
     *  @param  caloHitList the list of calo hits from which to extract a projection axis
     *  @param  origin the output origin of the axis
     *  @param  dir the output direction of the axis
     */
    void GetProjectionAxis(const pandora::CaloHitList &caloHitList, pandora::CartesianVector &origin, pandora::CartesianVector &dir) const;

    /**
     *  @brief  Get the expected longitudinal energy profile for a photon
     *
     *  @param  e0 the energy of the photon (MeV)
     *  @param  positions the array of longitudinal positions (interaction lengths)
     *  @param  binSize the size of the position bins (interaction lengths)
     *  @param  fractionalEnergies the output fractional energies (same length as tt)
     */
    void GetPhotonLongitudinalEnergyProfile(const float e0, const pandora::FloatVector &positions, const float binSize, pandora::FloatVector &fractionalEnergies) const;

    /**
     *  @brief  Get the expected longitudinal energy profile for an electron
     *
     *  @param  e0 the energy of the electron (MeV)
     *  @param  positions the array of longitudinal positions (interaction lengths)
     *  @param  binSize the size of the position bins (interaction lengths)
     *  @param  fractionalEnergies the output fractional energies (same length as tt)
     */
    void GetElectronLongitudinalEnergyProfile(const float e0, const pandora::FloatVector &positions, const float binSize, pandora::FloatVector &fractionalEnergies) const;

    /**
     *  @brief  Get the expected longitudinal energy profile for a particle
     *
     *  @param  e0 the energy of the particle (MeV)
     *  @param  positions the array of longitudinal positions (interaction lengths)
     *  @param  cj length scale adjustment factor, -0.5 for electrons, +0.5 for photons
     *  @param  binSize the size of the position bins (interaction lengths)
     *  @param  fractionalEnergies the output fractional energies (same length as tt)
     */
    void GetLongitudinalEnergyProfile(const float e0, const pandora::FloatVector &positions, const float cj, const float binSize, pandora::FloatVector &fractionalEnergies) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputClusterListNames; ///< The names of the input cluster lists
    unsigned int m_minCaloHitsForSeed; ///< The minimum number of calo hits per seed cluster
    float m_radiationLength; ///< The radiation length in liquid argon
    float m_moliereRadius; ///< The Moliere radius in liquid argon
    bool m_visualize; ///< Whether or not to visualize the algorithm steps
};

} // namespace lar_content

#endif // #ifndef LAR_CALO_SHOWER_GROWING_ALGORITHM_H
