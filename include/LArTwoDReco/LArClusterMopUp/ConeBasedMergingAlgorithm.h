/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.h
 * 
 *  @brief  Header file for the cone based merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CONE_BASED_MERGING_ALGORITHM_H
#define LAR_CONE_BASED_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterMopUp/ClusterMopUpAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ConeBasedMergingAlgorithm class
 */
class ConeBasedMergingAlgorithm : public ClusterMopUpAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    ConeBasedMergingAlgorithm();

private:
    /**
     *  @brief  ConeParameters class
     */
    class ConeParameters
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster address of the cluster
         *  @param  slidingFitWindow the layer window to use in sliding shower fit for cone construction
         *  @param  slidingFitLayerPitch the layer pitch, in cm, to use in sliding shower fit for cone construction
         *  @param  coneAngleCentile the cone angle centile
         */
        ConeParameters(pandora::Cluster *pCluster, const unsigned int slidingFitWindow, const float slidingFitLayerPitch, const float coneAngleCentile);

        /**
         *  @brief  Get the address of the cluster
         * 
         *  @return address of the cluster
         */
        pandora::Cluster *GetCluster() const;

        /**
         *  @brief  Get the cone direction
         * 
         *  @return the cone direction
         */
        const pandora::CartesianVector &GetDirection() const;

        /**
         *  @brief  Get the position vector of the cone apex
         * 
         *  @return the position vector of the cone apex
         */
        const pandora::CartesianVector &GetApex() const;

        /**
         *  @brief  Get the position vector of the centre of the cone base
         * 
         *  @return the position vector of the centre of the cone base
         */
        const pandora::CartesianVector &GetBaseCentre() const;

        /**
         *  @brief  Get the cone length
         * 
         *  @return the cone length
         */
        float GetConeLength() const;

        /**
         *  @brief  Get the cone cos half angle
         * 
         *  @return the cone cos half angle
         */
        float GetConeCosHalfAngle() const;

        /**
         *  @brief  Whether the cone points forward wrt to increasing sliding fit layers (i.e. apex at min layer)
         * 
         *  @return boolean
         */
        bool IsForward() const;

        /**
         *  @brief  Get the cos half angle for a specified position
         * 
         *  @param  position the specified position
         * 
         *  @return the position cos half angle
         */
        float GetPositionCosHalfAngle(const pandora::CartesianVector &position) const;

        /**
         *  @brief  Get the projection of a specified position onto the cone axis
         * 
         *  @param  position the specified position
         * 
         *  @return the projection onto the cone axis
         */
        float GetConeAxisProjection(const pandora::CartesianVector &position) const;

    private:
        /**
         *  @brief  Get the cone direction estimate, using a provided sliding linear fit
         * 
         *  @param  fitResult the sliding linear fit result
         *  @param  isForward to receive the estimate of whether cone points forward wrt to increasing sliding fit layers
         *  @param  direction to receive the direction estimate
         */
        void GetDirectionEstimate(const TwoDSlidingFitResult &fitResult, bool &isForward, pandora::CartesianVector &direction) const;

        /**
         *  @brief  Get the cone cos half angle estimate for a specified cluster
         * 
         *  @param  position the specified position
         *  @param  direction the cone direction
         *  @param  apex the cone apex
         *  @param  coneAngleCentile the cone angle centile
         * 
         *  @return the cone cos half angle estimate
         */
        float GetCosHalfAngleEstimate(const pandora::Cluster *const pCluster, const pandora::CartesianVector &direction,
            const pandora::CartesianVector &apex, const float coneAngleCentile) const;

        typedef std::map<int, unsigned int> HitsPerLayerMap;

        pandora::Cluster           *m_pCluster;             ///< The address of the cluster
        pandora::CartesianVector    m_direction;            ///< The cone direction
        pandora::CartesianVector    m_apex;                 ///< The position vector of the cone apex
        pandora::CartesianVector    m_baseCentre;           ///< The position vector of the centre of the cone base
        float                       m_coneLength;           ///< The cone length
        float                       m_coneCosHalfAngle;     ///< The cone cos half angle
        bool                        m_isForward;            ///< The cone points forward wrt to increasing sliding fit layers (i.e. apex at min layer)
    };

    void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters, const ClusterToListNameMap &clusterToListNameMap) const;

    /**
     *  @brief  Get the fraction of hits in a cluster bounded by a specified cluster cone
     * 
     *  @param  pCluster address of the cluster
     *  @param  coneParameters the cone parameters
     * 
     *  @return the fraction of bounded hits
     */
    float GetBoundedFraction(const pandora::Cluster *const pCluster, const ConeParameters &coneParameters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_slidingFitWindow;         ///< The layer window for the sliding linear fits
    float           m_coneAngleCentile;         ///< Cluster cone angle is defined using specified centile of distribution of hit cos half angles
    float           m_maxConeLengthMultiplier;  ///< Consider hits as bound if inside cone, with projected distance less than N times cone length
    float           m_minBoundedFraction;       ///< The minimum cluster bounded fraction for merging
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ConeBasedMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ConeBasedMergingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *ConeBasedMergingAlgorithm::ConeParameters::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ConeBasedMergingAlgorithm::ConeParameters::GetDirection() const
{
    return m_direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ConeBasedMergingAlgorithm::ConeParameters::GetApex() const
{
    return m_apex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ConeBasedMergingAlgorithm::ConeParameters::GetBaseCentre() const
{
    return m_baseCentre;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::ConeParameters::GetConeLength() const
{
    return m_coneLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::ConeParameters::GetConeCosHalfAngle() const
{
    return m_coneCosHalfAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ConeBasedMergingAlgorithm::ConeParameters::IsForward() const
{
    return m_isForward;
}

} // namespace lar_content

#endif // #ifndef LAR_CONE_BASED_MERGING_ALGORITHM_H
