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

namespace lar
{

/**
 *  @brief  ConeBasedMergingAlgorithm class
 */
class ConeBasedMergingAlgorithm : public pandora::Algorithm
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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

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
         */
        ConeParameters(pandora::Cluster *pCluster);

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
         *  @brief  Whether the cone points forwards in Z (i.e. apex at low Z)
         * 
         *  @return boolean
         */
        bool IsForwardInZ() const;

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
        pandora::Cluster           *m_pCluster;                 ///< 
        pandora::CartesianVector    m_direction;                ///< 
        pandora::CartesianVector    m_apex;                     ///< 
        pandora::CartesianVector    m_baseCentre;               ///< 
        float                       m_coneLength;               ///< 
        float                       m_coneCosHalfAngle;         ///< 
        bool                        m_isForwardInZ;             ///< 
    };

    typedef std::vector<ConeParameters> ConeParametersList;

    /**
     *  @brief  MergeParameters class
     */
    class MergeParameters
    {
    public:
        /**
         *  @brief  Default constructor
         */
        MergeParameters();

        /**
         *  @brief  Constructor
         * 
         *  @param  pDaughterCluster address of the daughter cluster
         *  @param  parentConeParameters the parent cone parameters
         */
        MergeParameters(pandora::Cluster *const pDaughterCluster, const ConeParameters &parentConeParameters);

        /**
         *  @brief  Whether the merge parameters are initialized
         * 
         *  @return boolean
         */
        bool IsInitialized() const;

        /**
         *  @brief  Get the address of the daughter cluster
         * 
         *  @return the address of the daughter cluster
         */
        pandora::Cluster *GetDaughterCluster() const;

        /**
         *  @brief  Get the address of the parent cluster
         * 
         *  @return the address of the parent cluster
         */
        pandora::Cluster *GetParentCluster() const;

        /**
         *  @brief  Get the cosine of the daughter inner vertex wrt the parent cone axis
         * 
         *  @return the cosine of the daughter inner vertex wrt the parent cone axis
         */
        float GetCosThetaInner() const;

        /**
         *  @brief  Get the cosine of the daughter outer vertex wrt the parent cone axis
         * 
         *  @return the cosine of the daughter outer vertex wrt the parent cone axis
         */
        float GetCosThetaOuter() const;

        /**
         *  @brief  Get the maximum of the inner and outer cosine theta measurements
         * 
         *  @return the maximum of the inner and outer cosine theta measurements
         */
        float GetCosThetaMax() const;

        /**
         *  @brief  Get the cosine of the daughter midpoint wrt the parent cone axis
         * 
         *  @return the cosine of the daughter midpoint wrt the parent cone axis
         */
        float GetCosThetaMidpoint() const;

        /**
         *  @brief  Get the daughter cone axis projection
         * 
         *  @return the daughter cone axis projection
         */
        float GetConeAxisProjection() const;

        /**
         *  @brief  Get the parent cluster cone length
         * 
         *  @return the parent cluster cone length
         */
        float GetParentConeLength() const;

        /**
         *  @brief  Get the cosine of the parent cluster cone half angle
         * 
         *  @return the cosine of the parent cluster cone half angle
         */
        float GetParentCosConeHalfAngle() const;

        /**
         *  @brief  operator<
         * 
         *  @param  rhs merge parameters for comparison
         */
        bool operator< (const MergeParameters &rhs) const;

    private:
        bool                        m_isInitialized;            ///< Whether the merge parameters are initialized
        pandora::Cluster           *m_pDaughterCluster;         ///< The address of the daughter cluster
        pandora::Cluster           *m_pParentCluster;           ///< The address of the daughter cluster
        float                       m_cosThetaInner;            ///< The cosine of the daughter inner vertex wrt. the parent cone axis
        float                       m_cosThetaOuter;            ///< The cosine of the daughter outer vertex wrt. the parent cone axis
        float                       m_cosThetaMidpoint;         ///< The cosine of the daughter midpoint wrt the parent cone axis
        float                       m_coneAxisProjection;       ///< The daughter cone axis projection
        float                       m_parentConeLength;         ///< The parent cluster cone length
        float                       m_parentCosConeHalfAngle;   ///< The cosine of the parent cluster cone half angle
    };

    typedef std::vector<MergeParameters> MergeParametersList;
    typedef std::map<pandora::Cluster*, MergeParametersList> ClusterMergeMap;

    /**
     *  @brief  Sort merge parameters by cos theta max
     * 
     *  @param  rhs merge parameters
     *  @param  lhs merge parameters
     */
    static bool SortByCosThetaMax(const MergeParameters &lhs, const MergeParameters &rhs);

    /**
     *  @brief  Sort merge parameters by cone axis projection
     * 
     *  @param  rhs merge parameters
     *  @param  lhs merge parameters
     */
    static bool SortByConeAxisProjection(const MergeParameters &lhs, const MergeParameters &rhs);

    std::string                     m_seedClusterListName;      ///< The seed cluster list name
    std::string                     m_nonSeedClusterListName;   ///< The non seed cluster list name
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

inline bool ConeBasedMergingAlgorithm::ConeParameters::IsForwardInZ() const
{
    return m_isForwardInZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ConeBasedMergingAlgorithm::MergeParameters::MergeParameters() :
     m_isInitialized(false),
     m_pDaughterCluster(NULL),
     m_pParentCluster(NULL),
     m_cosThetaInner(-1.f),
     m_cosThetaOuter(-1.f),
     m_cosThetaMidpoint(-1.f),
     m_coneAxisProjection(std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ConeBasedMergingAlgorithm::MergeParameters::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *ConeBasedMergingAlgorithm::MergeParameters::GetDaughterCluster() const
{
    return m_pDaughterCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *ConeBasedMergingAlgorithm::MergeParameters::GetParentCluster() const
{
    return m_pParentCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::MergeParameters::GetCosThetaInner() const
{
    return m_cosThetaInner;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::MergeParameters::GetCosThetaOuter() const
{
    return m_cosThetaOuter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::MergeParameters::GetCosThetaMax() const
{
    return std::max(m_cosThetaInner, m_cosThetaOuter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::MergeParameters::GetCosThetaMidpoint() const
{
    return m_cosThetaMidpoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::MergeParameters::GetConeAxisProjection() const
{
    return m_coneAxisProjection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::MergeParameters::GetParentConeLength() const
{
    return m_parentConeLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ConeBasedMergingAlgorithm::MergeParameters::GetParentCosConeHalfAngle() const
{
    return m_parentCosConeHalfAngle;
}

} // namespace lar

#endif // #ifndef LAR_CONE_BASED_MERGING_ALGORITHM_H
