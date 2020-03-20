/**
 *  @file   larpandoracontent/LArHelpers/LArHitWidthHelper.h
 *
 *  @brief  Header file for the lar hit width helper class.
 *
 *  $Log: $
 */

#ifndef LAR_HIT_WIDTH_HELPER_H
#define LAR_HIT_WIDTH_HELPER_H 1


#include "Objects/Cluster.h"

namespace lar_content
{
  

/**
 *  @brief  LArHitWidthHelper class
 */
class LArHitWidthHelper
{
public:

    class ConstituentHit
    {
    public:
        ConstituentHit(const pandora::CartesianVector &positionVector, const float hitWidth, const pandora::Cluster *const parentClusterAddress);

        pandora::CartesianVector GetPositionVector() const;
        float GetHitWidth() const;
        const pandora::Cluster* GetParentClusterAddress() const;

        struct  SortByDistanceToPoint
        {
	        SortByDistanceToPoint(const pandora::CartesianVector referencePoint) : m_referencePoint(referencePoint) {}
            bool operator() (const ConstituentHit &lhs, const ConstituentHit &rhs);
        
            const pandora::CartesianVector m_referencePoint;
        };

        
    private:
        pandora::CartesianVector m_positionVector;
        float m_hitWidth;
        const pandora::Cluster *m_parentClusterAddress; 
    };

    typedef std::vector<ConstituentHit> ConstituentHitVector;


   class ClusterParameters
    {
    public:
        ClusterParameters(const pandora::Cluster *const pCluster, const float maxConsituentHitWidth, const bool isUniformHits, const float hitWidthScalingFactor);
        ClusterParameters(const pandora::Cluster *const pCluster, const unsigned int numCaloHits, const float totalWeight, const ConstituentHitVector &constituentHitVector, const pandora::CartesianVector &lowerXExtrema, const pandora::CartesianVector &higherXExtrema);

        const pandora::Cluster* GetClusterAddress() const;
        unsigned int GetNumCaloHits() const;
        float GetTotalWeight() const;
        ConstituentHitVector GetConstituentHitVector() const;
        pandora::CartesianVector GetLowerXExtrema() const;
        pandora::CartesianVector GetHigherXExtrema() const;

    private:
        const pandora::Cluster *m_pCluster; 
        unsigned int m_numCaloHits;
        ConstituentHitVector m_constituentHitVector;
        float m_totalWeight;
        pandora::CartesianVector m_lowerXExtrema;
        pandora::CartesianVector m_higherXExtrema;
    };

    typedef std::map<const pandora::Cluster*, const ClusterParameters> ClusterToParametersMap;


    class ClusterToParametersMapStore
    {

    public:
        static ClusterToParametersMapStore* Instance();
        ClusterToParametersMap* GetMap();

    private:
        static ClusterToParametersMapStore* m_instance;
        ClusterToParametersMap m_clusterToParametersMap;
    };
    

    static ConstituentHitVector GetConstituentHits(const pandora::Cluster *const pCluster, const float maxConstituentHitWidth, const float hitWidthScalingFactor);
    static ConstituentHitVector GetUniformConstituentHits(const pandora::Cluster *const pCluster, const float constituentHitWidth, const float hitWidthScalingFactor);
    static pandora::CartesianPointVector GetConstituentHitPositionVector(const ConstituentHitVector &constituentHitVector);
    static pandora::CartesianVector GetExtremalCoordinatesLowerX(const ConstituentHitVector &constituentHitVector);
    static pandora::CartesianVector GetExtremalCoordinatesHigherX(const ConstituentHitVector &constituentHitVector);
    static float GetTotalClusterWeight(const ConstituentHitVector &constituentHitVector);
    static float GetOriginalTotalClusterWeight(const pandora::Cluster *const pCluster);
    static bool SortByHigherXExtrema(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);
    static void GetExtremalCoordinatesX(const ConstituentHitVector &constituentHitVector, pandora::CartesianVector &lowerXCoordinate, pandora::CartesianVector &higherXCoordinate);
    static void GetExtremalCoordinatesX(const pandora::CartesianPointVector &constituentHitPositionVector, pandora::CartesianVector &lowerXCoordinate, pandora::CartesianVector &higherXCoordinate);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector LArHitWidthHelper::ConstituentHit::GetPositionVector() const
{
    return m_positionVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArHitWidthHelper::ConstituentHit::GetHitWidth() const
{
    return m_hitWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArHitWidthHelper::ClusterParameters::GetNumCaloHits() const
{
    return m_numCaloHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArHitWidthHelper::ClusterParameters::GetTotalWeight() const
{
    return m_totalWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArHitWidthHelper::ConstituentHitVector LArHitWidthHelper::ClusterParameters::GetConstituentHitVector() const 
{
    return m_constituentHitVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector LArHitWidthHelper::ClusterParameters::GetLowerXExtrema() const 
{
    return m_lowerXExtrema;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector LArHitWidthHelper::ClusterParameters::GetHigherXExtrema() const 
{
    return m_higherXExtrema;
}

} // namespace lar_content


#endif // #ifndef LAR_HIT_WIDTH_HELPER_H
