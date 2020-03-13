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

    private:
        pandora::CartesianVector m_positionVector;
        float m_hitWidth;
        const pandora::Cluster *m_parentClusterAddress; 
    };

    typedef std::vector<ConstituentHit> ConstituentHitVector;


   class ClusterParameters
    {
    public:
      ClusterParameters(const pandora::Cluster *const pCluster, const float maxConsituentHitWidth, bool isUniformHits);

        const pandora::Cluster* GetClusterAddress() const;
	unsigned int GetNumCaloHits() const;
	float GetTotalWeight() const;
	ConstituentHitVector GetConstituentHitVector() const;
	pandora::CartesianVector GetLowerXExtrema() const;
	pandora::CartesianVector GetHigherXExtrema() const;

    private:
        const pandora::Cluster *m_pCluster;                 
        unsigned int m_numCaloHits;
        float m_totalWeight;
        ConstituentHitVector m_constituentHitVector;
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
    

    static ConstituentHitVector GetConstituentHits(const pandora::Cluster *const pCluster, const float maxConstituentHitWidth);
    static ConstituentHitVector GetUniformConstituentHits(const pandora::Cluster *const pCluster, const float constituentHitWidth);
    static pandora::CartesianPointVector GetConstituentHitPositionVector(const ConstituentHitVector &constituentHitVector);
    static pandora::CartesianVector GetExtremalCoordinatesLowerX(const ConstituentHitVector &constituentHitVector);
    static pandora::CartesianVector GetExtremalCoordinatesHigherX(const ConstituentHitVector &constituentHitVector);
    static float GetTotalClusterWeight(const ConstituentHitVector &constituentHitVector);
    static float GetTotalClusterWeight(const pandora::Cluster *const pCluster);
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
   
inline const pandora::Cluster* LArHitWidthHelper::ConstituentHit::GetParentClusterAddress() const
{
    return m_parentClusterAddress;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster* LArHitWidthHelper::ClusterParameters::GetClusterAddress() const
{
    return m_pCluster;
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
