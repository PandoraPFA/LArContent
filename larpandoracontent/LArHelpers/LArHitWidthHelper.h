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

        pandora::CartesianVector m_positionVector;
        float m_hitWidth;
        const pandora::Cluster *m_parentClusterAddress; 
    };

    typedef std::vector<ConstituentHit> ConstituentHitVector;


   class ClusterParameters
    {
    public:
        ClusterParameters(const pandora::Cluster *const pCluster, const float maxConsituentHitWidth);

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




};








    static ConstituentHitVector GetConstituentHits(const pandora::Cluster *const pCluster, const float maxConstituentHitWidth);
    static ConstituentHitVector GetUniformConstituentHits(const pandora::Cluster *const pCluster, const float constituentHitWidth);
    static pandora::CartesianVector GetExtremalCoordinatesLowerX(const ConstituentHitVector &constituentHitVector);
    static pandora::CartesianVector GetExtremalCoordinatesHigherX(const ConstituentHitVector &constituentHitVector);
    static float GetTotalClusterWeight(const ConstituentHitVector &constituentHitVector);
    static float GetTotalClusterWeight(const pandora::Cluster *const pCluster);
    static void GetExtremalCoordinatesX(const ConstituentHitVector &constituentHitVector, pandora::CartesianVector &lowerXCoordinate, pandora::CartesianVector &higherXCoordinate);
    static void GetExtremalCoordinatesX(const pandora::CartesianPointVector &constituentHitPositionVector, pandora::CartesianVector &lowerXCoordinate, pandora::CartesianVector &higherXCoordinate);
};

} // namespace lar_content


#endif // #ifndef LAR_HIT_WIDTH_HELPER_H
