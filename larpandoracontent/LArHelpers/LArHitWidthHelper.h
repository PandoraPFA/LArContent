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
    /**
     *  @brief  ConstituentHit class
     */
    class ConstituentHit
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  positionVector the central position of the constituent hit
         *  @param  hitWidth the hit width of the constituent hit
         *  @param  pParentClusterAddress the address of the original, unbroken hit to which it belongs
         */
        ConstituentHit(const pandora::CartesianVector &positionVector, const float hitWidth, const pandora::Cluster *const pParentClusterAddress);

        /**
         *  @brief  Returns the constituent hit central position
         */
        const pandora::CartesianVector &GetPositionVector() const;

        /**
         *  @brief  Returns the constituent hit width
         */
        float GetHitWidth() const;

        /**
         *  @brief  Returns the address of the parent cluster
         */
        const pandora::Cluster *GetParentClusterAddress() const;

        /**
          *  @brief  SortByDistanceToPoint class
          */
        class SortByDistanceToPoint
        {
        public:
            /**
             *  @brief  Constructor
             *
             *  @param  referencePoint the point relative to which constituent hits are ordered
             */
            SortByDistanceToPoint(const pandora::CartesianVector referencePoint) :
                m_referencePoint(referencePoint)
            {
            }

            /**
             *  @brief  Sort constituent hits by their position relative to a referencePoint
             *
             *  @param  lhs first constituent hit
             *  @param  rhs second constituent hit
             *
             *  @return  whether lhs hit is closer to the referencePoint than the rhs hit
             */
            bool operator()(const ConstituentHit &lhs, const ConstituentHit &rhs);

        private:
            const pandora::CartesianVector m_referencePoint; ///< The point relative to which constituent hits are ordered
        };

    private:
        pandora::CartesianVector m_positionVector;       ///< The central position of the consituent hit
        float m_hitWidth;                                ///< The width of the constituent hit
        const pandora::Cluster *m_pParentClusterAddress; ///< The address of the cluster the constituent hit belongs to
    };

    typedef std::vector<ConstituentHit> ConstituentHitVector;

    /**
     *  @brief  ClusterParameters class
     */
    class ClusterParameters
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pCluster from which the parameters will be obtained
         *  @param  maxConstituentHitWidth the maximum width of a constituent hit
         *  @param  isUniform whether to break up the hit into uniform constituent hits (and pad the hit) or not
         *          in the non-uniform case constituent hits from different hits may have different weights
         *  @param  hitWidthScalingFactor the constituent hit width scaling factor
         */
        ClusterParameters(const pandora::Cluster *const pCluster, const float maxConsituentHitWidth, const bool isUniformHits,
            const float hitWidthScalingFactor);

        /**
         *  @brief  Constructor
         *
         *  @param  pCluster from which the parameters will be obtained
         *  @param  numCaloHits the number of calo hits within the cluster
         *  @param  totalWeight the total weight of the constituent hits
         *  @param  constituentHitVector the vector of constituent hits
         *  @param  lowerXExtrema the lower x extremal point of the constituent hits
         *  @param  higherXExtrema the higher x extremal point of the constituent hits
         */
        ClusterParameters(const pandora::Cluster *const pCluster, const unsigned int numCaloHits, const float totalWeight,
            const ConstituentHitVector &constituentHitVector, const pandora::CartesianVector &lowerXExtrema,
            const pandora::CartesianVector &higherXExtrema);

        /**
         *  @brief  Returns the address of the cluster
         */
        const pandora::Cluster *GetClusterAddress() const;

        /**
         *  @brief  Returns the number of calo hits within the cluster
         */
        unsigned int GetNumCaloHits() const;

        /**
         *  @brief  Returns the total weight of the constituent hits
         */
        float GetTotalWeight() const;

        /**
         *  @brief  Returns the vector of constituent hits
         */
        const ConstituentHitVector &GetConstituentHitVector() const;

        /**
         *  @brief  Returns the lower x extremal point of the constituent hits
         */
        const pandora::CartesianVector &GetLowerXExtrema() const;

        /**
         *  @brief  Returns the higher x extremal point of the constituent hits
         */
        const pandora::CartesianVector &GetHigherXExtrema() const;

    private:
        const pandora::Cluster *m_pCluster;                ///< The address of the cluster
        const unsigned int m_numCaloHits;                  ///< The number of calo hits within the cluster
        const ConstituentHitVector m_constituentHitVector; ///< The vector of constituent hits
        const float m_totalWeight;                         ///< The total hit weight of the contituent hits
        const pandora::CartesianVector m_lowerXExtrema;    ///< The lower x extremal point of the constituent hits
        const pandora::CartesianVector m_higherXExtrema;   ///< The higher x extremal point of the constituent hits
    };

    typedef std::unordered_map<const pandora::Cluster *, const ClusterParameters> ClusterToParametersMap;

    /**
     *  @brief  SortByHigherExtrema class
     */
    class SortByHigherXExtrema
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param clusterToParametersMap the map [cluster -> cluster parameters]
         */
        SortByHigherXExtrema(const ClusterToParametersMap &clusterToParametersMap);

        /**
         *  @brief  Sort clusters by the higher x extremal point of their constituent hits
         *
         *  @param  pLhs first cluster
         *  @param  pRhs second cluster
         *
         *  @return  whether the pLhs cluster has a lower higherXExtrema than the pRhs cluster
         */
        bool operator()(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    private:
        const ClusterToParametersMap &m_clusterToParametersMap; ///< The map [cluster -> cluster parameters]
    };

    /**
     *  @brief  Return the cluster parameters of a given cluster, exception thrown if not found in map [cluster -> cluster parameter] or if map is empty
     *
     *  @param  pCluster the input cluster
     *  @param  clusterToParametersMap the map [cluster -> cluster parameter]
     *
     *  @return  the cluster parameters of the input cluster
     */
    static const ClusterParameters &GetClusterParameters(const pandora::Cluster *const pCluster, const ClusterToParametersMap &clusterToParametersMap);

    /**
     *  @brief  Return the number of constituent hits that a given cluster would be broken into
     *
     *  @param  pCluster the input cluster
     *  @param  maxConstituentHitWidth the maximum width of a constituent hit
     *  @param  hitWidthScalingFactor the constituent hit width scaling factor
     *
     *  @return  the number of constituent hits the cluster would be broken into
     */
    static unsigned int GetNProposedConstituentHits(
        const pandora::Cluster *const pCluster, const float maxConstituentHitWidth, const float hitWidthScalingFactor);

    /**
     *  @brief  Break up the cluster hits into constituent hits
     *
     *  @param  pCluster the input cluster
     *  @param  maxConstituentHitWidth the maximum width of a constituent hit
     *  @param  hitWidthScalingFactor the constituent hit width scaling factor
     *  @param  isUniform whether to break up the hit into uniform constituent hits (and pad the hit) or not
     *          in the non-uniform case constituent hits from different hits may have different weights
     *
     *  @return  the vector of constituent hits
     */
    static ConstituentHitVector GetConstituentHits(
        const pandora::Cluster *const pCluster, const float maxConstituentHitWidth, const float hitWidthScalingFactor, const bool isUniform);

    /**
     *  @brief  Break up the calo hit into constituent hits
     *
     *  @param  pCaloHit the input calo hit
     *  @param  pCluster the parent cluster
     *  @param  numberOfConstituentHits the number of constituent hits the hit will be broken into
     *  @param  constituentHitWidth the hit width of the constituent hits
     *  @param  constituentHitVector the input vector to which to add the contituent hits
     *
     */
    static void SplitHitIntoConstituents(const pandora::CaloHit *const pCaloHit, const pandora::Cluster *const pCluster,
        const unsigned int numberOfConstituentHits, const float constituentHitWidth, ConstituentHitVector &constituentHitVector);

    /**
     *  @brief  Obtain a vector of the contituent hit central positions
     *
     *  @param  constituentHitVector the input vector of contituent hits
     *
     *  @return  a vector of constituent hit central positions
     */
    static pandora::CartesianPointVector GetConstituentHitPositionVector(const ConstituentHitVector &constituentHitVector);

    /**
     *  @brief  Sum the widths of constituent hits
     *
     *  @param  constituentHitVector the input vector of contituent hits
     *
     *  @return  the total weight sum
     */
    static float GetTotalClusterWeight(const ConstituentHitVector &constituentHitVector);

    /**
     *  @brief  Sum the widths of the original, unscaled hits contained within a cluster
     *
     *  @param  pCluster the input cluster
     *
     *  @return  the total weight sum
     */
    static float GetOriginalTotalClusterWeight(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Return the lower x extremal point of the constituent hits
     *
     *  @param  constituentHitVector the input vector of contituent hits
     *
     *  @return  the lower x extremal point of the constituent hits
     */
    static pandora::CartesianVector GetExtremalCoordinatesLowerX(const ConstituentHitVector &constituentHitVector);

    /**
     *  @brief  Return the higher x extremal point of the constituent hits
     *
     *  @param  constituentHitVector the input vector of contituent hits
     *
     *  @return  the higher x extremal point of the constituent hits
     */
    static pandora::CartesianVector GetExtremalCoordinatesHigherX(const ConstituentHitVector &constituentHitVector);

    /**
     *  @brief  Calculate the higher and lower x extremal points of the constituent hits
     *
     *  @param  constituentHitVector the input vector of contituent hits
     *  @param  lowerXCoordinate the lower x extremal point
     *  @param  higherXCoordinate the higher x extremal point
     */
    static void GetExtremalCoordinatesX(const ConstituentHitVector &constituentHitVector, pandora::CartesianVector &lowerXCoordinate,
        pandora::CartesianVector &higherXCoordinate);

    /**
     *  @brief  Consider the hit width to find the closest position of a calo hit to a specified line
     *
     *  @param  lineStart the start position of the line
     *  @param  lineDirection the direction of the line
     *  @param  pCaloHit the input calo hit
     *
     *  @return  the closest position
     */
    static pandora::CartesianVector GetClosestPointToLine2D(
        const pandora::CartesianVector &lineStart, const pandora::CartesianVector &lineDirection, const pandora::CaloHit *const pCaloHit);

    /**
     *  @brief  Consider the hit width to find the smallest distance between a calo hit and a given point
     *
     *  @param  pCaloHit the input calo hit
     *  @param  point2D the position
     *
     *  @return  the smallest distance
     */
    static float GetClosestDistanceToPoint2D(const pandora::CaloHit *const pCaloHit, const pandora::CartesianVector &point2D);

    /**
     *  @brief  Find the smallest separation between a hit and a list of hits, with the consideration of their hit widths
     *
     *  @param  pThisCaloHit the input calo hit
     *  @param  caloHitList the input calo hit list
     *
     *  @return the smallest separation
     */
    static float GetClosestDistance(const pandora::CaloHit *const pThisCaloHit, const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Find the smallest separation between two hits, with the consideration of their hit widths
     *
     *  @param  pCaloHit1 the first calo hit
     *  @param  pCaloHit2 the second calo hit
     *
     *  @return the smallest separation
     */
    static float GetClosestDistance(const pandora::CaloHit *const pCaloHit1, const pandora::CaloHit *const pCaloHit2);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArHitWidthHelper::ConstituentHit::GetPositionVector() const
{
    return m_positionVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArHitWidthHelper::ConstituentHit::GetHitWidth() const
{
    return m_hitWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *LArHitWidthHelper::ConstituentHit::GetParentClusterAddress() const
{
    return m_pParentClusterAddress;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *LArHitWidthHelper::ClusterParameters::GetClusterAddress() const
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

inline const LArHitWidthHelper::ConstituentHitVector &LArHitWidthHelper::ClusterParameters::GetConstituentHitVector() const
{
    return m_constituentHitVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArHitWidthHelper::ClusterParameters::GetLowerXExtrema() const
{
    return m_lowerXExtrema;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArHitWidthHelper::ClusterParameters::GetHigherXExtrema() const
{
    return m_higherXExtrema;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArHitWidthHelper::SortByHigherXExtrema::SortByHigherXExtrema(const ClusterToParametersMap &clusterToParametersMap) :
    m_clusterToParametersMap(clusterToParametersMap)
{
}

} // namespace lar_content

#endif // #ifndef LAR_HIT_WIDTH_HELPER_H
