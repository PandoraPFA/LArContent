/**
 *  @file   LArContent/include/LArTwoDSeed/SeedConsolidationAlgorithm.h
 * 
 *  @brief  Header file for the seed consolidation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_CONSOLIDATION_ALGORITHM_H
#define LAR_SEED_CONSOLIDATION_ALGORITHM_H 1

#include "Objects/OrderedCaloHitList.h"

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  SeedConsolidationAlgorithm class
 */
class SeedConsolidationAlgorithm : public pandora::Algorithm
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
     *  @brief  ConeAssociation class
     */
    class ConeAssociation
    {
    public:
        /**
         *  @brief  Constructor (note asymmetric response to choice of parent and daughter clusters)
         * 
         *  @param  pParentCluster address of parent cluster
         *  @param  pDaughterCluster address of daughter cluster
         */
        ConeAssociation(const pandora::Cluster *const pParentCluster, const pandora::Cluster *const pDaughterCluster);

        /**
         *  @brief  Get the cosine of the relative angle between parent and daughter clusters
         * 
         *  @return the cosine of the relative angle between parent and daughter clusters
         */
        float GetCosRelativeAngle() const;

        /**
         *  @brief  Get the daughter start distance
         * 
         *  @return the daughter start distance
         */
        float GetDaughterStartDistance() const;

        /**
         *  @brief  Get the daughter end distance
         * 
         *  @return the daughter end distance
         */
        float GetDaughterEndDistance() const;

        /**
         *  @brief  Get the apex to daughter distance
         * 
         *  @return the apex to daughter distance
         */
        float GetApexToDaughterDistance() const;

        /**
         *  @brief  Get the apex to daughter distance as a fraction of the cone length
         * 
         *  @return the apex to daughter distance fraction
         */
        float GetApexToDaughterDistanceFraction() const;

        /**
         *  @brief  Get the parent to daughter distance
         * 
         *  @return the parent to daughter distance
         */
        float GetParentToDaughterDistance() const;

        /**
         *  @brief  Get the cone apex
         * 
         *  @return the cone apex
         */
        const pandora::CartesianVector &GetConeApex() const;

        /**
         *  @brief  Get the cone direction
         * 
         *  @return the cone direction
         */
        const pandora::CartesianVector &GetConeDirection() const;

        /**
         *  @brief  Whether the cone direction is forwards in z
         * 
         *  @return boolean
         */
        bool IsForwardsInZ() const;

        /**
         *  @brief  Get cone cos half angle of parent cluster
         * 
         *  @return the cone cos half angle of parent cluster
         */
        float GetConeCosHalfAngleParent() const;

        /**
         *  @brief  Get cone cos half angle of daughter cluster
         * 
         *  @return the cone cos half angle of daughter cluster
         */
        float GetConeCosHalfAngleDaughter() const;

        /**
         *  @brief  Get the enclosed fraction of daughter hits bounded by the parent
         * 
         *  @return the enclosed hit fraction
         */
        float GetEnclosedDaughterHitFraction() const; 

        /**
         *  @brief  Get the enclosed fraction of parent hits bounded by the daughter
         * 
         *  @return the enclosed hit fraction
         */
        float GetEnclosedParentHitFraction() const;

        /**
         *  @brief  Get the cone length
         * 
         *  @return the cone length
         */
        float GetConeLength() const;

    private:
        /**
         *  @brief  Calculate association properties (note asymmetric response to choice of parent and daughter clusters)
         * 
         *  @param  pParentCluster address of parent cluster
         *  @param  pDaughterCluster address of daughter cluster
         */
        void CalculateProperties(const pandora::Cluster *const pParentCluster, const pandora::Cluster *const pDaughterCluster);

        float                       m_cosRelativeAngle;             ///< 
        float                       m_daughterStartDistance;        ///< 
        float                       m_daughterEndDistance;          ///< 
        float                       m_apexToDaughterDistance;       ///< 
        float                       m_parentToDaughterDistance;     ///< 

        pandora::CartesianVector    m_coneApex;                     ///< 
        pandora::CartesianVector    m_coneDirection;                ///< 
        float                       m_coneLength;                   ///<
        bool                        m_isForwardsInZ;                ///< 
        float                       m_coneCosHalfAngleParent;       ///< 
        float                       m_coneCosHalfAngleDaughter;     ///< 
        float                       m_enclosedParentHitFraction;    ///< 
        float                       m_enclosedDaughterHitFraction;  ///< 
    };

    /**
     *  @brief  OverlapAssociation class
     */
    class OverlapAssociation
    {
    public:
        /**
         *  @brief  Constructor (note asymmetric response to choice of parent and daughter clusters)
         * 
         *  @param  pParentCluster address of parent cluster
         *  @param  pDaughterCluster address of daughter cluster
         */
        OverlapAssociation(const pandora::Cluster *const pParentCluster, const pandora::Cluster *const pDaughterCluster);

        /**
         *  @brief  Get inner overlap distance
         * 
         *  @return the inner overlap distance
         */
        float GetInnerOverlapDistance() const;

        /**
         *  @brief  Get outer overlap distance
         * 
         *  @return the outer overlap distance
         */
        float GetOuterOverlapDistance() const;

        /**
         *  @brief  Get inner overlap layer
         * 
         *  @return the inner overlap layer
         */
        unsigned int GetInnerOverlapLayer() const;

        /**
         *  @brief  Get outer overlap layer
         * 
         *  @return the outer overlap layer
         */
        unsigned int GetOuterOverlapLayer() const;

        /**
         *  @brief  Get number of overlap layers
         * 
         *  @return the number of overlap layers
         */
        unsigned int GetNOverlapLayers() const;

        /**
         *  @brief  Get number of contact layers
         * 
         *  @return the number of contact layers
         */
        unsigned int GetNContactLayers() const;

        /**
         *  @brief  Get number of non-contact layers
         * 
         *  @return the number of non-contact layers
         */
        unsigned int GetNNonContactLayers() const;

        /**
         *  @brief  Get number of contact groups
         * 
         *  @return the number of contact groups
         */
        unsigned int GetNContactGroups() const;

        /**
         *  @brief  Get number of non-vertex contact groups
         * 
         *  @return the number of non-vertex contact groups
         */
        unsigned int GetNNonVtxContactGroups() const;

        /**
         *  @brief  Get number of enclosed layers
         * 
         *  @return the number of enclosed layers
         */
        unsigned int GetNEnclosedLayers() const;

        /**
         *  @brief  Get layer distance 50
         * 
         *  @return the layer distance 50
         */
        float GetLayerDistance50() const;

    private:
        /**
         *  @brief  Calculate association properties (note asymmetric response to choice of parent and daughter clusters)
         * 
         *  @param  pParentCluster address of parent cluster
         *  @param  pDaughterCluster address of daughter cluster
         */
        void CalculateProperties(const pandora::Cluster *const pParentCluster, const pandora::Cluster *const pDaughterCluster);

        /**
         *  @brief  Whether hits in specified daughter layer are in contact with those in parent cluster
         * 
         *  @param  daughterListIter iterator specifying daughter layer and hits in this layer
         *  @param  parentList the parent ordered calo hit list
         * 
         *  @return boolean
         */
        void EvaluateLayerContact(const pandora::OrderedCaloHitList::const_iterator daughterListIter, const pandora::OrderedCaloHitList &parentList,
             float &closestHitDistance, bool &enclosedHits) const;

        float                       m_innerOverlapDistance;     ///< 
        float                       m_outerOverlapDistance;     ///< 
        unsigned int                m_innerOverlapLayer;        ///< 
        unsigned int                m_outerOverlapLayer;        ///< 
        unsigned int                m_nOverlapLayers;           ///< 

        unsigned int                m_nContactLayers;          ///< 
        unsigned int                m_nNonContactLayers;       ///< 
        unsigned int                m_nContactGroups;          ///< 
        unsigned int                m_nNonVtxContactGroups;    ///< 

        unsigned int                m_nEnclosedLayers;          ///< 
        float                       m_layerDistance50;          ///< 
    };

    /**
     *  @brief  Whether the cone association details are sufficient to make cluster merge
     * 
     *  @param  coneAssociation the cone association details
     * 
     *  @return boolean
     */
    bool PassesConeAssociation(const ConeAssociation &coneAssociation) const;

    /**
     *  @brief  Whether the clusters pass recovery cuts (similar, but more generous, than those in particle seed algorithm)
     * 
     *  @param  coneAssociation the cone association details
     * 
     *  @return boolean
     */
    bool PassesRecoveryCuts(const ConeAssociation &coneAssociation) const;

    /**
     *  @brief  Whether the overlap association details are sufficient to make cluster merge
     * 
     *  @param  overlapAssociation the overlap association details
     * 
     *  @return boolean
     */
    bool PassesOverlapAssociation(const OverlapAssociation &overlapAssociation) const;

    /**
     *  @brief  Whether a pair of clusters are both associated with the vertex
     * 
     *  @param  pClusterParent address of parent cluster
     *  @param  pClusterDaughter address of daughter cluster
     * 
     *  @return boolean
     */
    bool PassesVertexAssociation(const pandora::Cluster* const pClusterParent, const pandora::Cluster* const pClusterDaughter);

    typedef std::vector<pandora::ClusterVector::iterator> ClusterIteratorList;
    typedef std::vector<pandora::CartesianVector> CartesianVectorList;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SeedConsolidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new SeedConsolidationAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetCosRelativeAngle() const
{
    return m_cosRelativeAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetDaughterStartDistance() const
{
    return m_daughterStartDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetDaughterEndDistance() const
{
    return m_daughterEndDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetApexToDaughterDistance() const
{
    return m_apexToDaughterDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetApexToDaughterDistanceFraction() const
{
  if ( m_coneLength>0.f ) return m_apexToDaughterDistance / m_coneLength;  else return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetParentToDaughterDistance() const
{
    return m_parentToDaughterDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &SeedConsolidationAlgorithm::ConeAssociation::GetConeApex() const
{
    return m_coneApex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &SeedConsolidationAlgorithm::ConeAssociation::GetConeDirection() const
{
    return m_coneDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool SeedConsolidationAlgorithm::ConeAssociation::IsForwardsInZ() const
{
    return m_isForwardsInZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetConeCosHalfAngleParent() const
{
    return m_coneCosHalfAngleParent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetConeCosHalfAngleDaughter() const
{
    return m_coneCosHalfAngleDaughter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetEnclosedParentHitFraction() const
{
    return m_enclosedParentHitFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetEnclosedDaughterHitFraction() const
{
    return m_enclosedDaughterHitFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::ConeAssociation::GetConeLength() const
{
    return m_coneLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::OverlapAssociation::GetInnerOverlapDistance() const
{
    return m_innerOverlapDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------


inline float SeedConsolidationAlgorithm::OverlapAssociation::GetOuterOverlapDistance() const
{
    return m_outerOverlapDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedConsolidationAlgorithm::OverlapAssociation::GetInnerOverlapLayer() const
{
    return m_innerOverlapLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedConsolidationAlgorithm::OverlapAssociation::GetOuterOverlapLayer() const
{
    return m_outerOverlapLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedConsolidationAlgorithm::OverlapAssociation::GetNOverlapLayers() const
{
    return m_nOverlapLayers;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedConsolidationAlgorithm::OverlapAssociation::GetNContactLayers() const
{
    return m_nContactLayers;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedConsolidationAlgorithm::OverlapAssociation::GetNNonContactLayers() const
{
    return m_nNonContactLayers;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedConsolidationAlgorithm::OverlapAssociation::GetNContactGroups() const
{
    return m_nContactGroups;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedConsolidationAlgorithm::OverlapAssociation::GetNNonVtxContactGroups() const
{
    return m_nNonVtxContactGroups;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedConsolidationAlgorithm::OverlapAssociation::GetNEnclosedLayers() const
{
    return m_nEnclosedLayers;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float SeedConsolidationAlgorithm::OverlapAssociation::GetLayerDistance50() const
{
    return m_layerDistance50;
}

} // namespace lar

#endif // #ifndef LAR_SEED_CONSOLIDATION_ALGORITHM_H
