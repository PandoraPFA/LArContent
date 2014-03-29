/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterCreation/ClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_CREATION_ALGORITHM_H
#define LAR_CLUSTER_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ClusterCreationAlgorithm class
 */
class ClusterCreationAlgorithm : public pandora::Algorithm
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
    /**
     *  @brief  HitAssociation class
     */
    class HitAssociation
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pPrimaryTarget address of the primary target hit
         *  @param  primaryDistanceSquared distance to the primary target hit squared
         */
        HitAssociation(pandora::CaloHit *const pPrimaryTarget, const float primaryDistanceSquared);

        /**
         *  @brief  Set secondary target
         *
         *  @param  pSecondaryTarget address of the secondary target hit
         *  @param  secondaryDistanceSquared distance to the primary target hit squared
         */
        void SetSecondaryTarget(pandora::CaloHit *const pSecondaryTarget, const float secondaryDistanceSquared);

        /**
         *  @brief  Get the primary target
         *
         *  @return the target distance
         */
        pandora::CaloHit *GetPrimaryTarget() const;

        /**
         *  @brief  Get the secondary target
         *
         *  @return the secondary target
         */
        pandora::CaloHit *GetSecondaryTarget() const;

        /**
         *  @brief  Get the primary distance squared
         *
         *  @return the primary distance squared
         */
        float GetPrimaryDistanceSquared() const;

        /**
         *  @brief  Get the secondary distance squared
         *
         *  @return the secondary distance squared
         */
        float GetSecondaryDistanceSquared() const;

        static float        m_maxSeparationSquared;         ///< Maximum separation squared
        static float        m_closeSeparationSquared;       ///< Close separation squared

    private:
        pandora::CaloHit   *m_pPrimaryTarget;               ///< the primary target
        pandora::CaloHit   *m_pSecondaryTarget;             ///< the secondary target
        float               m_primaryDistanceSquared;       ///< the primary distance squared
        float               m_secondaryDistanceSquared;     ///< the secondary distance squared
    };

    typedef std::map<pandora::CaloHit*, HitAssociation> HitAssociationMap;
    typedef std::map<pandora::CaloHit*, pandora::CaloHit*> HitJoinMap;
    typedef std::map<pandora::CaloHit*, pandora::Cluster*> HitToClusterMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Filter out low pulse height hits in close proximity to high pulse height hits
     *
     *  @param  pCaloHitList input hit list
     *  @param  selectedCaloHitList the output selected list of selected hits
     *  @param  rejectedCaloHitList the output rejected list of rejected hits
     */
    pandora::StatusCode FilterCaloHits(const  pandora::CaloHitList *pCaloHitList, pandora::OrderedCaloHitList &selectedCaloHitList, pandora::OrderedCaloHitList& rejectedCaloHitList) const;

    /**
     *  @brief  Merge previously filtered hits back into their associated clusters
     *
     *  @param  selectedCaloHitList the ordered list of selected hits
     *  @param  rejectedCaloHitList the ordered list of rejected hits
     *  @param  hitToClusterMap the mapping between hits and their clusters
     */
    pandora::StatusCode AddFilteredCaloHits(const pandora::OrderedCaloHitList &selectedCaloHitList, const pandora::OrderedCaloHitList& rejectedCaloHitList, HitToClusterMap& hitToClusterMap) const;

    /**
     *  @brief  Control primary association formation
     *
     *  @param  orderedCaloHitList the ordered calo hit list
     *  @param  forwardHitAssociationMap the forward hit association map
     *  @param  backwardHitAssociationMap the backward hit association map
     */
    void MakePrimaryAssociations(const pandora::OrderedCaloHitList &orderedCaloHitList, HitAssociationMap &forwardHitAssociationMap,
        HitAssociationMap &backwardHitAssociationMap) const;

    /**
     *  @brief  Control secondary association formation
     *
     *  @param  orderedCaloHitList the ordered calo hit list
     *  @param  forwardHitAssociationMap the forward hit association map
     *  @param  backwardHitAssociationMap the backward hit association map
     */
    void MakeSecondaryAssociations(const pandora::OrderedCaloHitList &orderedCaloHitList, HitAssociationMap &forwardHitAssociationMap,
        HitAssociationMap &backwardHitAssociationMap) const;

    /**
     *  @brief  Identify final hit joins for use in cluster formation
     *
     *  @param  orderedCaloHitList the ordered calo hit list
     *  @param  forwardHitAssociationMap the forward hit association map
     *  @param  backwardHitAssociationMap the backward hit association map
     *  @param  hitJoinMap to receive the hit join map
     */
    void IdentifyJoins(const pandora::OrderedCaloHitList &orderedCaloHitList, const HitAssociationMap &forwardHitAssociationMap,
        const HitAssociationMap &backwardHitAssociationMap, HitJoinMap &hitJoinMap) const;

    /**
     *  @brief  Final cluster formation
     *
     *  @param  orderedCaloHitList the ordered calo hit list
     *  @param  hitJoinMap the hit join map
     *  @param  hitToClusterMap the mapping between hits and their clusters
     */
    void CreateClusters(const pandora::OrderedCaloHitList &orderedCaloHitList, const HitJoinMap &hitJoinMap, HitToClusterMap& hitToClusterMap) const;

    /**
     *  @brief  Create primary association if appropriate, hitI<->hitJ
     *
     *  @param  pCaloHitI address of calo hit I
     *  @param  pCaloHitJ address of calo hit J
     *  @param  forwardHitAssociationMap the forward hit association map
     *  @param  backwardHitAssociationMap the backward hit association map
     */
    void CreatePrimaryAssociation(pandora::CaloHit *pCaloHitI, pandora::CaloHit *pCaloHitJ, HitAssociationMap &forwardHitAssociationMap,
        HitAssociationMap &backwardHitAssociationMap) const;

    /**
     *  @brief  Create secondary association if appropriate, hitI<->hitJ
     *
     *  @param  pCaloHitI address of calo hit I
     *  @param  pCaloHitJ address of calo hit J
     *  @param  forwardHitAssociationMap the forward hit association map
     *  @param  backwardHitAssociationMap the backward hit association map
     */
    void CreateSecondaryAssociation(pandora::CaloHit *pCaloHitI, pandora::CaloHit *pCaloHitJ, HitAssociationMap &forwardHitAssociationMap,
        HitAssociationMap &backwardHitAssociationMap) const;

    /**
     *  @brief  Get hit to join by tracing associations via map I, checking via map J
     *
     *  @param  pCaloHit the initial calo hit
     *  @param  hitAssociationMapI hit association map I
     *  @param  hitAssociationMapJ hit association map J
     *
     *  @return the hit to join
     */
    pandora::CaloHit *GetJoinHit(pandora::CaloHit *pCaloHit, const HitAssociationMap &hitAssociationMapI, const HitAssociationMap &hitAssociationMapJ) const;

    /**
     *  @brief  Get last hit obtained by tracing associations via map I, checking via map J
     *
     *  @param  pCaloHit the initial calo hit
     *  @param  hitAssociationMapI hit association map I
     *  @param  hitAssociationMapJ hit association map J
     *  @param  nSteps to receive the number of association steps
     *
     *  @return the last hit obtained in the chain of associations
     */
    pandora::CaloHit *TraceHitAssociation(pandora::CaloHit *pCaloHit, const HitAssociationMap &hitAssociationMapI, const HitAssociationMap &hitAssociationMapJ,
        unsigned int &nSteps) const;

    std::string          m_inputCaloHitListName;         ///< The input calo hit list name
    std::string          m_outputClusterListName;        ///< The output cluster list name

    bool                 m_mergeBackFilteredHits;        ///< Merge rejected hits into their associated clusters
    unsigned int         m_maxGapLayers;                 ///< Maximum number of layers for a gap
    float                m_minCaloHitSeparationSquared;  ///< Square of minimum calo hit separation
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterCreationAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ClusterCreationAlgorithm::HitAssociation::HitAssociation(pandora::CaloHit *const pPrimaryTarget, const float primaryDistanceSquared) :
    m_pPrimaryTarget(pPrimaryTarget),
    m_pSecondaryTarget(NULL),
    m_primaryDistanceSquared(primaryDistanceSquared),
    m_secondaryDistanceSquared(m_closeSeparationSquared)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ClusterCreationAlgorithm::HitAssociation::SetSecondaryTarget(pandora::CaloHit *const pSecondaryTarget, const float secondaryDistanceSquared)
{
    m_pSecondaryTarget = pSecondaryTarget;
    m_secondaryDistanceSquared = secondaryDistanceSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CaloHit *ClusterCreationAlgorithm::HitAssociation::GetPrimaryTarget() const
{
    return m_pPrimaryTarget;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CaloHit *ClusterCreationAlgorithm::HitAssociation::GetSecondaryTarget() const
{
    return m_pSecondaryTarget;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ClusterCreationAlgorithm::HitAssociation::GetPrimaryDistanceSquared() const
{
    return m_primaryDistanceSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ClusterCreationAlgorithm::HitAssociation::GetSecondaryDistanceSquared() const
{
    return m_secondaryDistanceSquared;
}

} // namespace lar

#endif // #ifndef LAR_CLUSTER_CREATION_ALGORITHM_H
