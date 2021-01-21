/**
 *  @file   larpandoracontent/LArHelpers/LArMuonLeadingHelper.h
 *
 *  @brief  Header file for the lar muon leading helper class.
 *
 *  $Log: $
 */
#ifndef LAR_MUON_LEADING_HELPER_H
#define LAR_MUON_LEADING_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  LArMuonLeadingHelper class
 */
class LArMuonLeadingHelper
{
public:
    /**
     *  @brief   MuonLeadingParameters class
     */
    class ValidationParameters : public LArMCParticleHelper::PrimaryParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        ValidationParameters();

        float m_maxBremsstrahlungSeparation;
    };

    typedef std::map<const pandora::MCParticle*, std::map<const pandora::MCParticle*, pandora::CaloHitList>> LeadingMCParticleToPostPhotonHitLists;
    
    /**
     *  @brief  Return true if passed a DR tagged MCParticle 
     */    
    static bool IsDeltaRay(const pandora::MCParticle *const pMCParticle);

    static bool IsMichel(const pandora::MCParticle *const pMCParticle);

    static bool IsLeading(const pandora::MCParticle *const pMCParticle);   
    
    static const pandora::MCParticle *GetLeadingParticle(const pandora::MCParticle *const pMCParticle);

    static void GetMCToLeadingMap(const pandora::MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToLeadingMap);

    static void SelectReconstructableLeadingParticles(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, const ValidationParameters &parameters,
       const pandora::CaloHitList &recoMuonHitList, LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap,  const pandora::Pandora &pandora);

    static void SelectCaloHits(const pandora::CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap,
        pandora::CaloHitList &selectedCaloHitList, const bool selectInputHits, const float minHitSharingFraction,
        const pandora::CaloHitList &recoMuonHitList, LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists); 

    static void GetPfoMatchContamination(const pandora::MCParticle *const pLeadingParticle, const pandora::CaloHitList &matchedPfoHitList,
        pandora::CaloHitList &parentTrackHits, pandora::CaloHitList &otherTrackHits, pandora::CaloHitList &otherShowerHits);

    static void GetMuonPfoContaminationContribution(const pandora::CaloHitList &cosmicRayPfoHitList, const pandora::CaloHitList &leadingMCHitList,
        pandora::CaloHitList &leadingHitsInParentCosmicRay);

    static bool RejectBremsstrahlungHits(const pandora::CaloHit *const pCaloHit, LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists);

    static void AddInReconstructablePostPhotonHits(const LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists, const float maxBremsstrahlungSeparation,
                                                   LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const pandora::Pandora &pandora);
    static void AddHits(const pandora::MCParticle *const pLeadingParticle, const LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists, const float maxBremsstrahlungSeparation,
                        LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const pandora::HitType &tpcView, const pandora::Pandora &pandora);

    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &referencePoint, const pandora::CartesianPointVector &cartesianPointVector,
        const pandora::Cluster *const pCluster);

    static float GetClosestDistance(const pandora::Cluster *const pCluster, const pandora::CartesianPointVector &cartesianPointVector);

    static float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CartesianPointVector &cartesianPointVector);

    static float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &caloHitList);

    static void GetClosestPositions(const pandora::CartesianPointVector &pCluster1, const pandora::Cluster *const pCluster2, pandora::CartesianVector &outputPosition1,
        pandora::CartesianVector &outputPosition2);

 private:
    static void SelectLeadingMCParticles(const pandora::MCParticleList *pMCParticleList, pandora::MCParticleVector &selectedParticles);

};

} // namespace lar_content

#endif // #ifndef LAR_MUON_LEADING_HELPER_H
