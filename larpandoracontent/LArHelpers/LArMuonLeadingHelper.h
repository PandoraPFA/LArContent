/**
 *  @file   larpandoracontent/LArHelpers/LArMuonLeadingHelper.h
 *
 *  @brief  Header file for the muon leading helper class.
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
     *  @brief   ValidationParameters class
     */
    class ValidationParameters : public LArMCParticleHelper::PrimaryParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        ValidationParameters();

        float m_maxBremsstrahlungSeparation; ///< The maximum separation of a reconstructable post-bremsstrahlung hit from the pre-radiation hits
    };

    typedef std::map<const pandora::MCParticle*, pandora::CaloHitList> LeadingMCParticleToPostBremsstrahlungHitList;
    
    /**
     *  @brief  Return true if input MCParticle is a child of a cosmic ray and has 'delta ray' process tag
     *
     *  @param  pMCParticle the input MCParticle
     */
    static bool IsDeltaRay(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Return true if input MCParticle is a child of a cosmic ray and has 'decay' process tag
     *
     *  @param  pMCParticle the input MCParticle
     */
    static bool IsMichel(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Return true if input MCParticle is in tier 1 of the cosmic ray hierarchy
     *
     *  @param  pMCParticle the input MCParticle
     */
    static bool IsMuonLeading(const pandora::MCParticle *const pMCParticle);   

    /**
     *  @brief  Return leading particle in the delta ray/michel hierarchy i.e. tier below cosmic ray
     *
     *  @param  pMCParticle the input MCParticle
     */
    static const pandora::MCParticle *GetLeadingParticle(const pandora::MCParticle *const pMCParticle);

    /**
     *  @brief  Select target, reconstructable mc particles in the cosmic ray hierarchy
     *
     *  @param  pMCParticleList the address of the list of MCParticles
     *  @param  pCaloHitList the address of the list of CaloHits
     *  @param  parameters validation parameters to decide when an MCParticle is considered reconstructable
     *  @param  recoMuonHitList the list of reconstructed cosmic ray hits 
     *  @param  selectedMCParticlesToHitsMap the output mapping from selected MCParticles to their hits
     */
    static void SelectReconstructableLeadingParticles(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, 
       const ValidationParameters &parameters, const pandora::CaloHitList &recoMuonHitList, LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap);

    /**
     *  @brief  Separate a leading pfo hit list according to the true owner of the hit e.g. other shower
     *
     *  @param  pLeadingParticle the address of the input leading MCParticle
     *  @param  matchedPfoHitList the input leading pfo hit list
     *  @param  parentTrackHits the output list of hits that belong to the parent cosmic ray 
     *  @param  otherTrackHits the output list of hits that belong to a cosmic ray that is not the parent
     *  @param  otherShowerHits the output list of hits that belong to a different shower hierarchy
     */
    static void GetPfoMatchContamination(const pandora::MCParticle *const pLeadingParticle, const pandora::CaloHitList &matchedPfoHitList,
        pandora::CaloHitList &parentTrackHits, pandora::CaloHitList &otherTrackHits, pandora::CaloHitList &otherShowerHits);

    /**
     *  @brief  Determine the leading MCParticle hits within a cosmic ray pfo hit list
     *
     *  @param  cosmicRayPfoHitList the cosmic ray pfo hit list
     *  @param  leadingMCHitList the list of hits that belong to the leading MCParticle
     *  @param  leadingHitsInParentCosmicRay the output list of 'stolen' leading MCParticle hits
     */
    static void GetMuonPfoContaminationContribution(const pandora::CaloHitList &cosmicRayPfoHitList, const pandora::CaloHitList &leadingMCHitList,
        pandora::CaloHitList &leadingHitsInParentCosmicRay);

    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &referencePoint, const pandora::CartesianPointVector &cartesianPointVector,
        const pandora::Cluster *const pCluster);

    static float GetClosestDistance(const pandora::Cluster *const pCluster, const pandora::CartesianPointVector &cartesianPointVector);

    static float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CartesianPointVector &cartesianPointVector);

    static float GetClosestDistance(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &caloHitList);

    static void GetClosestPositions(const pandora::CartesianPointVector &pCluster1, const pandora::Cluster *const pCluster2, pandora::CartesianVector &outputPosition1,
        pandora::CartesianVector &outputPosition2);

 private:
    /**
     *  @brief  Construct the hierarchy folding map (cosmic rays folded to themselves, delta ray/michel hierarchy folded to leading particle)
     *
     *  @param  pMCParticleList the address of the list of MCParticles
     *  @param  mcToLeadingMap the hierarchy folding map
     */
    static void GetMCToLeadingMap(const pandora::MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToLeadingMap);

    /**
     *  @brief  Select a subset of calo hits representing those that represent "reconstructable" regions of the event
     *
     *  @param  pCaloHitList the address of the input calo hit list
     *  @param  mcToTargetMCMap the MCParticle to leading MCParticle map
     *  @param  selectedCaloHitList to receive the populated selected calo hit list
     *  @param  selectInputHits whether to select input hits
     *  @param  minHitSharingFraction the minimum required charge share of the hit
     *  @param  recoMuonHitList the list of reconstructed cosmic ray hits 
     *  @param  leadingMCParticleToPostBremsstrahlungHitList the mapping of leading MCParticles to post-bremsstrahlung hits
     */
    static void SelectCaloHits(const pandora::CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap,
        pandora::CaloHitList &selectedCaloHitList, const bool selectInputHits, const float minHitSharingFraction,
        const pandora::CaloHitList &recoMuonHitList, LeadingMCParticleToPostBremsstrahlungHitList &leadingMCParticleToPostBremsstrahlungHitList); 

    /**
     *  @brief  Identify and record the hits that are post-bremsstralung radiation in a cosmic ray hierarchy
     *
     *  @param  pCaloHit the address of the input calo hit
     *  @param  leadingMCParticleToPostBremsstrahlungHitList the mapping of leading MCParticles to post-bremsstrahlung hits
     *
     *  @return whether the hit lies post-bremsstrahlung radiation in a cosmic ray hierarchy
     */
    static bool RejectBremsstrahlungHits(const pandora::CaloHit *const pCaloHit, LeadingMCParticleToPostBremsstrahlungHitList &leadingMCParticleToPostBremsstrahlungHitList);

    /**
     *  @brief  Identify the reconstructable post-bremsstrahlung radiation hits
     *
     *  @param  leadingMCParticleToPostBremsstrahlungHitList the mapping of leading MCParticles to post-bremsstrahlung hits
     *  @param  maxBremsstrahlungSeparation the maximum separation of a reconstructable post-bremsstrahlung hit from the pre-radiation hits
     *  @param  leadingMCToTrueHitListMap the mapping of cosmic ray and leading MCParticles to their reconstructable hits
     */
    static void AddInPostBremsstrahlungHits(const LeadingMCParticleToPostBremsstrahlungHitList &leadingMCParticleToPostBremsstrahlungHitList, const float maxBremsstrahlungSeparation,
        LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap);

    /**
     *  @brief  Identify the reconstructable post-bremsstrahlung radiation hits
     *
     *  @param  pLeadingParticle the address of the input leading MCParticle
     *  @param  leadingMCParticleToPostBremsstrahlungHitList the mapping of leading MCParticles to post-bremsstrahlung hits
     *  @param  maxBremsstrahlungSeparation the maximum separation of a reconstructable post-bremsstrahlung hit from the pre-radiation hits
     *  @param  leadingMCToTrueHitListMap the mapping of cosmic ray and leading MCParticles to their reconstructable hits
     *  @param  tpcView the TPC view of the hits to be collected
     */
    static void AddInPostBremsstrahlungHits(const pandora::MCParticle *const pLeadingParticle, const LeadingMCParticleToPostBremsstrahlungHitList &leadingMCParticleToPostBremsstrahlungHitList, 
        const float maxBremsstrahlungSeparation, LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const pandora::HitType &tpcView);

    /**
     *  @brief  Select all tier 0 and tier 1 MCParticles in cosmic ray hierarchies from an input list
     *
     *  @param  pMCParticle the address of the input MCParticle list
     *  @param  selectedParticles the output vector of selected MCParticles
     */
    static void SelectLeadingMCParticles(const pandora::MCParticleList *pMCParticleList, pandora::MCParticleVector &selectedParticles);
};

} // namespace lar_content

#endif // #ifndef LAR_MUON_LEADING_HELPER_H
