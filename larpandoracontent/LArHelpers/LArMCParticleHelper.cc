/**
 *  @file   larpandoracontent/LArHelpers/LArMCParticleHelper.cc
 *
 *  @brief  Implementation of the lar monte carlo particle helper class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/MCParticle.h"

#include "Pandora/PdgTable.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <algorithm>
#include <cstdlib>

namespace lar_content
{

using namespace pandora;

LArMCParticleHelper::PrimaryParameters::PrimaryParameters() :
    m_minPrimaryGoodHits(15),
    m_minHitsForGoodView(5),
    m_minPrimaryGoodViews(2),
    m_selectInputHits(true),
    m_foldToPrimaries(true),
    m_maxPhotonPropagation(2.5f),
    m_minHitSharingFraction(0.9f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
  void LArMCParticleHelper::GetMCToPfoCompletenessPurityMaps(const MCContributionMap& mcParticleToHitsMap, const PfoContributionMap& pfoToHitsMap, const MCParticleToPfoHitSharingMap& mcParticleToPfoHitSharingMap, MCParticleToPfoCompletenessPurityMap& mcParticleToPfoCompletenessMap, MCParticleToPfoCompletenessPurityMap& mcParticleToPfoPurityMap) {

    for(const MCParticleToPfoHitSharingMap::value_type &entry : mcParticleToPfoHitSharingMap)
	{
		const PfoToSharedHitsVector pfoToSharedHitsVector(entry.second);
      	PfoToCompletenessPurityVector pfoToCompletenessVector;
      	PfoToCompletenessPurityVector pfoToPurityVector;
      	
		if(mcParticleToHitsMap.find(entry.first) == mcParticleToHitsMap.end()) continue;

      	const unsigned int totMCParticleHits(mcParticleToHitsMap.at(entry.first).size());

		for(const PfoCaloHitListPair &pfoSharedHits : pfoToSharedHitsVector) 
		{
			if(pfoToHitsMap.find(pfoSharedHits.first) == pfoToHitsMap.end()) continue;
			const unsigned int sharedHits(pfoSharedHits.second.size());
			const unsigned int totPfoHits(pfoToHitsMap.at(pfoSharedHits.first).size()); 
			const double completeness(static_cast<double>(sharedHits)/static_cast<double>(totMCParticleHits));
			const double purity(static_cast<double>(sharedHits)/static_cast<double>(totPfoHits));
			pfoToCompletenessVector.push_back({pfoSharedHits.first, completeness});
			pfoToPurityVector.push_back({pfoSharedHits.first, purity});
      	}

		mcParticleToPfoCompletenessMap[entry.first] = pfoToCompletenessVector;
      	mcParticleToPfoPurityMap[entry.first] = pfoToPurityVector;
    	}
		return;
  	}


//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsBeamNeutrinoFinalState(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && LArMCParticleHelper::IsNeutrino(pParentMCParticle));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsDownstreamOfBeamNeutrino(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
    return (LArMCParticleHelper::IsNeutrino(pParentMCParticle));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsTriggeredBeamParticle(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle)); 
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && (nuance == 2001));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsBeamParticle(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && ((nuance == 2000) || (nuance == 2001)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsLeadingBeamParticle(const MCParticle *const pMCParticle)
{
    // ATTN: Only the parent triggered beam particle has nuance code 2001
    const int parentNuance(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCParticle)));
    return (LArMCParticleHelper::IsLeading(pMCParticle) && (parentNuance == 2001 || parentNuance == 2000));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsCosmicRay(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && ((nuance == 3000) || ((nuance == 0) && !LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArMCParticleHelper::GetNuanceCode(const MCParticle *const pMCParticle)
{
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));
    if (pLArMCParticle)
        return pLArMCParticle->GetNuanceCode();

    std::cout << "LArMCParticleHelper::GetNuanceCode - Error: Can't cast to LArMCParticle" << std::endl;
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrino(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    if ((nuance == 0) || (nuance == 2000) || (nuance == 2001) || (nuance == 3000))
        return false;

    const int absoluteParticleId(std::abs(pMCParticle->GetParticleId()));
    return ((NU_E == absoluteParticleId) || (NU_MU == absoluteParticleId) || (NU_TAU == absoluteParticleId));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsPrimary(const pandora::MCParticle *const pMCParticle)
{
    try
    {
        return (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle);
    }
    catch (const StatusCodeException &) {}

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsLeading(const pandora::MCParticle *const pMCParticle)
{
    try
    {
        return (LArMCParticleHelper::GetLeadingMCParticle(pMCParticle) == pMCParticle);
    }
    catch (const StatusCodeException &) {}

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArMCParticleHelper::GetHierarchyTier(const pandora::MCParticle *const pMCParticle)
{
    const MCParticle *pParentMCParticle = pMCParticle;
    int tier(0);

    while (pParentMCParticle->GetParentList().empty() == false)
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        ++tier;
    }

    return tier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsVisible(const MCParticle *const pMCParticle)
{
    const int absoluteParticleId(std::abs(pMCParticle->GetParticleId()));

    if ((E_MINUS == absoluteParticleId) || (MU_MINUS == absoluteParticleId) ||
        (PI_PLUS == absoluteParticleId) || (K_PLUS == absoluteParticleId) ||
        (SIGMA_MINUS == absoluteParticleId) || (SIGMA_PLUS == absoluteParticleId) || (HYPERON_MINUS == absoluteParticleId) ||
        (PROTON == absoluteParticleId) || (PHOTON == absoluteParticleId) ||
        (NEUTRON == absoluteParticleId))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTrueNeutrinos(const MCParticleList *const pMCParticleList, MCParticleVector &trueNeutrinos)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(pMCParticle))
            trueNeutrinos.push_back(pMCParticle);
    }

    std::sort(trueNeutrinos.begin(), trueNeutrinos.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTrueTestBeamParticles(const MCParticleList *const pMCParticleList, MCParticleVector &trueTestBeamParticles)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsTriggeredBeamParticle(pMCParticle))
            trueTestBeamParticles.push_back(pMCParticle);
    }

    std::sort(trueTestBeamParticles.begin(), trueTestBeamParticles.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetParentMCParticle(const MCParticle *const pMCParticle)
{
    const MCParticle *pParentMCParticle = pMCParticle;

    while (pParentMCParticle->GetParentList().empty() == false)
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
    }

    return pParentMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPrimaryMCParticleList(const MCParticleList *const pMCParticleList, MCParticleVector &mcPrimaryVector)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsPrimary(pMCParticle))
            mcPrimaryVector.push_back(pMCParticle);
    }

    std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetLeadingMCParticleList(const MCParticleList *const pMCParticleList, MCParticleVector &mcLeadingVector)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        const bool isBeamParticle(LArMCParticleHelper::IsBeamParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle)));

        if ((isBeamParticle && LArMCParticleHelper::IsLeading(pMCParticle)) || (!isBeamParticle && LArMCParticleHelper::IsPrimary(pMCParticle)))
        {
            mcLeadingVector.push_back(pMCParticle);
        }
    }

    std::sort(mcLeadingVector.begin(), mcLeadingVector.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetPrimaryMCParticle(const MCParticle *const pMCParticle)
{
    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcVector.push_back(pParentMCParticle);
    }

    // Navigate downward through MC parent/daughter links - return the first long-lived charged particle
    for (MCParticleVector::const_reverse_iterator iter = mcVector.rbegin(), iterEnd = mcVector.rend(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pNextParticle = *iter;

        if (LArMCParticleHelper::IsVisible(pNextParticle))
            return pNextParticle;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetLeadingMCParticle(const MCParticle *const pMCParticle, const int hierarchyTierLimit)
{
    // ATTN: If not beam particle return primary particle
    const bool isBeamParticle(LArMCParticleHelper::IsBeamParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle)));

    if (!isBeamParticle)
        return LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);

    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcVector.push_back(pParentMCParticle);
    }

    int hierarchyTier(0);
    const MCParticle *pLeadingMCParticle(nullptr);

    // Navigate downward through MC parent/daughter links - return the first long-lived charged particle
    for (MCParticleVector::const_reverse_iterator iter = mcVector.rbegin(), iterEnd = mcVector.rend(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pNextParticle = *iter;

        // ATTN: Squash any invisible particles (e.g. pizero)
        if (!LArMCParticleHelper::IsVisible(pNextParticle))
            continue;

        if (hierarchyTier <= hierarchyTierLimit)
            pLeadingMCParticle = pNextParticle;

        hierarchyTier++;
    }

    if (!pLeadingMCParticle)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pLeadingMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCPrimaryMap(const MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        try
        {
            const MCParticle *const pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
            mcPrimaryMap[pMCParticle] = pPrimaryMCParticle;
        }
        catch (const StatusCodeException &) {}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------



void LArMCParticleHelper::GetMCLeadingMap(const MCParticleList *const pMCParticleList, MCRelationMap &mcLeadingMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        try
        {
            const MCParticle *const pLeadingMCParticle = LArMCParticleHelper::GetLeadingMCParticle(pMCParticle);
            mcLeadingMap[pMCParticle] = pLeadingMCParticle;
        }
        catch (const StatusCodeException &) {}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

  void LArMCParticleHelper::GetMCToSelfMap(const MCParticleList *const pMCParticleList, MCRelationMap &mcToSelfMap)
{
  for(const MCParticle *const pMCParticle : *pMCParticleList)
  {
    mcToSelfMap[pMCParticle] = pMCParticle;
  }    

}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetMainMCParticle(const ParticleFlowObject *const pPfo)
{
    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
    return MCParticleHelper::GetMainMCParticle(&clusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::SortByMomentum(const MCParticle *const pLhs, const MCParticle *const pRhs)
{
    // Sort by momentum (prefer higher momentum)
    const float momentumLhs(pLhs->GetMomentum().GetMagnitudeSquared());
    const float momentumRhs(pRhs->GetMomentum().GetMagnitudeSquared());

    if (std::fabs(momentumLhs - momentumRhs) > std::numeric_limits<float>::epsilon())
        return (momentumLhs > momentumRhs);

    // Sort by energy (prefer lighter particles)
    if (std::fabs(pLhs->GetEnergy() - pRhs->GetEnergy()) > std::numeric_limits<float>::epsilon())
        return (pLhs->GetEnergy() < pRhs->GetEnergy());

    // Sort by PDG code (prefer smaller numbers)
    if (pLhs->GetParticleId() != pRhs->GetParticleId())
        return (pLhs->GetParticleId() < pRhs->GetParticleId());

    // Sort by vertex position (tie-breaker)
    const float positionLhs(pLhs->GetVertex().GetMagnitudeSquared());
    const float positionRhs(pRhs->GetVertex().GetMagnitudeSquared());

    return (positionLhs < positionRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCParticleToCaloHitMatches(const CaloHitList *const pCaloHitList, const MCRelationMap &mcToPrimaryMCMap,
    CaloHitToMCMap &hitToMCMap, MCContributionMap &mcToTrueHitListMap)
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            const MCParticle *pPrimaryParticle(pHitParticle);

            // ATTN Do not map back to mc primaries if mc to primary mc map not provided
            if (!mcToPrimaryMCMap.empty())
            {
                MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

                if (mcToPrimaryMCMap.end() == mcIter)
                    continue;
		

                pPrimaryParticle = mcIter->second;
            }

            mcToTrueHitListMap[pPrimaryParticle].push_back(pCaloHit);
            hitToMCMap[pCaloHit] = pPrimaryParticle;
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectReconstructableMCParticles(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, const PrimaryParameters &parameters,
    std::function<bool(const MCParticle *const)> fCriteria, MCContributionMap &selectedMCParticlesToHitsMap)
{

    // Obtain map: [MC particle -> self] (to prevent folding to primary MC particle)
    // or [MC particle -> primary mc particle] (to fold to primary MC particle)
    LArMCParticleHelper::MCRelationMap mcToTargetMCMap;
    parameters.m_foldToPrimaries ? LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToTargetMCMap) : LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, mcToTargetMCMap);


    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    // Both cases of m_foldToPrimaries need to input [MC -> MC primary] to complete the above
    /**
    CaloHitList selectedCaloHitList;
    if(parameters.m_foldToPrimaries) {
      LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToTargetMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation); 
    } else {
      LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
      LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);
      LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToPrimaryMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);
    }
    */

    // Obtain maps: [hit -> MC particle (either primary or a downs/storage/epp2/phrwdg/Dune/newVarPandora/PandoraPFA/LArContent-feature/Daresbury/larpandoracontent/LArHelpers/LArMCParticleHelper.cc:378:9: error: expected ‘;’ before ‘}’ token
    CaloHitToMCMap hitToTargetMCMap;
    MCContributionMap targetMCToTrueHitListMap;
    //LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToTargetMCMap, hitToTargetMCMap, targetMCToTrueHitListMap);
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToTargetMCMap, hitToTargetMCMap, targetMCToTrueHitListMap);


    // Obtain vector: all or primary mc particles
    MCParticleVector targetMCVector;
    if(parameters.m_foldToPrimaries){ 
      //std::cout << "TRUE" << std::endl;
      LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, targetMCVector);
    } else {
      //std::cout << "FALSE" << std::endl;
      std::copy(pMCParticleList->begin(), pMCParticleList->end(), std::back_inserter(targetMCVector));
    }

    // Select MCParticles matching criteria
    MCParticleVector candidateTargets;
    LArMCParticleHelper::SelectParticlesMatchingCriteria(targetMCVector, fCriteria, candidateTargets);
    
    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(candidateTargets, targetMCToTrueHitListMap, mcToTargetMCMap, parameters, selectedMCParticlesToHitsMap);

}


//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectReconstructableTestBeamHierarchyMCParticles(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, const PrimaryParameters &parameters,
    std::function<bool(const MCParticle *const)> fCriteria, MCContributionMap &selectedMCParticlesToHitsMap)
{
    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCLeadingMap(pMCParticleList, mcToPrimaryMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList selectedCaloHitList;
    LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToPrimaryMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);

    // Obtain maps: [hit -> primary mc particle], [primary mc particle -> list of hits]
    CaloHitToMCMap hitToPrimaryMCMap;
    MCContributionMap mcToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToPrimaryMCMap, hitToPrimaryMCMap, mcToTrueHitListMap);

    // Obtain vector: primary mc particles
    MCParticleVector mcPrimaryVector;
    LArMCParticleHelper::GetLeadingMCParticleList(pMCParticleList, mcPrimaryVector);

    // Select MCParticles matching criteria
    MCParticleVector candidateTargets;
    LArMCParticleHelper::SelectParticlesMatchingCriteria(mcPrimaryVector, fCriteria, candidateTargets);

    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(candidateTargets, mcToTrueHitListMap, mcToPrimaryMCMap, parameters, selectedMCParticlesToHitsMap);

    // ATTN: The parent particle must be in the hierarchy map, event if not reconstructable
    for (const MCParticle *const pMCParticle : candidateTargets)
    {
        if (!LArMCParticleHelper::IsBeamParticle(pMCParticle))
            continue;

        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
        if (selectedMCParticlesToHitsMap.find(pParentMCParticle) == selectedMCParticlesToHitsMap.end())
        {
            CaloHitList caloHitList;
            selectedMCParticlesToHitsMap.insert(MCContributionMap::value_type(pParentMCParticle, caloHitList));

        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetUnfoldedPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMap &selectedMCParticleToHitsMap, PfoContributionMap &pfoToReconstructable2DHitsMap)
{

  for(const ParticleFlowObject *const pPfo : pfoList) {
    CaloHitList pfoHitList;
    LArMCParticleHelper::CollectReconstructable2DHits(PfoList{pPfo}, {selectedMCParticleToHitsMap}, pfoHitList);

    if (!pfoToReconstructable2DHitsMap.insert(PfoContributionMap::value_type(pPfo, pfoHitList)).second)
      throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
  }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMap &selectedMCParticleToHitsMap,
    PfoContributionMap &pfoToReconstructable2DHitsMap)
{
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(pfoList, MCContributionMapVector({selectedMCParticleToHitsMap}), pfoToReconstructable2DHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTestBeamHierarchyPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMap &selectedMCParticleToHitsMap,
    PfoContributionMap &pfoToReconstructable2DHitsMap)
{
    LArMCParticleHelper::GetTestBeamHierarchyPfoToReconstructable2DHitsMap(pfoList, MCContributionMapVector({selectedMCParticleToHitsMap}), pfoToReconstructable2DHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
    PfoContributionMap &pfoToReconstructable2DHitsMap)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        CaloHitList pfoHitList;
        LArMCParticleHelper::CollectReconstructable2DHits(pPfo, selectedMCParticleToHitsMaps, pfoHitList);

        if (!pfoToReconstructable2DHitsMap.insert(PfoContributionMap::value_type(pPfo, pfoHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTestBeamHierarchyPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
    PfoContributionMap &pfoToReconstructable2DHitsMap)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        CaloHitList pfoHitList;
        LArMCParticleHelper::CollectReconstructableTestBeamHierarchy2DHits(pPfo, selectedMCParticleToHitsMaps, pfoHitList);

        if (!pfoToReconstructable2DHitsMap.insert(PfoContributionMap::value_type(pPfo, pfoHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(const PfoContributionMap &pfoToReconstructable2DHitsMap, const MCContributionMapVector &selectedMCParticleToHitsMaps,
    PfoToMCParticleHitSharingMap &pfoToMCParticleHitSharingMap, MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap)
{
    PfoVector sortedPfos;
    for (const auto &mapEntry : pfoToReconstructable2DHitsMap) sortedPfos.push_back(mapEntry.first);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo : sortedPfos)
    {
        for (const MCContributionMap &mcParticleToHitsMap : selectedMCParticleToHitsMaps)
        {
            MCParticleVector sortedMCParticles;
            for (const auto &mapEntry : mcParticleToHitsMap) sortedMCParticles.push_back(mapEntry.first);
            std::sort(sortedMCParticles.begin(), sortedMCParticles.end(), PointerLessThan<MCParticle>());

            for (const MCParticle *const pMCParticle : sortedMCParticles)
            {
                // Add map entries for this Pfo & MCParticle if required
                if (pfoToMCParticleHitSharingMap.find(pPfo) == pfoToMCParticleHitSharingMap.end())
                    if (!pfoToMCParticleHitSharingMap.insert(PfoToMCParticleHitSharingMap::value_type(pPfo, MCParticleToSharedHitsVector())).second)
                        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT); // ATTN maybe overkill

                if (mcParticleToPfoHitSharingMap.find(pMCParticle) == mcParticleToPfoHitSharingMap.end())
                    if (!mcParticleToPfoHitSharingMap.insert(MCParticleToPfoHitSharingMap::value_type(pMCParticle, PfoToSharedHitsVector())).second)
                        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                // Check this Pfo & MCParticle pairing hasn't already been checked
                MCParticleToSharedHitsVector &mcHitPairs(pfoToMCParticleHitSharingMap.at(pPfo));
                PfoToSharedHitsVector &pfoHitPairs(mcParticleToPfoHitSharingMap.at(pMCParticle));

                if (std::any_of(mcHitPairs.begin(), mcHitPairs.end(), [&] (const MCParticleCaloHitListPair &pair) { return (pair.first == pMCParticle); }))
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                if (std::any_of(pfoHitPairs.begin(), pfoHitPairs.end(), [&] (const PfoCaloHitListPair &pair) { return (pair.first == pPfo); }))
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                // Add records to maps if there are any shared hits
                const CaloHitList sharedHits(LArMCParticleHelper::GetSharedHits(pfoToReconstructable2DHitsMap.at(pPfo), mcParticleToHitsMap.at(pMCParticle)));

                if (!sharedHits.empty())
                {
                    mcHitPairs.push_back(MCParticleCaloHitListPair(pMCParticle, sharedHits));
                    pfoHitPairs.push_back(PfoCaloHitListPair(pPfo, sharedHits));

                    std::sort(mcHitPairs.begin(), mcHitPairs.end(), [] (const MCParticleCaloHitListPair &a, const MCParticleCaloHitListPair &b) -> bool {
                        return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size() : LArMCParticleHelper::SortByMomentum(a.first, b.first)); });

                    std::sort(pfoHitPairs.begin(), pfoHitPairs.end(), [] (const PfoCaloHitListPair &a, const PfoCaloHitListPair &b) -> bool {
                        return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size() : LArPfoHelper::SortByNHits(a.first, b.first)); });
                }
            }
        }
    }
}

// private
//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CollectReconstructable2DHits(const ParticleFlowObject *const pPfo, const MCContributionMapVector &selectedMCParticleToHitsMaps,
    pandora::CaloHitList &reconstructableCaloHitList2D)
{
    // Collect all 2D calo hits in pfo hierarchy
    PfoList pfoList;
    LArPfoHelper::GetAllDownstreamPfos(pPfo, pfoList);

    LArMCParticleHelper::CollectReconstructable2DHits(pfoList, selectedMCParticleToHitsMaps, reconstructableCaloHitList2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CollectReconstructableTestBeamHierarchy2DHits(const ParticleFlowObject *const pPfo, const MCContributionMapVector &selectedMCParticleToHitsMaps,
    pandora::CaloHitList &reconstructableCaloHitList2D)
{
    // Collect all 2D calo hits in pfo hierarchy
    PfoList pfoList;

    // ATTN: Only collect downstream pfos for daughter test beam particles & cosmics
    if (pPfo->GetParentPfoList().empty() && LArPfoHelper::IsTestBeam(pPfo))
    {
        pfoList.push_back(pPfo);
    }
    else
    {
        LArPfoHelper::GetAllDownstreamPfos(pPfo, pfoList);
    }

    LArMCParticleHelper::CollectReconstructable2DHits(pfoList, selectedMCParticleToHitsMaps, reconstructableCaloHitList2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CollectReconstructable2DHits(const PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToHitsMaps,
    pandora::CaloHitList &reconstructableCaloHitList2D)
{
    CaloHitList caloHitList2D;
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_U, caloHitList2D);
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_V, caloHitList2D);
    LArPfoHelper::GetCaloHits(pfoList, TPC_VIEW_W, caloHitList2D);
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_U, caloHitList2D); // TODO check isolated usage throughout
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_V, caloHitList2D);
    LArPfoHelper::GetIsolatedCaloHits(pfoList, TPC_VIEW_W, caloHitList2D);

    // Filter for only reconstructable hits
    for (const CaloHit *const pCaloHit : caloHitList2D)
    {
        bool isTargetHit(false);
        for (const MCContributionMap &mcParticleToHitsMap : selectedMCParticleToHitsMaps)
        { 
            // ATTN This map is unordered, but this does not impact search for specific target hit
            for (const MCContributionMap::value_type &mapEntry : mcParticleToHitsMap)
            {
                if (std::find(mapEntry.second.begin(), mapEntry.second.end(), pCaloHit) != mapEntry.second.end())
                {
                    isTargetHit = true;
                    break;
                }
            }
            if (isTargetHit) break;
        }

        if (isTargetHit)
            reconstructableCaloHitList2D.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    CaloHitList &selectedCaloHitList, const bool selectInputHits, const float maxPhotonPropagation)
{
    if (!selectInputHits)
    {
        selectedCaloHitList.insert(selectedCaloHitList.end(), pCaloHitList->begin(), pCaloHitList->end());
        return;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

            if (mcToPrimaryMCMap.end() == mcIter)
                continue;
	    

            const MCParticle *const pPrimaryParticle = mcIter->second;

            if (PassMCParticleChecks(pPrimaryParticle, pPrimaryParticle, pHitParticle, maxPhotonPropagation))
                selectedCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectGoodCaloHits(const CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    CaloHitList &selectedGoodCaloHitList, const bool selectInputHits, const float minHitSharingFraction)
{
    if (!selectInputHits)
    {
        selectedGoodCaloHitList.insert(selectedGoodCaloHitList.end(), pSelectedCaloHitList->begin(), pSelectedCaloHitList->end());
        return;
    }

    for (const CaloHit *const pCaloHit : *pSelectedCaloHitList)
    {
        MCParticleVector mcParticleVector;
        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap()) mcParticleVector.push_back(mapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        MCParticleWeightMap primaryWeightMap;

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pMCParticle);

            if (mcToPrimaryMCMap.end() != mcIter)
                primaryWeightMap[mcIter->second] += weight;
        }

        MCParticleVector mcPrimaryVector;
        for (const auto &mapEntry : primaryWeightMap) mcPrimaryVector.push_back(mapEntry.first);
        std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), PointerLessThan<MCParticle>());

        const MCParticle *pBestPrimaryParticle(nullptr);
        float bestPrimaryWeight(0.f), primaryWeightSum(0.f);

        for (const MCParticle *const pPrimaryMCParticle : mcPrimaryVector)
        {
            const float primaryWeight(primaryWeightMap.at(pPrimaryMCParticle));
            primaryWeightSum += primaryWeight;

            if (primaryWeight > bestPrimaryWeight)
            {
                bestPrimaryWeight = primaryWeight;
                pBestPrimaryParticle = pPrimaryMCParticle;
            }
        }

        if (!pBestPrimaryParticle || (primaryWeightSum < std::numeric_limits<float>::epsilon()) || ((bestPrimaryWeight / primaryWeightSum) < minHitSharingFraction))
            continue;

        selectedGoodCaloHitList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectParticlesMatchingCriteria(const MCParticleVector &inputMCParticles, std::function<bool(const MCParticle *const)> fCriteria,
    MCParticleVector &selectedParticles)
{
    for (const MCParticle *const pMCParticle : inputMCParticles)
    {
        if (fCriteria(pMCParticle))
            selectedParticles.push_back(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectParticlesByHitCount(const MCParticleVector &candidateTargets, const MCContributionMap &mcToTrueHitListMap,
    const MCRelationMap &mcToPrimaryMCMap, const PrimaryParameters &parameters, MCContributionMap &selectedMCParticlesToHitsMap)
{
    // Apply restrictions on the number of good hits associated with the MCParticles
    for (const MCParticle * const pMCTarget : candidateTargets)
    {
        MCContributionMap::const_iterator trueHitsIter = mcToTrueHitListMap.find(pMCTarget);
        if (mcToTrueHitListMap.end() == trueHitsIter)
            continue;

        const CaloHitList &caloHitList(trueHitsIter->second);

        // Remove shared hits where target particle deposits below threshold energy fraction
        CaloHitList goodCaloHitList;
        LArMCParticleHelper::SelectGoodCaloHits(&caloHitList, mcToPrimaryMCMap, goodCaloHitList, parameters.m_selectInputHits, parameters.m_minHitSharingFraction);

        if (goodCaloHitList.size() < parameters.m_minPrimaryGoodHits)
            continue;

        unsigned int nGoodViews(0);
        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, goodCaloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, goodCaloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, goodCaloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (nGoodViews < parameters.m_minPrimaryGoodViews)
            continue;

        if (!selectedMCParticlesToHitsMap.insert(MCContributionMap::value_type(pMCTarget, caloHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::PassMCParticleChecks(const MCParticle *const pOriginalPrimary, const MCParticle *const pThisMCParticle,
    const MCParticle *const pHitMCParticle, const float maxPhotonPropagation)
{
  //if (NEUTRON == std::abs(pThisMCParticle->GetParticleId()))
  //return false;

  if ((PHOTON == pThisMCParticle->GetParticleId()) && (PHOTON != GetPrimaryMCParticle(pThisMCParticle)->GetParticleId()) && (E_MINUS != std::abs(GetPrimaryMCParticle(pThisMCParticle)->GetParticleId())))
    {
      if ((pThisMCParticle->GetEndpoint() - pThisMCParticle->GetVertex()).GetMagnitude() > maxPhotonPropagation)
	return false;
    }

    if (pThisMCParticle == pHitMCParticle)
        return true;

    for (const MCParticle *const pDaughterMCParticle : pThisMCParticle->GetDaughterList())
    {
        if (PassMCParticleChecks(pOriginalPrimary, pDaughterMCParticle, pHitMCParticle, maxPhotonPropagation))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArMCParticleHelper::GetSharedHits(const CaloHitList &hitListA, const CaloHitList &hitListB)
{
    CaloHitList sharedHits;

    for (const CaloHit *const pCaloHit : hitListA)
    {
        if (std::find(hitListB.begin(), hitListB.end(), pCaloHit) != hitListB.end())
            sharedHits.push_back(pCaloHit);
    }

    return sharedHits;
}

} // namespace lar_content
