/**
 *  @file   larpandoracontent/LArMonitoring/MuonLeadingEventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the muon leading event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"

#include "larpandoracontent/LArMonitoring/MuonLeadingEventValidationAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

MuonLeadingEventValidationAlgorithm::MuonLeadingEventValidationAlgorithm() :
    m_removeRecoMuonHits(false),    
    m_deltaRayMode(false),
    m_michelMode(false),
    m_muonsToSkip(0),
    m_visualize(false),
    m_ignoreIncorrectMuons(false),
    m_writeRawMatchesToTree(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MuonLeadingEventValidationAlgorithm::~MuonLeadingEventValidationAlgorithm()
{ 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{    
    this->FillDeltaRayValidationInfo(pMCParticleList, pCaloHitList, pPfoList, validationInfo);

    if (m_ignoreIncorrectMuons)
    {
        MCParticleList incorrectlyReconstructedCosmicRays;
        this->DetermineIncorrectlyReconstructedMuons(pMCParticleList, pCaloHitList, pPfoList, incorrectlyReconstructedCosmicRays);

        LArMCParticleHelper::MCContributionMap allMCToHitsMap(validationInfo.GetAllMCParticleToHitsMap());
        LArMCParticleHelper::MCContributionMap targetMCToHitsMap(validationInfo.GetTargetMCParticleToHitsMap());
        LArMCParticleHelper::MCParticleToPfoHitSharingMap matchingMap(validationInfo.GetMCToPfoHitSharingMap());
        LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMatchingMap(validationInfo.GetInterpretedMCToPfoHitSharingMap());

        for (const MCParticle *const pIncorrectCosmicRay : incorrectlyReconstructedCosmicRays)
        {
            auto allIter(allMCToHitsMap.find(pIncorrectCosmicRay));

            if (allIter != allMCToHitsMap.end())
                allMCToHitsMap.erase(allIter);

            auto targetIter(targetMCToHitsMap.find(pIncorrectCosmicRay));

            if (targetIter != targetMCToHitsMap.end())
                targetMCToHitsMap.erase(targetIter);

            auto matchingIter(matchingMap.find(pIncorrectCosmicRay));

            if (matchingIter != matchingMap.end())
                matchingMap.erase(matchingIter);
            
            auto interpretedMatchingIter(interpretedMatchingMap.find(pIncorrectCosmicRay));

            if (interpretedMatchingIter != interpretedMatchingMap.end())
                interpretedMatchingMap.erase(interpretedMatchingIter);
        }

        validationInfo.SetAllMCParticleToHitsMap(allMCToHitsMap);
        validationInfo.SetTargetMCParticleToHitsMap(targetMCToHitsMap);
        validationInfo.SetMCToPfoHitSharingMap(matchingMap);
        validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMatchingMap);
    }
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::FillDeltaRayValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    CaloHitList recoMuonHitList;
    if (m_removeRecoMuonHits)
        this->GetRecoMuonHits(pMCParticleList, pCaloHitList, pPfoList, recoMuonHitList);
    
    if (pMCParticleList && pCaloHitList)
    {
        // Get reconstructable MCParticle hit ownership map (non-muon leading hierarchy is folded whilst muon is unfolded)
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;        
        LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, m_validationParameters, recoMuonHitList,
            targetMCParticleToHitsMap, this->GetPandora());    

        // Do not change the hit share fraction (these hits are classed as 'ambiguous' and should not be in the metrics)
        LArMuonLeadingHelper::ValidationParameters allValidationParams(m_validationParameters);
        allValidationParams.m_minPrimaryGoodHits = 0;
        allValidationParams.m_minHitsForGoodView = 0;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, allValidationParams, recoMuonHitList,
            allMCParticleToHitsMap, this->GetPandora());
        
        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }

    if (pPfoList)
    {        
        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, false);

        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }
    
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::DetermineIncorrectlyReconstructedMuons(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, MCParticleList &incorrectlyReconstructedCosmicRays) const
{
    // Perform cosmic ray matching
    ValidationInfo validationInfo;
    
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_validationParameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        // Since muon reco, set sharing hit fraction to 0 (ask why john did this for uboone?)
        LArMCParticleHelper::PrimaryParameters parameters(m_validationParameters);
        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);

        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }

    if (pPfoList)
    {
        PfoList allConnectedPfos;
        LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

        PfoList finalStatePfos;
        for (const ParticleFlowObject *const pPfo : allConnectedPfos)
        {
            if (pPfo->GetParentPfoList().empty())
                finalStatePfos.push_back(pPfo);
        }

        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, m_validationParameters.m_foldBackHierarchy);

        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);

    // Determine correctly reconstructed cosmic rays
    MCParticleVector targetCosmicRayVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, targetCosmicRayVector);

    for (const MCParticle *const pTargetCosmicRay : targetCosmicRayVector)
    {
        if (mcToPfoHitSharingMap.at(pTargetCosmicRay).empty())
        {
            incorrectlyReconstructedCosmicRays.push_back(pTargetCosmicRay);
            continue;
        }

        const CaloHitList &mcHitList(validationInfo.GetAllMCParticleToHitsMap().at(pTargetCosmicRay));

        unsigned int nAboveThresholdMatches(0);
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pTargetCosmicRay))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));
            
            const bool isGoodMatch(this->IsGoodMatch(mcHitList, pfoHitList, sharedHitList));

            if (isGoodMatch)
                ++nAboveThresholdMatches;
        }

        if (nAboveThresholdMatches != 1)
            incorrectlyReconstructedCosmicRays.push_back(pTargetCosmicRay);
    }

    /////////////////////////
    /*
    for (const MCParticle *const pCosmicRay : targetCosmicRayVector)
    {
        bool isCorrect(std::find(correctlyReconstructedCosmicRays.begin(), correctlyReconstructedCosmicRays.end(), pCosmicRay) != correctlyReconstructedCosmicRays.end());

        if (mcToPfoHitSharingMap.at(pCosmicRay).empty())
        {
            std::cout << "No Matches: " << std::endl;
            std::cout << "Correct: " << isCorrect << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }

        const CaloHitList &mcHitList(validationInfo.GetAllMCParticleToHitsMap().at(pCosmicRay));

        unsigned int nAboveThresholdMatches(0);
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pCosmicRay))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));
            
            const bool isGoodMatch(this->IsGoodMatch(mcHitList, pfoHitList, sharedHitList));

            if (isGoodMatch)
                ++nAboveThresholdMatches;

            if (isGoodMatch)
            {
                PfoList matchedPfo({pfoToSharedHits.first});
                PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &matchedPfo, "Matched Pfo", BLUE);
            }
        }

        std::cout << nAboveThresholdMatches << " above threshold matches" << std::endl;
        std::cout << "Correct: " << isCorrect << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
    */
    /////////////////////////    
}       

//------------------------------------------------------------------------------------------------------------------------------------------
    
void MuonLeadingEventValidationAlgorithm::ProcessOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const
{
    // Folded hit ownership/sharing maps for leading muon ionisation particles
    const LArMCParticleHelper::MCContributionMap &foldedAllMCToHitsMap(validationInfo.GetAllMCParticleToHitsMap());
    const LArMCParticleHelper::MCContributionMap &foldedTargetMCToHitsMap(validationInfo.GetTargetMCParticleToHitsMap());
    const LArMCParticleHelper::PfoContributionMap &foldedPfoToHitsMap(validationInfo.GetPfoToHitsMap());
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &foldedMCToPfoHitSharingMap((fillTree && m_writeRawMatchesToTree) ? validationInfo.GetMCToPfoHitSharingMap() :
        useInterpretedMatching ? validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());
    
    // Consider only delta rays from reconstructable CR muons
    MCParticleVector mcCRVector;
    for (auto &entry : foldedTargetMCToHitsMap)
    {
        if (LArMCParticleHelper::IsCosmicRay(entry.first))
            mcCRVector.push_back(entry.first);
    }
    std::sort(mcCRVector.begin(), mcCRVector.end(), LArMCParticleHelper::SortByMomentum);

    // Process matches
    int muonCount(0);
    int totalReconstructableCRLs(0);
    std::stringstream stringStream;
    for (const MCParticle *const pCosmicRay : mcCRVector)
    {
        // Cosmic ray parameters
        int ID_CR, nReconstructableChildCRLs(0);
        float mcE_CR, mcPX_CR, mcPY_CR, mcPZ_CR;
        int nMCHitsTotal_CR, nMCHitsU_CR, nMCHitsV_CR, nMCHitsW_CR;
        float mcVertexX_CR, mcVertexY_CR, mcVertexZ_CR, mcEndX_CR, mcEndY_CR, mcEndZ_CR;
        int nCorrectChildCRLs(0);
    
        // Leading particle parameters
        FloatVector mcE_CRL, mcPX_CRL, mcPY_CRL, mcPZ_CRL;
        IntVector ID_CRL;        
        IntVector nMCHitsTotal_CRL, nMCHitsU_CRL, nMCHitsV_CRL, nMCHitsW_CRL;
        FloatVector mcVertexX_CRL, mcVertexY_CRL, mcVertexZ_CRL, mcEndX_CRL, mcEndY_CRL, mcEndZ_CRL;
        IntVector nAboveThresholdMatches_CRL, isCorrect_CRL, isCorrectParentLink_CRL;
        IntVector bestMatchNHitsTotal_CRL, bestMatchNHitsU_CRL, bestMatchNHitsV_CRL, bestMatchNHitsW_CRL;
        IntVector bestMatchNSharedHitsTotal_CRL, bestMatchNSharedHitsU_CRL, bestMatchNSharedHitsV_CRL, bestMatchNSharedHitsW_CRL;
        IntVector bestMatchNParentTrackHitsTotal_CRL, bestMatchNParentTrackHitsU_CRL, bestMatchNParentTrackHitsV_CRL, bestMatchNParentTrackHitsW_CRL;
        IntVector bestMatchNOtherTrackHitsTotal_CRL, bestMatchNOtherTrackHitsU_CRL, bestMatchNOtherTrackHitsV_CRL, bestMatchNOtherTrackHitsW_CRL;
        IntVector bestMatchNOtherShowerHitsTotal_CRL, bestMatchNOtherShowerHitsU_CRL, bestMatchNOtherShowerHitsV_CRL, bestMatchNOtherShowerHitsW_CRL;
        IntVector totalCRLHitsInBestMatchParentCR_CRL, uCRLHitsInBestMatchParentCR_CRL, vCRLHitsInBestMatchParentCR_CRL, wCRLHitsInBestMatchParentCR_CRL;
        FloatVector bestMatchVertexX_CRL, bestMatchVertexY_CRL, bestMatchVertexZ_CRL;

        // Contamination parameters
        IntVector  bestMatchOtherShowerHitsID_CRL, bestMatchOtherTrackHitsID_CRL, bestMatchParentTrackHitsID_CRL, bestMatchCRLHitsInCRID_CRL;
        FloatVector bestMatchOtherShowerHitsDistance_CRL, bestMatchOtherTrackHitsDistance_CRL, bestMatchParentTrackHitsDistance_CRL, bestMatchCRLHitsInCRDistance_CRL;

        // Move on if cosmic ray has not been reconstructed
        if (foldedMCToPfoHitSharingMap.at(pCosmicRay).empty())
            continue;

        // Obtain reconstructable leading particles
        MCParticleVector childLeadingParticles;
        for (const MCParticle *const pMuonChild : pCosmicRay->GetDaughterList())
        {
            if (!(m_michelMode || m_deltaRayMode))
                continue;
                        
            if (m_deltaRayMode)
            {
                if (!LArMuonLeadingHelper::IsDeltaRay(pMuonChild))
                    continue;
            }

            if (m_michelMode)
            {
                if (!LArMuonLeadingHelper::IsMichel(pMuonChild))
                    continue;
            }

            // Move on if leading particle is not reconstructable
            LArMCParticleHelper::MCContributionMap::const_iterator iter(foldedTargetMCToHitsMap.find(pMuonChild));
            if (iter == foldedTargetMCToHitsMap.end())
                continue;

            childLeadingParticles.push_back(pMuonChild);
        }

        // Move on if cosmic ray has no leading delta ray child particles
        if (childLeadingParticles.empty())
            continue;
        
        std::sort(childLeadingParticles.begin(), childLeadingParticles.end(), LArMCParticleHelper::SortByMomentum);

        ++muonCount;

        if (muonCount < (m_muonsToSkip))
            continue;

        // Pull cosmic ray info
        const CaloHitList &cosmicRayHitList(foldedAllMCToHitsMap.at(pCosmicRay));

        ///////////////////////////////
        if(m_visualize && useInterpretedMatching)
        {
            std::cout << "MC COSMIC RAY HITS" << std::endl;
            this->PrintHits(cosmicRayHitList, "MC_CR", BLUE);
            
            const CartesianVector vertex(pCosmicRay->GetVertex()), endpoint(pCosmicRay->GetEndpoint());

            //const float x0(cosmicRayHitList.front()->GetX0());
            //const CartesianVector shiftedVertex(vertex.GetX() - x0, vertex.GetY(), vertex.GetZ());
            //const CartesianVector shiftedEndpoint(endpoint.GetX() - x0, vertex.GetY(), vertex.GetZ());
            
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vertex, "vertex ", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint, "endpoint", BLACK, 2);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }
        ///////////////////////////////
        
        ID_CR = muonCount;
        mcE_CR = pCosmicRay->GetEnergy();
        mcPX_CR = pCosmicRay->GetMomentum().GetX();
        mcPY_CR = pCosmicRay->GetMomentum().GetY();
        mcPZ_CR = pCosmicRay->GetMomentum().GetZ();        
        mcVertexX_CR = pCosmicRay->GetVertex().GetX();
        mcVertexY_CR = pCosmicRay->GetVertex().GetY();
        mcVertexZ_CR = pCosmicRay->GetVertex().GetZ();                
        mcEndX_CR = pCosmicRay->GetEndpoint().GetX();
        mcEndY_CR = pCosmicRay->GetEndpoint().GetY();
        mcEndZ_CR = pCosmicRay->GetEndpoint().GetZ();
        nMCHitsTotal_CR = cosmicRayHitList.size();
        nMCHitsU_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, cosmicRayHitList);
        nMCHitsV_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, cosmicRayHitList);        
        nMCHitsW_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, cosmicRayHitList);
        nReconstructableChildCRLs = childLeadingParticles.size();
        
        stringStream << "\033[34m" << "(Parent CR: " << muonCount << ") " << "\033[0m" 
            << "Energy " << pCosmicRay->GetEnergy()
            << ", Dist. " << (pCosmicRay->GetEndpoint() - pCosmicRay->GetVertex()).GetMagnitude()
            << ", nMCHits " << cosmicRayHitList.size()
            << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, cosmicRayHitList)
            << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, cosmicRayHitList)
            << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, cosmicRayHitList) << ")"
            << ", nReconstructableCRLs " << nReconstructableChildCRLs << std::endl;   
        
        // Pull delta ray data
        int leadingCount(0);
        for (const MCParticle *const pLeadingParticle : childLeadingParticles)
        {
            // Pull delta ray MC info
            const CaloHitList &leadingParticleHitList(foldedAllMCToHitsMap.at(pLeadingParticle));
            ++leadingCount;
            ++totalReconstructableCRLs;

            ///////////////////////////////
            if (m_visualize && useInterpretedMatching)
            {
                std::cout << "MC DELTA RAY HITS" << std::endl;
                this->PrintHits(leadingParticleHitList, "MC_DR", RED);

                const CartesianVector vertex(pLeadingParticle->GetVertex()), endpoint(pLeadingParticle->GetEndpoint());

                //const float x0(leadingParticleHitList.front()->GetX0());
                //const CartesianVector shiftedVertex(vertex.GetX() - x0, vertex.GetY(), vertex.GetZ());
                //const CartesianVector shiftedEndpoint(endpoint.GetX() - x0, vertex.GetY(), vertex.GetZ());
                
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vertex, "vertex", BLACK, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint, "endpoint", BLACK, 2);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }
            ///////////////////////////////

            mcE_CRL.push_back(pLeadingParticle->GetEnergy());
            ID_CRL.push_back(leadingCount);
            mcPX_CRL.push_back(pLeadingParticle->GetMomentum().GetX());
            mcPY_CRL.push_back(pLeadingParticle->GetMomentum().GetY());
            mcPZ_CRL.push_back(pLeadingParticle->GetMomentum().GetZ());
            mcVertexX_CRL.push_back(pLeadingParticle->GetVertex().GetX());
            mcVertexY_CRL.push_back(pLeadingParticle->GetVertex().GetY());
            mcVertexZ_CRL.push_back(pLeadingParticle->GetVertex().GetZ());        
            mcEndX_CRL.push_back(pLeadingParticle->GetEndpoint().GetX());
            mcEndY_CRL.push_back(pLeadingParticle->GetEndpoint().GetY());
            mcEndZ_CRL.push_back(pLeadingParticle->GetEndpoint().GetZ());   
            nMCHitsTotal_CRL.push_back(leadingParticleHitList.size());
            nMCHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitList));
            nMCHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitList));   
            nMCHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitList));

            stringStream << "\033[33m"  << "(Child "<< (m_deltaRayMode ? "DR: " : "Michel: ") << leadingCount << ")  " << "\033[0m"
                << "Energy " << pLeadingParticle->GetEnergy()
                << ", Dist. " << (pLeadingParticle->GetEndpoint() - pLeadingParticle->GetVertex()).GetMagnitude()
                << ", nMCHits " << leadingParticleHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitList) << ")" << std::endl;        

            // Look at the pfo matches
            int nMatches(0), nAboveThresholdMatches(0);
            bool isCorrectParentLink(false);
            for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : foldedMCToPfoHitSharingMap.at(pLeadingParticle))
            {
                const ParticleFlowObject *const pMatchedPfo(pfoToSharedHits.first);
                const CaloHitList &pfoHitList(foldedPfoToHitsMap.at(pMatchedPfo));
                const CaloHitList &sharedHitList(pfoToSharedHits.second);

                const bool isGoodMatch(this->IsGoodMatch(leadingParticleHitList, pfoHitList, sharedHitList));

                ++nMatches;

                if (isGoodMatch)
                    ++nAboveThresholdMatches;

                CaloHitList parentTrackHits, otherTrackHits, otherShowerHits;
                LArMuonLeadingHelper::GetPfoMatchContamination(pLeadingParticle, pfoHitList, parentTrackHits, otherTrackHits, otherShowerHits);

                // Check whether the reconstructed pfo has the correct parent-child link
                CaloHitList mcParentMatchedPfoHits;                
                bool isMatchedToCorrectCosmicRay(false);
                const ParticleFlowObject *const pParentPfo(LArPfoHelper::GetParentPfo(pMatchedPfo));

                for (LArMCParticleHelper::PfoToSharedHitsVector::const_iterator cosmicRayMatchedPfoPair = foldedMCToPfoHitSharingMap.at(pCosmicRay).begin();
                    cosmicRayMatchedPfoPair != foldedMCToPfoHitSharingMap.at(pCosmicRay).end(); ++cosmicRayMatchedPfoPair)
                {
                    const ParticleFlowObject *const pCosmicRayPfo(cosmicRayMatchedPfoPair->first);
                    
                    mcParentMatchedPfoHits.insert(mcParentMatchedPfoHits.end(), foldedPfoToHitsMap.at(pCosmicRayPfo).begin(), foldedPfoToHitsMap.at(pCosmicRayPfo).end());
                    
                    if (pCosmicRayPfo == pParentPfo)
                        isMatchedToCorrectCosmicRay = true;
                }                    

                CaloHitList leadingParticleHitsInParentCosmicRay;
                LArMuonLeadingHelper::GetMuonPfoContaminationContribution(mcParentMatchedPfoHits, leadingParticleHitList, leadingParticleHitsInParentCosmicRay);
                
                if ((nAboveThresholdMatches == 1) && isGoodMatch)
                {
                    isCorrectParentLink = isMatchedToCorrectCosmicRay;
                        
                    if (isMatchedToCorrectCosmicRay)
                    {
                        isCorrectParentLink_CRL.push_back(1);
                    }
                    else
                    {
                        isCorrectParentLink_CRL.push_back(0);
                    }
                    
                    bestMatchNHitsTotal_CRL.push_back(pfoHitList.size());
                    bestMatchNHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList));
                    bestMatchNHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList));
                    bestMatchNHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList));
                    
                    bestMatchNSharedHitsTotal_CRL.push_back(sharedHitList.size());
                    bestMatchNSharedHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList));
                    bestMatchNSharedHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList));
                    bestMatchNSharedHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList));

                    bestMatchNParentTrackHitsTotal_CRL.push_back(parentTrackHits.size());
                    bestMatchNParentTrackHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, parentTrackHits));
                    bestMatchNParentTrackHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, parentTrackHits));
                    bestMatchNParentTrackHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, parentTrackHits));

                    bestMatchNOtherTrackHitsTotal_CRL.push_back(otherTrackHits.size());
                    bestMatchNOtherTrackHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, otherTrackHits));
                    bestMatchNOtherTrackHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, otherTrackHits));
                    bestMatchNOtherTrackHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, otherTrackHits));

                    bestMatchNOtherShowerHitsTotal_CRL.push_back(otherShowerHits.size());
                    bestMatchNOtherShowerHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, otherShowerHits));
                    bestMatchNOtherShowerHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, otherShowerHits));
                    bestMatchNOtherShowerHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, otherShowerHits));
                    
                    totalCRLHitsInBestMatchParentCR_CRL.push_back(leadingParticleHitsInParentCosmicRay.size());
                    uCRLHitsInBestMatchParentCR_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitsInParentCosmicRay));
                    vCRLHitsInBestMatchParentCR_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitsInParentCosmicRay));
                    wCRLHitsInBestMatchParentCR_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitsInParentCosmicRay));

                    bestMatchOtherShowerHitsID_CRL.insert(bestMatchOtherShowerHitsID_CRL.end(), otherShowerHits.size(), leadingCount);
                    this->FillContaminationHitsDistance(otherShowerHits, leadingParticleHitList, bestMatchOtherShowerHitsDistance_CRL);

                    bestMatchOtherTrackHitsID_CRL.insert(bestMatchOtherTrackHitsID_CRL.end(), otherTrackHits.size(), leadingCount);
                    this->FillContaminationHitsDistance(otherTrackHits, leadingParticleHitList, bestMatchOtherTrackHitsDistance_CRL);

                    bestMatchParentTrackHitsID_CRL.insert(bestMatchParentTrackHitsID_CRL.end(), parentTrackHits.size(), leadingCount);
                    this->FillContaminationHitsDistance(parentTrackHits, leadingParticleHitList, bestMatchParentTrackHitsDistance_CRL);

                    bestMatchCRLHitsInCRID_CRL.insert(bestMatchCRLHitsInCRID_CRL.end(), leadingParticleHitsInParentCosmicRay.size(), leadingCount);
                    this->FillContaminationHitsDistance(leadingParticleHitsInParentCosmicRay, cosmicRayHitList, bestMatchCRLHitsInCRDistance_CRL);

                    const VertexList &deltaRayVertexList(pMatchedPfo->GetVertexList());

                    if (deltaRayVertexList.size() > 1)
                    {
                        std::cout << "ISOBEL: DELTA RAY RECO VERTEX LIST LARGER THAN 1" << std::endl;
                        throw;
                    }

                    bestMatchVertexX_CRL.push_back(deltaRayVertexList.empty() ? 1000001 : deltaRayVertexList.front()->GetPosition().GetX());
                    bestMatchVertexY_CRL.push_back(deltaRayVertexList.empty() ? 1000001 : deltaRayVertexList.front()->GetPosition().GetY());
                    bestMatchVertexZ_CRL.push_back(deltaRayVertexList.empty() ? 1000001 : deltaRayVertexList.front()->GetPosition().GetZ());                    
                }
                                      
                stringStream << "-" << (!isGoodMatch ? "(Below threshold) " : "")
                    << "nPfoHits " << pfoHitList.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")"
                    << ", nMatchedHits " << sharedHitList.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                    << ",  nCRLHitsInParentCR " << leadingParticleHitsInParentCosmicRay.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitsInParentCosmicRay)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitsInParentCosmicRay)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitsInParentCosmicRay) << ")" << std::endl                 
                    <<  (!isGoodMatch ? "                   " : " ")
                    << "nParentTrackHits " << parentTrackHits.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, parentTrackHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, parentTrackHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, parentTrackHits) << ")"
                    << ", nOtherTrackHits " << otherTrackHits.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, otherTrackHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, otherTrackHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, otherTrackHits) << ")"
                    << ", nOtherShowerHits " << otherShowerHits.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, otherShowerHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, otherShowerHits)
                             << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, otherShowerHits) << ")" << std::endl                    
                    <<  (!isGoodMatch ? "                   " : " ")                    
                    << (isMatchedToCorrectCosmicRay ? "Correct" :  "Incorrect") << "\033[0m" << " parent link" << std::endl;

                ///////////////////////////////
                if (m_visualize && useInterpretedMatching)
                {
                    std::cout << stringStream.str() << std::endl;
                    std::cout << "DELTA RAY PFO HITS" << std::endl;

                    const VertexList &deltaRayVertexList(pMatchedPfo->GetVertexList());

                    if (deltaRayVertexList.size() > 1)
                    {
                        std::cout << "ISOBEL: DELTA RAY RECO VERTEX LIST LARGER THAN 1" << std::endl;
                        throw;
                    }

                    if (deltaRayVertexList.empty())
                        std::cout << "No reconstructed vertex" << std::endl;
                    
                    this->PrintHits(pfoHitList, otherShowerHits, otherTrackHits, parentTrackHits, "DR_PFO");

                    if (!deltaRayVertexList.empty())
                    {

                        PfoList deltaRayList({pMatchedPfo});
                        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &deltaRayList, "Delta Ray Pfo", BLACK);
                        
                        const CartesianVector &vertex(deltaRayVertexList.front()->GetPosition());
                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vertex, "vertex", RED, 2);
                        PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    }

                    if (pParentPfo != pMatchedPfo)
                    {
                        std::cout << "PARENT PFO" << std::endl;
                        
                        const CaloHitList &parentCRHits(foldedPfoToHitsMap.at(pParentPfo));
                        this->PrintHits(parentCRHits, leadingParticleHitList, "DR_PARENT_PFO");
                    }
                }
                ///////////////////////////////
             }
                
            nAboveThresholdMatches_CRL.push_back(nAboveThresholdMatches);

            const bool isCorrect((nAboveThresholdMatches == 1) && isCorrectParentLink);
            
            if (isCorrect)
            {
                ++nCorrectChildCRLs;
                isCorrect_CRL.push_back(1);
            }
            else
            {
                isCorrect_CRL.push_back(0);
            }

            if (foldedMCToPfoHitSharingMap.at(pLeadingParticle).empty())
            {
                stringStream << "-" << "No matched pfo" << std::endl;

                if (m_visualize && useInterpretedMatching)
                {
                    std::cout << stringStream.str() << std::endl;
                }
            }

            if (nAboveThresholdMatches == 0)
            {
                isCorrectParentLink_CRL.push_back(0);                
                bestMatchNHitsTotal_CRL.push_back(0); bestMatchNHitsU_CRL.push_back(0); bestMatchNHitsV_CRL.push_back(0); bestMatchNHitsW_CRL.push_back(0);
                bestMatchNSharedHitsTotal_CRL.push_back(0); bestMatchNSharedHitsU_CRL.push_back(0); bestMatchNSharedHitsV_CRL.push_back(0); bestMatchNSharedHitsW_CRL.push_back(0);
                bestMatchNParentTrackHitsTotal_CRL.push_back(0); bestMatchNParentTrackHitsU_CRL.push_back(0); bestMatchNParentTrackHitsV_CRL.push_back(0); bestMatchNParentTrackHitsW_CRL.push_back(0);
                bestMatchNOtherTrackHitsTotal_CRL.push_back(0); bestMatchNOtherTrackHitsU_CRL.push_back(0), bestMatchNOtherTrackHitsV_CRL.push_back(0); bestMatchNOtherTrackHitsW_CRL.push_back(0);
                bestMatchNOtherShowerHitsTotal_CRL.push_back(0); bestMatchNOtherShowerHitsU_CRL.push_back(0); bestMatchNOtherShowerHitsV_CRL.push_back(0), bestMatchNOtherShowerHitsW_CRL.push_back(0);
                totalCRLHitsInBestMatchParentCR_CRL.push_back(0);
                uCRLHitsInBestMatchParentCR_CRL.push_back(0); vCRLHitsInBestMatchParentCR_CRL.push_back(0); wCRLHitsInBestMatchParentCR_CRL.push_back(0);
            }
            
            stringStream << nAboveThresholdMatches << " above threshold matches" << std::endl
                         << "Reconstruction is "
                         << (isCorrect ?  "\033[32m" : "\033[31m")
                         << (isCorrect ?  "CORRECT" : "INCORRECT") << "\033[0m" << std::endl;

            if (m_visualize && useInterpretedMatching)
                std::cout << stringStream.str() << std::endl;
        }

        if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ID_CR", ID_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcE_CR", mcE_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPX_CR", mcPX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPY_CR", mcPY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPZ_CR", mcPZ_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsTotal_CR", nMCHitsTotal_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsU_CR", nMCHitsU_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsV_CR", nMCHitsV_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsW_CR", nMCHitsW_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexX_CR", mcVertexX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexY_CR", mcVertexY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexZ_CR", mcVertexZ_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndX_CR", mcEndX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndY_CR", mcEndY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndZ_CR", mcEndZ_CR));            
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nReconstructableChildCRLs", nReconstructableChildCRLs));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nCorrectChildCRLs", nCorrectChildCRLs));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ID_CRL", &ID_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcE_CRL", &mcE_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPX_CRL", &mcPX_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPY_CRL", &mcPY_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPZ_CRL", &mcPZ_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsTotal_CRL", &nMCHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsU_CRL", &nMCHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsV_CRL", &nMCHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsW_CRL", &nMCHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexX_CRL", &mcVertexX_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexY_CRL", &mcVertexY_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexZ_CRL", &mcVertexZ_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndX_CRL", &mcEndX_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndY_CRL", &mcEndY_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndZ_CRL", &mcEndZ_CRL));             
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nAboveThresholdMatches_CRL", &nAboveThresholdMatches_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrect_CRL", &isCorrect_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectParentLink_CRL", &isCorrectParentLink_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsTotal_CRL", &bestMatchNHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsU_CRL", &bestMatchNHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsV_CRL", &bestMatchNHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsW_CRL", &bestMatchNHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsTotal_CRL", &bestMatchNSharedHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsU_CRL", &bestMatchNSharedHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsV_CRL", &bestMatchNSharedHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsW_CRL", &bestMatchNSharedHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNParentTrackHitsTotal_CRL", &bestMatchNParentTrackHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNParentTrackHitsU_CRL", &bestMatchNParentTrackHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNParentTrackHitsV_CRL", &bestMatchNParentTrackHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNParentTrackHitsW_CRL", &bestMatchNParentTrackHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherTrackHitsTotal_CRL", &bestMatchNOtherTrackHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherTrackHitsU_CRL", &bestMatchNOtherTrackHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherTrackHitsV_CRL", &bestMatchNOtherTrackHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherTrackHitsW_CRL", &bestMatchNOtherTrackHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherShowerHitsTotal_CRL", &bestMatchNOtherShowerHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherShowerHitsU_CRL", &bestMatchNOtherShowerHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherShowerHitsV_CRL", &bestMatchNOtherShowerHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherShowerHitsW_CRL", &bestMatchNOtherShowerHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "totalCRLHitsInBestMatchParentCR_CRL", &totalCRLHitsInBestMatchParentCR_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "uCRLHitsInBestMatchParentCR_CRL", &uCRLHitsInBestMatchParentCR_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vCRLHitsInBestMatchParentCR_CRL", &vCRLHitsInBestMatchParentCR_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "wCRLHitsInBestMatchParentCR_CRL", &wCRLHitsInBestMatchParentCR_CRL));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchOtherShowerHitsID_CRL", &bestMatchOtherShowerHitsID_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchOtherShowerHitsDistance_CRL", &bestMatchOtherShowerHitsDistance_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchOtherTrackHitsID_CRL", &bestMatchOtherTrackHitsID_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchOtherTrackHitsDistance_CRL", &bestMatchOtherTrackHitsDistance_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchParentTrackHitsID_CRL", &bestMatchParentTrackHitsID_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchParentTrackHitsDistance_CRL", &bestMatchParentTrackHitsDistance_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchCRLHitsInCRID_CRL", &bestMatchCRLHitsInCRID_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchCRLHitsInCRDistance_CRL", &bestMatchCRLHitsInCRDistance_CRL));

            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }

        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;      
        stringStream << nCorrectChildCRLs << " / " << nReconstructableChildCRLs << " CRLs correctly reconstructed" << std::endl;
        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;
        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;
    }
    
    if (printToScreen && !(m_visualize && useInterpretedMatching))
    {
        std::cout << stringStream.str() << std::endl;
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    //std::cout << "Muon Count: " << muonCount << std::endl;
    //std::cout << "Reconstructable CRL Count: " << totalReconstructableCRLs << std::endl;
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());    
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonLeadingEventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodHits", m_validationParameters.m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForGoodView", m_validationParameters.m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodViews", m_validationParameters.m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_validationParameters.m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_validationParameters.m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_validationParameters.m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToPrimaries", m_validationParameters.m_foldBackHierarchy));

    //////////////////////////////////////////////////////////

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DeltaRayMode", m_deltaRayMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MichelMode", m_michelMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonsToSkip", m_muonsToSkip));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Visualize", m_visualize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RemoveRecoMuonHits", m_removeRecoMuonHits));                                                

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IgnoreIncorrectMuons", m_ignoreIncorrectMuons));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteRawMatchesToTree", m_writeRawMatchesToTree));    
    
    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::GetRecoMuonHits(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, CaloHitList &recoMuonHitList) const
{
    recoMuonHitList.clear();
    
    ValidationInfo validationInfo;
    
    if (pMCParticleList && pCaloHitList)
    {
        // Get reconstructable MCParticle hit ownership map (non-muon leading hierarchy is folded whilst muon is unfolded)
        LArMuonLeadingHelper::ValidationParameters recoValidationParams(m_validationParameters);
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;        
        LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, recoValidationParams, recoMuonHitList,
            targetMCParticleToHitsMap, this->GetPandora());    
        
        LArMuonLeadingHelper::ValidationParameters allValidationParams(m_validationParameters);
        allValidationParams.m_minPrimaryGoodHits = 0;
        allValidationParams.m_minHitsForGoodView = 0;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, allValidationParams, recoMuonHitList,
            allMCParticleToHitsMap, this->GetPandora());
        
        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }

    if (pPfoList)
    {        
        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, false);

        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }
    
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);

    const LArMCParticleHelper::MCContributionMap &foldedAllMCToHitsMap(validationInfo.GetAllMCParticleToHitsMap());
    const LArMCParticleHelper::PfoContributionMap &foldedPfoToHitsMap(validationInfo.GetPfoToHitsMap());
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &foldedMCToPfoHitSharingMap(validationInfo.GetInterpretedMCToPfoHitSharingMap());
    
    // Consider only delta rays from reconstructable CR muons
    MCParticleVector mcCRVector;
    for (auto &entry : foldedAllMCToHitsMap)
    {
        if (LArMCParticleHelper::IsCosmicRay(entry.first))
            mcCRVector.push_back(entry.first);
    }
    std::sort(mcCRVector.begin(), mcCRVector.end(), LArMCParticleHelper::SortByMomentum);

   
    for (const MCParticle *const pCosmicRay : mcCRVector)
    {
        for (LArMCParticleHelper::PfoToSharedHitsVector::const_iterator cosmicRayMatchedPfoPair = foldedMCToPfoHitSharingMap.at(pCosmicRay).begin();
             cosmicRayMatchedPfoPair != foldedMCToPfoHitSharingMap.at(pCosmicRay).end(); ++cosmicRayMatchedPfoPair)
        {
            const ParticleFlowObject *const pCosmicRayPfo(cosmicRayMatchedPfoPair->first);
                    
            recoMuonHitList.insert(recoMuonHitList.end(), foldedPfoToHitsMap.at(pCosmicRayPfo).begin(), foldedPfoToHitsMap.at(pCosmicRayPfo).end());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::PrintHits(const CaloHitList caloHitList, const std::string &stringTag, const Color &colour) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());                        

        if (pCaloHit->GetHitType() == TPC_VIEW_U)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, colour, 2);
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());                
        
        if (pCaloHit->GetHitType() == TPC_VIEW_V)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, colour, 2);            
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());        

        if (pCaloHit->GetHitType() == TPC_VIEW_W)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, colour, 2);            
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::PrintHits(const CaloHitList totalCaloHitList, const CaloHitList otherShowerCaloHitList,
    const CaloHitList otherTrackCaloHitList, const CaloHitList parentTrackCaloHitList, const std::string &stringTag) const
{
    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_U)
            continue;

        Color color(BLACK);
        std::string newStringTag(stringTag);
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());                       
        
        if (std::find(otherShowerCaloHitList.begin(), otherShowerCaloHitList.end(), pCaloHit) != otherShowerCaloHitList.end())
            newStringTag += "_OTHER_SHOWER"; color = VIOLET;

        if (std::find(otherTrackCaloHitList.begin(), otherTrackCaloHitList.end(), pCaloHit) != otherTrackCaloHitList.end())
            newStringTag += "_OTHER_TRACK"; color = RED;

        if (std::find(parentTrackCaloHitList.begin(), parentTrackCaloHitList.end(), pCaloHit) != parentTrackCaloHitList.end())
            newStringTag += "_PARENT_TRACK"; color = BLUE;
        
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());  

    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_V)
            continue;

        Color color(BLACK);
        std::string newStringTag(stringTag);        
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());                   
            
        if (std::find(otherShowerCaloHitList.begin(), otherShowerCaloHitList.end(), pCaloHit) != otherShowerCaloHitList.end())
            newStringTag += "_OTHER_SHOWER"; color = VIOLET;

        if (std::find(otherTrackCaloHitList.begin(), otherTrackCaloHitList.end(), pCaloHit) != otherTrackCaloHitList.end())
            newStringTag += "_OTHER_TRACK"; color = RED;

        if (std::find(parentTrackCaloHitList.begin(), parentTrackCaloHitList.end(), pCaloHit) != parentTrackCaloHitList.end())
            newStringTag += "_PARENT_TRACK"; color = BLUE;

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());  

    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_W)
            continue;

        Color color(BLACK);
        std::string newStringTag(stringTag);        
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());           
            
        if (std::find(otherShowerCaloHitList.begin(), otherShowerCaloHitList.end(), pCaloHit) != otherShowerCaloHitList.end())
            newStringTag += "_OTHER_SHOWER"; color = VIOLET;

        if (std::find(otherTrackCaloHitList.begin(), otherTrackCaloHitList.end(), pCaloHit) != otherTrackCaloHitList.end())
            newStringTag += "_OTHER_TRACK"; color = RED;

        if (std::find(parentTrackCaloHitList.begin(), parentTrackCaloHitList.end(), pCaloHit) != parentTrackCaloHitList.end())
            newStringTag += "_PARENT_TRACK"; color = BLUE;

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    
    PandoraMonitoringApi::ViewEvent(this->GetPandora());  
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::PrintHits(const CaloHitList totalCaloHitList, const CaloHitList leadingCaloHitList,
    const std::string &stringTag) const
{
    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_U)
            continue;

        Color color(DARKGREEN);
        std::string newStringTag(stringTag);
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
        
        if (std::find(leadingCaloHitList.begin(), leadingCaloHitList.end(), pCaloHit) != leadingCaloHitList.end())
            newStringTag += "_LEADING"; color = RED;

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());  

    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_V)
            continue;

        Color color(DARKGREEN);
        std::string newStringTag(stringTag);        
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());        
            
        if (std::find(leadingCaloHitList.begin(), leadingCaloHitList.end(), pCaloHit) != leadingCaloHitList.end())
            newStringTag += "_LEADING"; color = RED;
        
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());  

    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_W)
            continue;

        Color color(DARKGREEN);
        std::string newStringTag(stringTag);        
        const CartesianVector hitPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());                
            
        if (std::find(leadingCaloHitList.begin(), leadingCaloHitList.end(), pCaloHit) != leadingCaloHitList.end())
            newStringTag += "_LEADING"; color = RED;

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());  
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::FillContaminationHitsDistance(const CaloHitList &contaminationHits, const CaloHitList &leadingMCHits,
    FloatVector &bestMatchContaminationHitsDistance) const
{
    CaloHitList leadingU, leadingV, leadingW;
    this->GetHitsOfType(leadingMCHits, TPC_VIEW_U, leadingU);
    this->GetHitsOfType(leadingMCHits, TPC_VIEW_V, leadingV);
    this->GetHitsOfType(leadingMCHits, TPC_VIEW_W, leadingW);    
    
    for (const CaloHit *const pContaminationHit : contaminationHits)
    {
        const CartesianVector &hitPosition(pContaminationHit->GetPositionVector());
        
        if (pContaminationHit->GetHitType() == TPC_VIEW_U)
            bestMatchContaminationHitsDistance.push_back(LArClusterHelper::GetClosestDistance(hitPosition, leadingU));

        if (pContaminationHit->GetHitType() == TPC_VIEW_V)
            bestMatchContaminationHitsDistance.push_back(LArClusterHelper::GetClosestDistance(hitPosition, leadingV));
            
        if (pContaminationHit->GetHitType() == TPC_VIEW_W)
            bestMatchContaminationHitsDistance.push_back(LArClusterHelper::GetClosestDistance(hitPosition, leadingW));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::GetHitsOfType(const CaloHitList &inputList, const HitType &hitType, CaloHitList &outputList) const
{
    for (const CaloHit *const pCaloHit : inputList)
    {
        if (pCaloHit->GetHitType() == hitType)
            outputList.push_back(pCaloHit);
    }   
}

} // namespace lar_content    
