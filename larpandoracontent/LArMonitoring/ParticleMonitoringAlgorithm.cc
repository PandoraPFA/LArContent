/**
 *  @file   LArContent/src/LArMonitoring/ParticleMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArContent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArContent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArContent/LArObjects/LArTrackPfo.h"
#include "larpandoracontent/LArContent/LArObjects/LArShowerPfo.h"

#include "larpandoracontent/LArContent/LArMonitoring/ParticleMonitoringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ParticleMonitoringAlgorithm::ParticleMonitoringAlgorithm() :
    m_primaryPfosOnly(true),
    m_collapseToPrimaryPfos(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleMonitoringAlgorithm::~ParticleMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleMonitoringAlgorithm::Run()
{
#ifdef MONITORING
    // Tree elements
    int nMCParticlesTotal(0), nPfosTotal(0);

    IntVector mcPdgVector, mcNuPdgVector, nMCHitsVector, pfoPdgVector, pfoNuPdgVector, nPfoHitsVector, nMatchedHitsVector,
        nMCHitsUVector, nPfoHitsUVector, nMatchedHitsUVector, nMCHitsVVector, nPfoHitsVVector, nMatchedHitsVVector, nMCHitsWVector,
        nPfoHitsWVector, nMatchedHitsWVector, nAvailableHitsVector, pfoVertexVector, pfoTrackVector, pfoShowerVector;

    FloatVector completenessVector, purityVector, mcPxVector, mcPyVector, mcPzVector, mcThetaVector, mcEnergyVector, mcPTotVector,
        mcVtxXVector, mcVtxYVector, mcVtxZVector, mcEndXVector, mcEndYVector, mcEndZVector,
        pfoPxVector, pfoPyVector, pfoPzVector, pfoPTotVector, pfoVtxXVector, pfoVtxYVector, pfoVtxZVector,
        trkVtxXVector, trkVtxYVector, trkVtxZVector, trkVtxDirXVector, trkVtxDirYVector, trkVtxDirZVector,
        trkEndXVector, trkEndYVector, trkEndZVector, trkEndDirXVector, trkEndDirYVector, trkEndDirZVector;

    try
    {
        // Load List of MC particles
        const MCParticleList *pMCParticleList = NULL;
        (void) PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);

        if (NULL == pMCParticleList)
        {
            std::cout << "ParticleMonitoringAlgorithm: cannot find mc particle list " << m_mcParticleListName << std::endl;
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }

        // Load List of Calo Hits
        const CaloHitList *pCaloHitList = NULL;
        (void) PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList);

        if (NULL == pCaloHitList)
        {
            std::cout << "ParticleMonitoringAlgorithm: cannot find calo hit list " << m_caloHitListName << std::endl;
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }

        // Load List of Pfos
        const PfoList *pPfoList = NULL;
        PfoList inputPfoList((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) ? PfoList(*pPfoList) : PfoList());

        PfoList pfoList;
        LArMonitoringHelper::ExtractTargetPfos(inputPfoList, m_primaryPfosOnly, pfoList);

        nPfosTotal = pfoList.size();

        PfoList recoNeutrinos;                                          // reco neutrinos
        LArPfoHelper::GetRecoNeutrinos(pPfoList, recoNeutrinos);

        MCParticleVector mcNeutrinoList;                                // true neutrinos
        LArMCParticleHelper::GetNeutrinoMCParticleList(pMCParticleList, mcNeutrinoList);

        MCParticleVector mcPrimaryList;                                 // primary mc particles
        LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryList);

        LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;            // [mc particles -> primary mc particle]
        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

        LArMonitoringHelper::MCContributionMap mcToTrueHitListMap;      // [primary mc particle -> true hit list]
        LArMonitoringHelper::CaloHitToMCMap hitToPrimaryMCMap;          // [hit -> primary mc particle]
        LArMonitoringHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToPrimaryMCMap, hitToPrimaryMCMap, mcToTrueHitListMap);

        LArMonitoringHelper::CaloHitToPfoMap hitToPfoMap;               // [hit -> pfo]
        LArMonitoringHelper::PfoContributionMap pfoToHitListMap;        // [pfo -> reco hit list]
        LArMonitoringHelper::GetPfoToCaloHitMatches(pCaloHitList, pfoList, m_collapseToPrimaryPfos, hitToPfoMap, pfoToHitListMap);

        LArMonitoringHelper::MCToPfoMap mcToBestPfoMap;                 // [mc particle -> best matched pfo]
        LArMonitoringHelper::MCContributionMap mcToBestPfoHitsMap;      // [mc particle -> list of hits included in best pfo]
        LArMonitoringHelper::MCToPfoMatchingMap mcToFullPfoMatchingMap; // [mc particle -> all matched pfos (and matched hits)]
        LArMonitoringHelper::GetMCParticleToPfoMatches(pCaloHitList, pfoToHitListMap, hitToPrimaryMCMap, mcToBestPfoMap, mcToBestPfoHitsMap, mcToFullPfoMatchingMap);

        LArMonitoringHelper::MCToPfoMap matchedNeutrinoMap;             // [neutrino mc particle -> neutrino pfo]
        LArMonitoringHelper::GetNeutrinoMatches(pCaloHitList, recoNeutrinos, hitToPrimaryMCMap, matchedNeutrinoMap);

        nMCParticlesTotal = mcToTrueHitListMap.size();

        // Write out reco/true neutrino information
        std::sort(mcNeutrinoList.begin(), mcNeutrinoList.end(), LArMCParticleHelper::SortByMomentum);

        for (MCParticleVector::const_iterator mcIter = mcNeutrinoList.begin(), mcIterEnd = mcNeutrinoList.end();
            mcIter != mcIterEnd; ++mcIter)
        {
            const MCParticle *const pTrueNeutrino(*mcIter);

            mcNuPdgVector.push_back(LArMCParticleHelper::GetParentNeutrinoId(pTrueNeutrino));
            mcPdgVector.push_back(pTrueNeutrino->GetParticleId());
            mcPxVector.push_back(pTrueNeutrino->GetMomentum().GetX());
            mcPyVector.push_back(pTrueNeutrino->GetMomentum().GetY());
            mcPzVector.push_back(pTrueNeutrino->GetMomentum().GetZ());
            mcPTotVector.push_back(pTrueNeutrino->GetMomentum().GetMagnitude());
            mcEnergyVector.push_back(pTrueNeutrino->GetEnergy());
            mcThetaVector.push_back(pTrueNeutrino->GetMomentum().GetOpeningAngle(CartesianVector(0.f, 0.f, 1.f)));

            mcVtxXVector.push_back(pTrueNeutrino->GetVertex().GetX());
            mcVtxYVector.push_back(pTrueNeutrino->GetVertex().GetY());
            mcVtxZVector.push_back(pTrueNeutrino->GetVertex().GetZ());
            mcEndXVector.push_back(pTrueNeutrino->GetEndpoint().GetX());
            mcEndYVector.push_back(pTrueNeutrino->GetEndpoint().GetY());
            mcEndZVector.push_back(pTrueNeutrino->GetEndpoint().GetZ());

            purityVector.push_back(0.f);
            completenessVector.push_back(0.f);
            nMCHitsVector.push_back(0);
            nPfoHitsVector.push_back(0);
            nMatchedHitsVector.push_back(0);
            nMCHitsUVector.push_back(0);
            nPfoHitsUVector.push_back(0);
            nMatchedHitsUVector.push_back(0);
            nMCHitsVVector.push_back(0);
            nPfoHitsVVector.push_back(0);
            nMatchedHitsVVector.push_back(0);
            nMCHitsWVector.push_back(0);
            nPfoHitsWVector.push_back(0);
            nMatchedHitsWVector.push_back(0);
            nAvailableHitsVector.push_back(0);

            int pfoVertex(0), pfoPdg(0), pfoNeutrinoPdg(0);
            float pfoPTot(0.f), pfoPx(0.f), pfoPy(0.f), pfoPz(0.f), pfoVtxX(0.f), pfoVtxY(0.f), pfoVtxZ(0.f);

            LArMonitoringHelper::MCToPfoMap::const_iterator pIter = matchedNeutrinoMap.find(pTrueNeutrino);
            if (matchedNeutrinoMap.end() != pIter)
            {
                const ParticleFlowObject *const pRecoNeutrino(pIter->second);

                pfoPdg = pRecoNeutrino->GetParticleId();
                pfoNeutrinoPdg = LArPfoHelper::GetPrimaryNeutrino(pRecoNeutrino);

                pfoPTot = pRecoNeutrino->GetMomentum().GetMagnitude();
                pfoPx = pRecoNeutrino->GetMomentum().GetX();
                pfoPy = pRecoNeutrino->GetMomentum().GetY();
                pfoPz = pRecoNeutrino->GetMomentum().GetZ();

                if (pRecoNeutrino->GetVertexList().size() > 0)
                {
                    if (pRecoNeutrino->GetVertexList().size() != 1)
                        throw StatusCodeException(STATUS_CODE_FAILURE);

                    const Vertex *const pVertex = *(pRecoNeutrino->GetVertexList().begin());

                    pfoVertex = 1;
                    pfoVtxX = pVertex->GetPosition().GetX();
                    pfoVtxY = pVertex->GetPosition().GetY();
                    pfoVtxZ = pVertex->GetPosition().GetZ();
                }
            }

            pfoPdgVector.push_back(pfoPdg);
            pfoNuPdgVector.push_back(pfoNeutrinoPdg);
            pfoPxVector.push_back(pfoPx);
            pfoPyVector.push_back(pfoPy);
            pfoPzVector.push_back(pfoPz);
            pfoPTotVector.push_back(pfoPTot);
            pfoVertexVector.push_back(pfoVertex);
            pfoVtxXVector.push_back(pfoVtxX);
            pfoVtxYVector.push_back(pfoVtxY);
            pfoVtxZVector.push_back(pfoVtxZ);

            pfoTrackVector.push_back(0);
            trkVtxXVector.push_back(0.f);
            trkVtxYVector.push_back(0.f);
            trkVtxZVector.push_back(0.f);
            trkVtxDirXVector.push_back(0.f);
            trkVtxDirYVector.push_back(0.f);
            trkVtxDirZVector.push_back(0.f);
            trkEndXVector.push_back(0.f);
            trkEndYVector.push_back(0.f);
            trkEndZVector.push_back(0.f);
            trkEndDirXVector.push_back(0.f);
            trkEndDirYVector.push_back(0.f);
            trkEndDirZVector.push_back(0.f);

            pfoShowerVector.push_back(0);
        }

        // Write out reco/true particle information
        std::sort(mcPrimaryList.begin(), mcPrimaryList.end(), LArMCParticleHelper::SortBySource);

        for (MCParticleVector::const_iterator iter1 = mcPrimaryList.begin(), iterEnd1 = mcPrimaryList.end(); iter1 != iterEnd1; ++iter1)
        {
            const MCParticle *const pMCParticle3D = *iter1;

            LArMonitoringHelper::MCContributionMap::const_iterator iter2 = mcToTrueHitListMap.find(pMCParticle3D);
            if (mcToTrueHitListMap.end() == iter2)
                continue;

            const CaloHitList &mcHitList(iter2->second);

            mcNuPdgVector.push_back(LArMCParticleHelper::GetParentNeutrinoId(pMCParticle3D));
            mcPdgVector.push_back(pMCParticle3D->GetParticleId());
            mcPxVector.push_back(pMCParticle3D->GetMomentum().GetX());
            mcPyVector.push_back(pMCParticle3D->GetMomentum().GetY());
            mcPzVector.push_back(pMCParticle3D->GetMomentum().GetZ());
            mcPTotVector.push_back(pMCParticle3D->GetMomentum().GetMagnitude());
            mcEnergyVector.push_back(pMCParticle3D->GetEnergy());
            mcThetaVector.push_back(pMCParticle3D->GetMomentum().GetOpeningAngle(CartesianVector(0.f, 0.f, 1.f)));

            mcVtxXVector.push_back(pMCParticle3D->GetVertex().GetX());
            mcVtxYVector.push_back(pMCParticle3D->GetVertex().GetY());
            mcVtxZVector.push_back(pMCParticle3D->GetVertex().GetZ());
            mcEndXVector.push_back(pMCParticle3D->GetEndpoint().GetX());
            mcEndYVector.push_back(pMCParticle3D->GetEndpoint().GetY());
            mcEndZVector.push_back(pMCParticle3D->GetEndpoint().GetZ());

            const int nMCHits(mcHitList.size());
            const int nMCHitsU(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcHitList));
            const int nMCHitsV(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcHitList));
            const int nMCHitsW(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcHitList));

            int pfoPdg(0), pfoNeutrinoPdg(0), nPfoHits(0), nPfoHitsU(0), nPfoHitsV(0), nPfoHitsW(0);
            int pfoVertex(0), pfoTrack(0), pfoShower(0);
            float pfoPTot(0.f), pfoPx(0.f), pfoPy(0.f), pfoPz(0.f), pfoVtxX(0.f), pfoVtxY(0.f), pfoVtxZ(0.f);
            float trkVtxX(0.f), trkVtxY(0.f), trkVtxZ(0.f), trkVtxDirX(0.f), trkVtxDirY(0.f), trkVtxDirZ(0.f);
            float trkEndX(0.f), trkEndY(0.f), trkEndZ(0.f), trkEndDirX(0.f), trkEndDirY(0.f), trkEndDirZ(0.f);

            LArMonitoringHelper::MCToPfoMap::const_iterator pIter1 = mcToBestPfoMap.find(pMCParticle3D);

            if (mcToBestPfoMap.end() != pIter1)
            {
                const ParticleFlowObject *const pPfo(pIter1->second);
                LArMonitoringHelper::PfoContributionMap::const_iterator pIter2 = pfoToHitListMap.find(pPfo);

                if (pfoToHitListMap.end() == pIter2)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                pfoPdg = pPfo->GetParticleId();
                pfoNeutrinoPdg = LArPfoHelper::GetPrimaryNeutrino(pPfo);

                pfoPTot = pPfo->GetMomentum().GetMagnitude();
                pfoPx = pPfo->GetMomentum().GetX();
                pfoPy = pPfo->GetMomentum().GetY();
                pfoPz = pPfo->GetMomentum().GetZ();

                if (pPfo->GetVertexList().size() > 0)
                {
                    if (pPfo->GetVertexList().size() != 1)
                        throw StatusCodeException(STATUS_CODE_FAILURE);

                    const Vertex *const pVertex = *(pPfo->GetVertexList().begin());

                    pfoVertex = 1;
                    pfoVtxX = pVertex->GetPosition().GetX();
                    pfoVtxY = pVertex->GetPosition().GetY();
                    pfoVtxZ = pVertex->GetPosition().GetZ();
                }

                const LArTrackPfo *const pLArTrackPfo = dynamic_cast<const LArTrackPfo*>(pPfo);
                if (NULL != pLArTrackPfo)
                {
                    pfoTrack = 1;
                    trkVtxX = pLArTrackPfo->GetVertexPosition().GetX();
                    trkVtxY = pLArTrackPfo->GetVertexPosition().GetY();
                    trkVtxZ = pLArTrackPfo->GetVertexPosition().GetZ();
                    trkVtxDirX = pLArTrackPfo->GetVertexDirection().GetX();
                    trkVtxDirY = pLArTrackPfo->GetVertexDirection().GetY();
                    trkVtxDirZ = pLArTrackPfo->GetVertexDirection().GetZ();
                    trkEndX = pLArTrackPfo->GetEndPosition().GetX();
                    trkEndY = pLArTrackPfo->GetEndPosition().GetY();
                    trkEndZ = pLArTrackPfo->GetEndPosition().GetZ();
                    trkEndDirX = pLArTrackPfo->GetEndDirection().GetX();
                    trkEndDirY = pLArTrackPfo->GetEndDirection().GetY();
                    trkEndDirZ = pLArTrackPfo->GetEndDirection().GetZ();
                }

                const LArShowerPfo *const pLArShowerPfo = dynamic_cast<const LArShowerPfo*>(pPfo);
                if (NULL != pLArShowerPfo)
                {
                    pfoShower = 1;
                }

                const CaloHitList &pfoHitList(pIter2->second);
                nPfoHits = pfoHitList.size();
                nPfoHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList);
                nPfoHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList);
                nPfoHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList);
            }

            int nMatchedHits(0), nMatchedHitsU(0), nMatchedHitsV(0), nMatchedHitsW(0);
            LArMonitoringHelper::MCContributionMap::const_iterator mIter = mcToBestPfoHitsMap.find(pMCParticle3D);

            if (mcToBestPfoHitsMap.end() != mIter)
            {
                const CaloHitList &matchedHitList(mIter->second);
                nMatchedHits = matchedHitList.size();
                nMatchedHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, matchedHitList);
                nMatchedHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, matchedHitList);
                nMatchedHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, matchedHitList);
            }

            int nRecoHits(0), nAvailableHits(nMCHits);
            LArMonitoringHelper::MCToPfoMatchingMap::const_iterator fIter = mcToFullPfoMatchingMap.find(pMCParticle3D);

            if (mcToFullPfoMatchingMap.end() != fIter)
            {
                for (const LArMonitoringHelper::PfoContributionMap::value_type contribution : fIter->second)
                    nRecoHits += contribution.second.size();

                nAvailableHits = nMCHits - nRecoHits;
            }

            pfoPdgVector.push_back(pfoPdg);
            pfoNuPdgVector.push_back(pfoNeutrinoPdg);
            pfoPxVector.push_back(pfoPx);
            pfoPyVector.push_back(pfoPy);
            pfoPzVector.push_back(pfoPz);
            pfoPTotVector.push_back(pfoPTot);
            pfoVtxXVector.push_back(pfoVtxX);
            pfoVtxYVector.push_back(pfoVtxY);
            pfoVtxZVector.push_back(pfoVtxZ);

            pfoVertexVector.push_back(pfoVertex);
            pfoTrackVector.push_back(pfoTrack);
            pfoShowerVector.push_back(pfoShower);

            trkVtxXVector.push_back(trkVtxX);
            trkVtxYVector.push_back(trkVtxY);
            trkVtxZVector.push_back(trkVtxZ);
            trkVtxDirXVector.push_back(trkVtxDirX);
            trkVtxDirYVector.push_back(trkVtxDirY);
            trkVtxDirZVector.push_back(trkVtxDirZ);
            trkEndXVector.push_back(trkEndX);
            trkEndYVector.push_back(trkEndY);
            trkEndZVector.push_back(trkEndZ);
            trkEndDirXVector.push_back(trkEndDirX);
            trkEndDirYVector.push_back(trkEndDirY);
            trkEndDirZVector.push_back(trkEndDirZ);

            purityVector.push_back((nPfoHits == 0) ? 0 : static_cast<float>(nMatchedHits) / static_cast<float>(nPfoHits));
            completenessVector.push_back((nPfoHits == 0) ? 0 : static_cast<float>(nMatchedHits) / static_cast<float>(nMCHits));
            nMCHitsVector.push_back(nMCHits);
            nPfoHitsVector.push_back(nPfoHits);
            nMatchedHitsVector.push_back(nMatchedHits);

            nMCHitsUVector.push_back(nMCHitsU);
            nPfoHitsUVector.push_back(nPfoHitsU);
            nMatchedHitsUVector.push_back(nMatchedHitsU);

            nMCHitsVVector.push_back(nMCHitsV);
            nPfoHitsVVector.push_back(nPfoHitsV);
            nMatchedHitsVVector.push_back(nMatchedHitsV);

            nMCHitsWVector.push_back(nMCHitsW);
            nPfoHitsWVector.push_back(nPfoHitsW);
            nMatchedHitsWVector.push_back(nMatchedHitsW);

            nAvailableHitsVector.push_back(nAvailableHits);
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_NOT_FOUND != statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCParticlesTotal", nMCParticlesTotal));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfosTotal", nPfosTotal));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPdg", &mcPdgVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNuPdg", &mcNuPdgVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPx", &mcPxVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPy", &mcPyVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPz", &mcPzVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPTot", &mcPTotVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEnergy", &mcEnergyVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcTheta", &mcThetaVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVtxX", &mcVtxXVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVtxY", &mcVtxYVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVtxZ", &mcVtxZVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndX", &mcEndXVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndY", &mcEndYVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndZ", &mcEndZVector));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoPdg", &pfoPdgVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoNuPdg", &pfoNuPdgVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoPx", &pfoPxVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoPy", &pfoPyVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoPz", &pfoPzVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoPTot", &pfoPTotVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoVtxX", &pfoVtxXVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoVtxY", &pfoVtxYVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoVtxZ", &pfoVtxZVector));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoVertex", &pfoVertexVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoTrack", &pfoTrackVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoShower", &pfoShowerVector));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkVtxX", &trkVtxXVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkVtxY", &trkVtxYVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkVtxZ", &trkVtxZVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkVtxDirX", &trkVtxDirXVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkVtxDirY", &trkVtxDirYVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkVtxDirZ", &trkVtxDirZVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkEndX", &trkEndXVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkEndY", &trkEndYVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkEndZ", &trkEndZVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkEndDirX", &trkEndDirXVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkEndDirY", &trkEndDirYVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkEndDirZ", &trkEndDirZVector));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "completeness", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "purity", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHits", &nMCHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHits", &nPfoHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMatchedHits", &nMatchedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nAvailableHits", &nAvailableHitsVector));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsU", &nMCHitsUVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsU", &nPfoHitsUVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMatchedHitsU", &nMatchedHitsUVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsV", &nMCHitsVVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsV", &nPfoHitsVVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMatchedHitsV", &nMatchedHitsVVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsW", &nMCHitsWVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsW", &nPfoHitsWVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMatchedHitsW", &nMatchedHitsWVector));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
#endif
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrimaryPfosOnly", m_primaryPfosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CollapseToPrimaryPfos", m_collapseToPrimaryPfos));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
