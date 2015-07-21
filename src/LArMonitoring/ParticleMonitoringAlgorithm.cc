/**
 *  @file   LArContent/src/LArMonitoring/ParticleMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArObjects/LArTrackPfo.h"
#include "LArObjects/LArShowerPfo.h"

#include "LArMonitoring/ParticleMonitoringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ParticleMonitoringAlgorithm::ParticleMonitoringAlgorithm() :
    m_useDaughterPfos(false),
    m_extractNeutrinoDaughters(true)
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
        PfoList pfoList((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) ? PfoList(*pPfoList) : PfoList());

        if (m_extractNeutrinoDaughters)
            this->ExtractNeutrinoDaughters(pfoList);

        nPfosTotal = pfoList.size();

        // Extract Reco Neutrinos
        PfoList recoNeutrinos;
        this->GetRecoNeutrinos(pPfoList, recoNeutrinos);

        // Extract True Neutrinos
        MCParticleVector mcNeutrinoList;
        LArMCParticleHelper::GetNeutrinoMCParticleList(pMCParticleList, mcNeutrinoList);

        // Extract Primary MC Particles
        MCParticleVector mcPrimaryList;
        LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryList);

        // Match Pfos to Primary MC Particles
        LArMCParticleHelper::MCRelationMap mcPrimaryMap;    // [particles -> primary particle]
        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

        MCContributionMap trueHitMap;    // [primary particle -> true hit list]
        CaloHitToMCMap mcHitMap;         // [hit -> parent primary]
        this->GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, mcHitMap, trueHitMap);

        PfoContributionMap recoHitMap;   // [pfo -> reco hit list]
        CaloHitToPfoMap pfoHitMap;       // [hit -> parent pfo]
        this->GetPfoToCaloHitMatches(pCaloHitList, pfoList, pfoHitMap, recoHitMap);

        MCContributionMap fullHitMap;    // [particle -> list of matched hits with any pfo]
        MCContributionMap matchedHitMap; // [particle -> list of matched hits with best pfo]
        MCToPfoMap matchedPfoMap;        // [particle -> best pfo]
        this->GetMCParticleToPfoMatches(pCaloHitList, pfoList, mcHitMap, matchedPfoMap, matchedHitMap, fullHitMap);

        MCToPfoMap matchedNeutrinoMap;   // [neutrino particle -> neutrino pfo]
        this->GetNeutrinoMatches(pCaloHitList, recoNeutrinos, mcHitMap, matchedNeutrinoMap);

        nMCParticlesTotal = trueHitMap.size();

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

            MCToPfoMap::const_iterator pIter = matchedNeutrinoMap.find(pTrueNeutrino);
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

            MCContributionMap::const_iterator iter2 = trueHitMap.find(pMCParticle3D);
            if (trueHitMap.end() == iter2)
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
            const int nMCHitsU(this->CountHitsByType(TPC_VIEW_U, mcHitList));
            const int nMCHitsV(this->CountHitsByType(TPC_VIEW_V, mcHitList));
            const int nMCHitsW(this->CountHitsByType(TPC_VIEW_W, mcHitList));

            int pfoPdg(0), pfoNeutrinoPdg(0), nPfoHits(0), nPfoHitsU(0), nPfoHitsV(0), nPfoHitsW(0);
            int pfoVertex(0), pfoTrack(0), pfoShower(0);
            float pfoPTot(0.f), pfoPx(0.f), pfoPy(0.f), pfoPz(0.f), pfoVtxX(0.f), pfoVtxY(0.f), pfoVtxZ(0.f);
            float trkVtxX(0.f), trkVtxY(0.f), trkVtxZ(0.f), trkVtxDirX(0.f), trkVtxDirY(0.f), trkVtxDirZ(0.f);
            float trkEndX(0.f), trkEndY(0.f), trkEndZ(0.f), trkEndDirX(0.f), trkEndDirY(0.f), trkEndDirZ(0.f);

            MCToPfoMap::const_iterator pIter1 = matchedPfoMap.find(pMCParticle3D);

            if (matchedPfoMap.end() != pIter1)
            {
                const ParticleFlowObject *const pPfo(pIter1->second);
                PfoContributionMap::const_iterator pIter2 = recoHitMap.find(pPfo);

                if (recoHitMap.end() == pIter2)
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
                nPfoHitsU = this->CountHitsByType(TPC_VIEW_U, pfoHitList);
                nPfoHitsV = this->CountHitsByType(TPC_VIEW_V, pfoHitList);
                nPfoHitsW = this->CountHitsByType(TPC_VIEW_W, pfoHitList);
            }

            int nMatchedHits(0), nMatchedHitsU(0), nMatchedHitsV(0), nMatchedHitsW(0);
            MCContributionMap::const_iterator mIter = matchedHitMap.find(pMCParticle3D);

            if (matchedHitMap.end() != mIter)
            {
                const CaloHitList &matchedHitList(mIter->second);
                nMatchedHits = matchedHitList.size();
                nMatchedHitsU = this->CountHitsByType(TPC_VIEW_U, matchedHitList);
                nMatchedHitsV = this->CountHitsByType(TPC_VIEW_V, matchedHitList);
                nMatchedHitsW = this->CountHitsByType(TPC_VIEW_W, matchedHitList);
            }

            int nRecoHits(0), nAvailableHits(nMCHits);
            MCContributionMap::const_iterator fIter = fullHitMap.find(pMCParticle3D);

            if (fullHitMap.end() != fIter)
            {
                const CaloHitList &fullHitList(fIter->second);
                nRecoHits = fullHitList.size();
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

void ParticleMonitoringAlgorithm::GetTrueNeutrinos(const MCParticleList *pMCParticleList, MCParticleList &trueNeutrinos) const
{
    if (NULL == pMCParticleList)
        return;

    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pMCParticle = *iter;

        if (pMCParticle->GetParentList().empty() && LArMCParticleHelper::IsNeutrino(pMCParticle))
            trueNeutrinos.insert(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetRecoNeutrinos(const PfoList *pPfoList, PfoList &recoNeutrinos) const
{
    if (NULL == pPfoList)
        return;

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pPfo(*iter);

        if (LArPfoHelper::IsNeutrino(pPfo))
            recoNeutrinos.insert(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::ExtractNeutrinoDaughters(PfoList &pfoList) const
{
    for (PfoList::iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; )
    {
        const ParticleFlowObject *const pPfo(*iter);

        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            pfoList.insert(pPfo->GetDaughterPfoList().begin(), pPfo->GetDaughterPfoList().end());
            pfoList.erase(iter);
            iter = pfoList.begin();
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetNeutrinoMatches(const CaloHitList *const pCaloHitList, const PfoList &recoNeutrinos,
    const CaloHitToMCMap &mcHitMap, MCToPfoMap &outputPrimaryMap) const
{
    for (PfoList::const_iterator pIter = recoNeutrinos.begin(), pIterEnd = recoNeutrinos.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pNeutrinoPfo = *pIter;

        if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PfoList pfoList;
        LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, pfoList);

        CaloHitList clusterHits;
        this->CollectCaloHits(pfoList, clusterHits);

        MCContributionMap inputContributionMap, outputContributionMap;

        for (CaloHitList::const_iterator hIter = clusterHits.begin(), hIterEnd = clusterHits.end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;

            if (pCaloHitList->count(pCaloHit) == 0)
                continue;

            CaloHitToMCMap::const_iterator mcIter = mcHitMap.find(pCaloHit);

            if (mcIter == mcHitMap.end())
                continue;

            const MCParticle *const pFinalStateParticle = mcIter->second;
            const MCParticle *const pNeutrinoParticle = LArMCParticleHelper::GetParentMCParticle(pFinalStateParticle);

            if (!LArMCParticleHelper::IsNeutrino(pNeutrinoParticle))
                continue;

            inputContributionMap[pNeutrinoParticle].insert(pCaloHit);
        }

        CaloHitList biggestHitList;
        const MCParticle *biggestContributor = NULL;

        for (MCContributionMap::const_iterator cIter = inputContributionMap.begin(), cIterEnd = inputContributionMap.end(); cIter != cIterEnd; ++cIter)
        {
            if (cIter->second.size() > biggestHitList.size())
            {
                biggestHitList = cIter->second;
                biggestContributor = cIter->first;
            }
        }

        if (biggestContributor)
        {
            MCContributionMap::const_iterator cIter = outputContributionMap.find(biggestContributor);

            if ((outputContributionMap.end() == cIter) || (biggestHitList.size() > cIter->second.size()))
            {
                outputPrimaryMap[biggestContributor] = pNeutrinoPfo;
                outputContributionMap[biggestContributor] = biggestHitList;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetMCParticleToPfoMatches(const CaloHitList *const pCaloHitList, const PfoList &pfoList,
    const CaloHitToMCMap &mcHitMap, MCToPfoMap &matchedPrimaryMap, MCContributionMap &matchedContributionMap,
    MCContributionMap &fullContributionMap) const
{
    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pPfo = *pIter;

        CaloHitList clusterHits;
        this->CollectCaloHits(pPfo, clusterHits);

        MCContributionMap mcContributionMap;
        for (CaloHitList::const_iterator hIter = clusterHits.begin(), hIterEnd = clusterHits.end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;

            if (TPC_3D == pCaloHit->GetHitType())
                throw StatusCodeException(STATUS_CODE_FAILURE);

            if (pCaloHitList->count(pCaloHit) == 0)
                continue;

            CaloHitToMCMap::const_iterator mcIter = mcHitMap.find(pCaloHit);

            if (mcIter == mcHitMap.end())
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;
            mcContributionMap[pPrimaryParticle].insert(pCaloHit);
            fullContributionMap[pPrimaryParticle].insert(pCaloHit);
        }

        CaloHitList biggestHitList;
        const MCParticle *biggestContributor = NULL;

        for (MCContributionMap::const_iterator cIter = mcContributionMap.begin(), cIterEnd = mcContributionMap.end(); cIter != cIterEnd; ++cIter)
        {
            if (cIter->second.size() > biggestHitList.size())
            {
                biggestHitList = cIter->second;
                biggestContributor = cIter->first;
            }
        }

        if (biggestContributor)
        {
            MCContributionMap::const_iterator cIter = matchedContributionMap.find(biggestContributor);

            if ((matchedContributionMap.end() == cIter) || (biggestHitList.size() > cIter->second.size()))
            {
                matchedPrimaryMap[biggestContributor] = pPfo;
                matchedContributionMap[biggestContributor] = biggestHitList;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetMCParticleToCaloHitMatches(const CaloHitList *const pCaloHitList, const MCRelationMap &mcPrimaryMap,
    CaloHitToMCMap &mcHitMap, MCContributionMap &mcContributionMap) const
{
    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            const CaloHit *const pCaloHit = *iter;
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            MCRelationMap::const_iterator mcIter = mcPrimaryMap.find(pHitParticle);

            if (mcPrimaryMap.end() == mcIter)
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;
            mcContributionMap[pPrimaryParticle].insert(pCaloHit);
            mcHitMap[pCaloHit] = pPrimaryParticle;
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetPfoToCaloHitMatches(const CaloHitList *const pCaloHitList, const PfoList &pfoList,
    CaloHitToPfoMap &pfoHitMap, PfoContributionMap &pfoContributionMap) const
{
    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pPfo = *pIter;

        ClusterList clusterList;
        LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);

        if (m_useDaughterPfos)
            this->CollectDaughterClusters(pPfo, clusterList);

        CaloHitList pfoHitList;

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pCluster = *cIter;

            CaloHitList clusterHits;
            pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterHits);
            clusterHits.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

            for (CaloHitList::const_iterator hIter = clusterHits.begin(), hIterEnd = clusterHits.end(); hIter != hIterEnd; ++hIter)
            {
                const CaloHit *const pCaloHit = *hIter;

                if (TPC_3D == pCaloHit->GetHitType())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                if (pCaloHitList->count(pCaloHit) == 0)
                    continue;

                pfoHitMap[pCaloHit] = pPfo;
                pfoHitList.insert(pCaloHit);
            }
        }

        pfoContributionMap[pPfo] = pfoHitList;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::CollectCaloHits(const ParticleFlowObject *const pParentPfo, CaloHitList &caloHitList) const
{
    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pParentPfo, clusterList);

    if (m_useDaughterPfos)
        this->CollectDaughterClusters(pParentPfo, clusterList);

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        caloHitList.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::CollectCaloHits(const PfoList &pfoList, CaloHitList &caloHitList) const
{
    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pPfo = *pIter;
        this->CollectCaloHits(pPfo, caloHitList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::CollectDaughterClusters(const ParticleFlowObject *const pParentPfo, ClusterList &clusterList) const
{
    const PfoList &daughterPfoList(pParentPfo->GetDaughterPfoList());
    for (PfoList::const_iterator dIter = daughterPfoList.begin(), dIterEnd = daughterPfoList.end(); dIter != dIterEnd; ++dIter)
    {
        const ParticleFlowObject *const pDaughterPfo = *dIter;

        for (ClusterList::const_iterator cIter = pDaughterPfo->GetClusterList().begin(), cIterEnd = pDaughterPfo->GetClusterList().end();
            cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pCluster = *cIter;

            if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
                continue;

            if (clusterList.count(pCluster))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            clusterList.insert(pCluster);
        }

        this->CollectDaughterClusters(pDaughterPfo, clusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int ParticleMonitoringAlgorithm::CountHitsByType(const HitType hitType, const CaloHitList &caloHitList) const
{
    unsigned int nHitsOfSpecifiedType(0);

    for (CaloHitList::const_iterator iter = caloHitList.begin(), iterEnd = caloHitList.end(); iter != iterEnd; ++iter)
    {
        if (hitType == (*iter)->GetHitType())
            ++nHitsOfSpecifiedType;
    }

    return nHitsOfSpecifiedType;
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
        "UseDaughterPfos", m_useDaughterPfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtractNeutrinoDaughters", m_extractNeutrinoDaughters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
