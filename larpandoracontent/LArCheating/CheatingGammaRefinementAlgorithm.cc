/**
 *  @file   larpandoracontent/LArCheating/CheatingGammaRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArCheating/CheatingGammaRefinementAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingGammaRefinementAlgorithm::CheatingGammaRefinementAlgorithm() :
    m_maxImpactT(9.f),
    m_minOpeningAngle(20.f),
    m_minGammaCompleteness(0.33f),
    m_removeHierarchy(false),
    m_truncateMode(false),
    m_creationCompletenessThreshold(0.6f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingGammaRefinementAlgorithm::Run()
{
    std::exception_ptr eptr;

    try
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

        if (!pCaloHitList || pCaloHitList->empty())
        {
            std::cout << "CheatingCCElectronRefinementAlgorithm: No calo hit list found, returning..." << std::endl;
            return STATUS_CODE_FAILURE;
        }

        const PfoList *pShowerPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_showerPfoListName, pShowerPfoList));

        if (!pShowerPfoList || pShowerPfoList->empty())
        {
            std::cout << "CheatingGammaRefinementAlgorithm: No shower pfo list found, returning..." << std::endl;
            return STATUS_CODE_FAILURE;
        }

        PfoVector showerVector(pShowerPfoList->begin(), pShowerPfoList->end());
        std::sort(showerVector.begin(), showerVector.end(), LArPfoHelper::SortByNHits);

        LArMCParticleHelper::MCContributionMap gammaHitMap;
        this->FillGammaHitMap(pCaloHitList, gammaHitMap);

        if (gammaHitMap.empty())
            return STATUS_CODE_SUCCESS;

        for (const ParticleFlowObject *const pPfo : showerVector)
        {
            if (LArPfoHelper::GetNumberOfTwoDHits(pPfo) < 100)
                continue;

            const MCParticle *const pMCGamma(this->FindMatchedGamma(pPfo, gammaHitMap));

            if (!pMCGamma)
                continue;

            //PfoList pfoDisplay;
            //pfoDisplay.push_back(pPfo);
            //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoDisplay, "matched gamma pfo - before", RED, true, true);        
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
            MCParticleList contaminantMCParticles;
            this->FindContaminantMCParticles(pMCGamma, pPfo, contaminantMCParticles);

            if (contaminantMCParticles.empty())
                continue;

            this->FilterContaminantMCParticleList(pMCGamma, pPfo, contaminantMCParticles);

            if (contaminantMCParticles.empty())
                continue;

            this->RemoveContaminantHits(pMCGamma, pPfo, contaminantMCParticles, gammaHitMap);

            //PfoList pfoDisplay2;
            //pfoDisplay2.push_back(pPfo);
            //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoDisplay2, "matched gamma pfo - after", RED, true, true);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }
    }
    catch (...)
    {
        eptr = std::current_exception();
        this->handle_eptr(eptr);
    }

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::FillGammaHitMap(const CaloHitList *const pCaloHitList, LArMCParticleHelper::MCContributionMap &gammaHitMap) const
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            
            //if (pMCParticle->GetParticleId() == 22)
            gammaHitMap[pMCParticle].push_back(pCaloHit);
        }
        catch (...)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *CheatingGammaRefinementAlgorithm::FindMatchedGamma(const ParticleFlowObject *const pPfo, const LArMCParticleHelper::MCContributionMap &gammaHitMap) const
{
    const MCParticle *pMainMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
    const int pdg(std::abs(pMainMCParticle->GetParticleId()));

    
    if (pdg == 22)
        return pMainMCParticle;
       
    if (pdg == 11)
        return nullptr;

    // but does it contain a large chunk of a gamma even if it is not the main owner?
    // find gamma with highest number of shared hits that has a completeness of >50 % in one view

    MCParticleVector mcGammaVector;

    for (auto &entry : gammaHitMap)
    {
        if (entry.first->GetParticleId() == 22)
            mcGammaVector.push_back(entry.first);
    }

    std::sort(mcGammaVector.begin(), mcGammaVector.end(), LArMCParticleHelper::SortByMomentum);

    CaloHitList pfoHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, pfoHitList);

    unsigned int highestSharedHits(0);
    const MCParticle *pBestGamma(nullptr);

    for (const MCParticle *const pMCGamma : mcGammaVector)
    {
        const CaloHitList &mcGammaHitList(gammaHitMap.at(pMCGamma));

        if (mcGammaHitList.size() < 100)
            continue;

        const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(pfoHitList, mcGammaHitList));

        if (sharedHitList.size() < highestSharedHits)
            continue;

        // get each view completeness
        for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            const float completeness(static_cast<float>(LArMonitoringHelper::CountHitsByType(hitType, sharedHitList)) / 
                static_cast<float>(LArMonitoringHelper::CountHitsByType(hitType, mcGammaHitList)));

            if (completeness > m_minGammaCompleteness)
            {
                highestSharedHits = sharedHitList.size();
                pBestGamma = pMCGamma;
            }
        }
    }

    return pBestGamma;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::FindContaminantMCParticles(const MCParticle *const pMCGamma, const ParticleFlowObject *const pGammaPfo, 
    MCParticleList &contaminantMCParticleList) const
{
    const CartesianVector &mcVertex(pMCGamma->GetVertex());
    const CartesianVector mcVertexU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_U));
    const CartesianVector mcVertexV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_V));
    const CartesianVector mcVertexW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_W));

    CartesianVector mcDirection(pMCGamma->GetMomentum());

    if (mcDirection.GetMagnitude() < std::numeric_limits<float>::epsilon())
        return;

    mcDirection = mcDirection * (1.f / mcDirection.GetMagnitude());
    const CartesianVector mcDirectionU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_U));
    const CartesianVector mcDirectionV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_V));
    const CartesianVector mcDirectionW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_W));

    ClusterList twoDClusterList;
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_U, twoDClusterList);
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_V, twoDClusterList);
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_W, twoDClusterList);

    //CartesianPointVector hitPositions;
    for (const Cluster *const pCluster : twoDClusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        const CaloHitList &isolatedHitList(pCluster->GetIsolatedCaloHitList());

        for (const CaloHit *const pIsolated : isolatedHitList)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                caloHitList.push_back(pIsolated);
        }

        const HitType hitType(caloHitList.front()->GetHitType());
        const CartesianVector viewMCVertex(hitType == TPC_VIEW_U ? mcVertexU : hitType == TPC_VIEW_V ? mcVertexV : mcVertexW);
        const CartesianVector viewMCDirection(hitType == TPC_VIEW_U ? mcDirectionU : hitType == TPC_VIEW_V ? mcDirectionV : mcDirectionW);

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            try
            {
                const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

                if (pMCParticle == pMCGamma)
                    continue;

                const CartesianVector &position(pCaloHit->GetPositionVector());
                const CartesianVector displacement(position - viewMCVertex);
                const float impactL(viewMCDirection.GetDotProduct(displacement));

                if (impactL > 0.f)
                    continue;

                //hitPositions.push_back(pCaloHit->GetPositionVector());

                if (std::find(contaminantMCParticleList.begin(), contaminantMCParticleList.end(), pMCParticle) != contaminantMCParticleList.end())
                    continue;

                contaminantMCParticleList.push_back(pMCParticle);
            }
            catch (const StatusCodeException &)
            {
            }
        }
    }

    
    //for (const CartesianVector jam : hitPositions)
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &jam, "remove", VIOLET, 1);

    
    //if (hitPositions.size())
    //{
    //PfoList pfoDisplay;
    //pfoDisplay.push_back(pGammaPfo);
    //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoDisplay, "matched gamma pfo", RED, true, true);

    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexU, "mcVertex", BLACK, 1);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexV, "mcVertex", BLACK, 1);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexW, "mcVertex", BLACK, 1);
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    //}
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::FilterContaminantMCParticleList(const MCParticle *const pMCGamma, const ParticleFlowObject *const pGammaPfo, 
   MCParticleList &mcContaminantList) const
{
    MCParticleVector mcContaminantVector(mcContaminantList.begin(), mcContaminantList.end());
    std::sort(mcContaminantVector.begin(), mcContaminantVector.end(), LArMCParticleHelper::SortByMomentum);

    mcContaminantList.clear();

    for (const MCParticle *const pMCContaminant : mcContaminantVector)
    {
        if (!this->AreParticlesSeparated(pMCGamma, pMCContaminant))
            continue;        

        if (std::find(mcContaminantList.begin(), mcContaminantList.end(), pMCContaminant) == mcContaminantList.end())
            mcContaminantList.push_back(pMCContaminant);

        if (m_removeHierarchy)
            this->GetDistinguishableChildren(pMCContaminant, pMCGamma, mcContaminantList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingGammaRefinementAlgorithm::IsShower(const MCParticle *const pMCParticle) const
{
  const int pdg(pMCParticle->GetParticleId());

  return ((E_MINUS == std::abs(pdg)) || (PHOTON == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingGammaRefinementAlgorithm::AreParticlesSeparated(const MCParticle *const pMCParticle1, const MCParticle *const pMCParticle2) const
{
    // Make sure particles are not inside a shower

    const CartesianVector &mcVertex1(pMCParticle1->GetVertex()), &mcVertex2(pMCParticle2->GetVertex());
    CartesianVector mcDirection1(pMCParticle1->GetMomentum()), mcDirection2(pMCParticle2->GetMomentum());

    if ((mcDirection1.GetMagnitude() < std::numeric_limits<float>::epsilon()) || (mcDirection2.GetMagnitude() < std::numeric_limits<float>::epsilon()))
        return false;

    mcDirection1 = mcDirection1 * (1.f / mcDirection1.GetMagnitude());
    mcDirection2 = mcDirection2 * (1.f / mcDirection2.GetMagnitude());

    const CartesianVector displacement1(mcVertex2 - mcVertex1), displacement2(mcVertex1 - mcVertex2);
    const float impactT1((mcDirection1.GetCrossProduct(displacement1)).GetMagnitude());
    const float impactL1(mcDirection1.GetDotProduct(displacement1));
    const float impactT2((mcDirection2.GetCrossProduct(displacement2)).GetMagnitude());
    const float impactL2(mcDirection2.GetDotProduct(displacement2));
    const float openingAngle(mcDirection1.GetOpeningAngle(mcDirection2) / 3.14 * 180.f);

    /*
    std::cout << "impactT1: " << impactT1 << std::endl;
    std::cout << "impactL1: " << impactL1 << std::endl;
    std::cout << "impactT2: " << impactT2 << std::endl;
    std::cout << "impactL2: " << impactL2 << std::endl;
    std::cout << "openingangle: " << openingAngle << std::endl;
    */

    if (this->IsShower(pMCParticle1) && (impactT1 < m_maxImpactT) && (impactL1 > 0.f))
    {          
        if (impactL1 < 1.f)
        {
            if (openingAngle < m_minOpeningAngle)
                return false;
        }
        else
        {
            return false;
        }
    }

    if (this->IsShower(pMCParticle2) && (impactT2 < m_maxImpactT) && (impactL2 > 0.f))
    {
        if (impactL2 < 1.f)
        {
            if (openingAngle < m_minOpeningAngle)
                return false;
        }
        else
        {
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::GetDistinguishableChildren(const MCParticle *const pMCContaminant, const MCParticle *const pMCGamma, 
   MCParticleList &mcHierarchyList) const
{
    const MCParticleList &mcChildrenList(pMCContaminant->GetDaughterList());

    for (const MCParticle *const pMCChild : mcChildrenList)
    {
        if (pMCChild == pMCGamma)
            continue;

        if (std::find(mcHierarchyList.begin(), mcHierarchyList.end(), pMCChild) != mcHierarchyList.end())
            continue;

        if (!this->AreParticlesSeparated(pMCGamma, pMCChild))
            continue;        

        mcHierarchyList.push_back(pMCChild);

        this->GetDistinguishableChildren(pMCChild, pMCGamma, mcHierarchyList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::RemoveContaminantHits(const MCParticle *const pMCGamma, const ParticleFlowObject *const pGammaPfo, const MCParticleList &mcContaminantList,
    const LArMCParticleHelper::MCContributionMap &mcParticleHitMap) const
{
    const CartesianVector &mcVertex(pMCGamma->GetVertex());
    const CartesianVector mcVertexU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_U));
    const CartesianVector mcVertexV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_V));
    const CartesianVector mcVertexW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_W));

    CartesianVector mcDirection(pMCGamma->GetMomentum());

    if (mcDirection.GetMagnitude() < std::numeric_limits<float>::epsilon())
        return;

    mcDirection = mcDirection * (1.f / mcDirection.GetMagnitude());

    const CartesianVector mcDirectionU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_U));
    const CartesianVector mcDirectionV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_V));
    const CartesianVector mcDirectionW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_W));

    ClusterList twoDClusterList;
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_U, twoDClusterList);
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_V, twoDClusterList);
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_W, twoDClusterList);

    MCParticleToHitListMap mcParticleToHitListMapU, mcParticleToHitListMapV, mcParticleToHitListMapW;

    for (const Cluster *const pCluster : twoDClusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        const CaloHitList &isolatedHitList(pCluster->GetIsolatedCaloHitList());

        for (const CaloHit *const pIsolated : isolatedHitList)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                caloHitList.push_back(pIsolated);
        }

        const HitType hitType(caloHitList.front()->GetHitType());
        MCParticleToHitListMap &viewMCParticleToHitListMap(hitType == TPC_VIEW_U ? mcParticleToHitListMapU : hitType == TPC_VIEW_V ? mcParticleToHitListMapV : mcParticleToHitListMapW);
        const CartesianVector &viewMCVertex(hitType == TPC_VIEW_U ? mcVertexU : hitType == TPC_VIEW_V ? mcVertexV : mcVertexW);
        const CartesianVector &viewMCDirection(hitType == TPC_VIEW_U ? mcDirectionU : hitType == TPC_VIEW_V ? mcDirectionV : mcDirectionW);
        const std::string &clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

        CaloHitList removedIsolatedHits;

        // add isolated hits here
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (std::find(removedIsolatedHits.begin(), removedIsolatedHits.end(), pCaloHit) != removedIsolatedHits.end())
                continue;

            const MCParticle *pMCParticle(nullptr);

            try
            {
                pMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);
            }
            catch (const StatusCodeException &)
            {
                continue;
            }

            if (std::find(mcContaminantList.begin(), mcContaminantList.end(), pMCParticle) == mcContaminantList.end())
                continue;

            if (m_truncateMode)
            {
                const CartesianVector &position(pCaloHit->GetPositionVector());
                const CartesianVector displacement(position - viewMCVertex);
                const float impactL(viewMCDirection.GetDotProduct(displacement));

                if (impactL > 0.f)
                    continue;
            }

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

            CaloHitList clusterNormalHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterNormalHitList);
            const CaloHitList clusterIsolatedHitList(pCluster->GetIsolatedCaloHitList());

            const bool isIsolated(std::find(isolatedHitList.begin(), isolatedHitList.end(), pCaloHit) != isolatedHitList.end());

            if (!isIsolated && (clusterNormalHitList.size() == 1) && !(clusterIsolatedHitList.empty()))
            {
                for (const CaloHit * const pIsolatedHit : clusterIsolatedHitList)
                {
                    removedIsolatedHits.push_back(pIsolatedHit);
                    const StatusCode isolatedStatusCode(PandoraContentApi::RemoveIsolatedFromCluster(*this, pCluster, pIsolatedHit));

                    if (isolatedStatusCode != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "ISOBEL CANNOT REMOVE ISOLATED HIT?" << std::endl;
                        throw;
                    }
                }
            }

            const StatusCode statusCodeCluster(isIsolated ? PandoraContentApi::RemoveIsolatedFromCluster(*this, pCluster, pCaloHit) :
                PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));

            //std::cout << StatusCodeToString(statusCodeCluster) << std::endl;

            if (statusCodeCluster != STATUS_CODE_SUCCESS)
            {
                if (statusCodeCluster != STATUS_CODE_NOT_ALLOWED)
                {
                    std::cout << "CheatingGammaRefinementAlgorithm: cluster jam" << std::endl;
                    throw StatusCodeException(statusCodeCluster);
                }

                const StatusCode statusCodePfo(PandoraContentApi::RemoveFromPfo(*this, pGammaPfo, pCluster));
                const unsigned int nHits(LArPfoHelper::GetNumberOfTwoDHits(pGammaPfo));

                if (nHits == 0)
                    std::cout << "CheatingGammaRefinementAlgorithm: ISOBEL - PFO HAS ZERO HITS" << std::endl;

                if (statusCodePfo != STATUS_CODE_SUCCESS)
                {
                    std::cout << "CheatingGammaRefinementAlgorithm: pfo jam" << std::endl;
                    throw;
                }

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pCluster));
            }

            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            {
                std::cout << "CheatingGammaRefinementAlgorithm: CALO HIT IS NOT AVAILABLE!!" << std::endl;
                throw;
            }
    
            viewMCParticleToHitListMap[pMCParticle].push_back(pCaloHit);
        }
    }
     
    MCParticleVector mcContaminantVector(mcContaminantList.begin(), mcContaminantList.end());
    std::sort(mcContaminantVector.begin(), mcContaminantVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcContaminantVector)
    {
        if (mcParticleHitMap.find(pMCParticle) == mcParticleHitMap.end())
            continue;

        const CaloHitList &mcParticleHitList(mcParticleHitMap.at(pMCParticle));

        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = (this->IsShower(pMCParticle) ? 11 : 13);
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = pMCParticle->GetEnergy();
        pfoParameters.m_momentum = pMCParticle->GetMomentum();

        CaloHitList allHits(0);

        if (mcParticleToHitListMapU.find(pMCParticle) == mcParticleToHitListMapU.end())
            continue;

        if (mcParticleToHitListMapV.find(pMCParticle) == mcParticleToHitListMapV.end())
            continue;

        if (mcParticleToHitListMapW.find(pMCParticle) == mcParticleToHitListMapW.end())
            continue;

        allHits.insert(allHits.end(), mcParticleToHitListMapU.at(pMCParticle).begin(), mcParticleToHitListMapU.at(pMCParticle).end());
        allHits.insert(allHits.end(), mcParticleToHitListMapV.at(pMCParticle).begin(), mcParticleToHitListMapV.at(pMCParticle).end());
        allHits.insert(allHits.end(), mcParticleToHitListMapW.at(pMCParticle).begin(), mcParticleToHitListMapW.at(pMCParticle).end());

        if (allHits.size() < 100)
            continue;

        // completeness cut
        const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(allHits, mcParticleHitList));
        const float completeness(static_cast<float>(sharedHitList.size()) / static_cast<float>(mcParticleHitList.size()));

        if (completeness < m_creationCompletenessThreshold)
            continue;

        pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapU.at(pMCParticle), TPC_VIEW_U));
        pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapV.at(pMCParticle), TPC_VIEW_V));
        pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapW.at(pMCParticle), TPC_VIEW_W));

        const PfoList *pTemporaryList(nullptr);
        std::string temporaryListName, currentListName;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<ParticleFlowObject>(*this, currentListName));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::CreateTemporaryListAndSetCurrent<PfoList>(*this, pTemporaryList, temporaryListName));

        const ParticleFlowObject *pPfo(nullptr);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
            PandoraContentApi::SaveList<ParticleFlowObject>(*this, temporaryListName, this->IsShower(pMCParticle) ? "ShowerParticles3D" : "TrackParticles3D"));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, currentListName));
      
        
        //PfoList createdDisplay;
        //createdDisplay.push_back(pPfo);
        //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &createdDisplay, "CREATED DISPLAY", BLUE, true, true);
        //CartesianVector vertexPosition(pMCParticle->GetVertex());
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vertexPosition, "CREATED VERTEX", BLUE, 2);
        //this->AreParticlesSeparated(pMCGamma, pMCParticle);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        //created.push_back(pPfo);
    }
    

    //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &created, "CREATED", BLUE, true, true);
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *CheatingGammaRefinementAlgorithm::CreateCluster(const MCParticle *const pMCParticle, const CaloHitList &caloHitList, const HitType hitType) const
{
    const ClusterList *pTemporaryList(nullptr);
    std::string temporaryListName, currentListName;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentListName));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*this, pTemporaryList, temporaryListName));

    const Cluster *pCluster(nullptr);
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = caloHitList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = this->IsShower(pMCParticle) ? 11 : 13;

    if (metadata.m_particleId.IsInitialized())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));

    const std::string &clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, temporaryListName, clusterListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentListName));

    return pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingGammaRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TruncateMode", m_truncateMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "RemoveHierarchy", m_removeHierarchy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MaxImpactT", m_maxImpactT));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinOpeningAngle", m_minOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinGammaCompleteness", m_minGammaCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "CreationCompletenessThreshold", m_creationCompletenessThreshold));


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::handle_eptr(std::exception_ptr eptr)
{
    try
    {
        if (eptr)
        {
            std::rethrow_exception(eptr);
        }
    }
    catch(const std::exception& e)
    {
        std::cout << "Caught exception \"" << e.what() << "\"\n";
    }
}

} // namespace lar_content


