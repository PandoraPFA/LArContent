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

#include "larpandoracontent/LArCheating/CheatingGammaRefinementAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingGammaRefinementAlgorithm::CheatingGammaRefinementAlgorithm() :
    m_truncateMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingGammaRefinementAlgorithm::Run()
{
    //PandoraMonitoringApi::Create(this->GetPandora());
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    const PfoList *pShowerPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_showerPfoListName, pShowerPfoList));

    if (!pShowerPfoList || pShowerPfoList->empty())
    {
        std::cout << "CheatingGammaRefinementAlgorithm: No shower pfo list found, returning..." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    PfoVector showerVector(pShowerPfoList->begin(), pShowerPfoList->end());
    std::sort(showerVector.begin(), showerVector.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo : showerVector)
    {
        if (!this->IsMatchedToGamma(pPfo))
            continue;
        
        PfoList pfoDisplay;
        pfoDisplay.push_back(pPfo);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoDisplay, "matched gamma pfo", RED, true, true);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
        MCParticleList contaminantMCParticles;
        this->FindContaminantMCParticles(pPfo, contaminantMCParticles);

        PandoraMonitoringApi::ViewEvent(this->GetPandora());

        if (contaminantMCParticles.empty())
        {
            //PandoraMonitoringApi::Delete(this->GetPandora());
            return STATUS_CODE_SUCCESS;
        }

        this->FilterContaminantMCParticleList(pPfo, contaminantMCParticles);

        if (contaminantMCParticles.empty())
        {
            //PandoraMonitoringApi::Delete(this->GetPandora());
            return STATUS_CODE_SUCCESS;
        }

        this->RemoveContaminantHits(pPfo, contaminantMCParticles);
    }

    //PandoraMonitoringApi::Delete(this->GetPandora());

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingGammaRefinementAlgorithm::IsMatchedToGamma(const ParticleFlowObject *const pPfo) const
{
    /*
    // could also make sure that the MC gamma is a primary?
    */

    const MCParticle *pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
    const int pdg(std::abs(pMCParticle->GetParticleId()));

    if (pdg != 22)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::FindContaminantMCParticles(const ParticleFlowObject *const pGammaPfo, MCParticleList &contaminantMCParticleList) const
{
    const MCParticle *pMCGamma(LArMCParticleHelper::GetMainMCParticle(pGammaPfo));
    const CartesianVector mcVertex(pMCGamma->GetVertex());
    const CartesianVector mcVertexU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_U));
    const CartesianVector mcVertexV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_V));
    const CartesianVector mcVertexW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_W));

    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexU, "mc vertex", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexV, "mc vertex", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexW, "mc vertex", BLACK, 2);

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

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

    CartesianPointVector hitPositions;
    for (const Cluster *const pCluster : twoDClusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList); 

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

                hitPositions.push_back(pCaloHit->GetPositionVector());

                if (std::find(contaminantMCParticleList.begin(), contaminantMCParticleList.end(), pMCParticle) != contaminantMCParticleList.end())
                    continue;

                std::cout << "UNREFINED MC LIST: " << pMCParticle->GetParticleId() << std::endl;
                contaminantMCParticleList.push_back(pMCParticle);
            }
            catch (const StatusCodeException &)
            {
            }
        }
    }

    
    for (const CartesianVector jam : hitPositions)
    {
        //std::cout << &jam << std::endl;
        std::cout << "1111111111111" << std::endl;
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &jam, "remove", VIOLET, 1);
        std::cout << "2222222222222" << std::endl;
    }
    
    //if (hitPositions.size())
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::FilterContaminantMCParticleList(const ParticleFlowObject *const pGammaPfo, MCParticleList &mcContaminantList) const
{
    const MCParticle *pMCGamma(LArMCParticleHelper::GetMainMCParticle(pGammaPfo));

    MCParticleVector mcContaminantVector(mcContaminantList.begin(), mcContaminantList.end());
    std::sort(mcContaminantVector.begin(), mcContaminantVector.end(), LArMCParticleHelper::SortByMomentum);

    mcContaminantList.clear();

    for (const MCParticle *const pMCContaminant : mcContaminantVector)
    {
        if (this->IsShower(pMCContaminant) && !this->AreShowersSeparated(pMCGamma, pMCContaminant))
            continue;

        std::cout << "REFINED MC LIST: " << pMCContaminant->GetParticleId() << std::endl;
        mcContaminantList.push_back(pMCContaminant);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingGammaRefinementAlgorithm::IsShower(const MCParticle *const pMCParticle) const
{
  const int pdg(pMCParticle->GetParticleId());

  return ((E_MINUS == std::abs(pdg)) || (PHOTON == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingGammaRefinementAlgorithm::AreShowersSeparated(const MCParticle *const pMCParticle1, const MCParticle *const pMCParticle2) const
{
    const CartesianVector &mcVertex1(pMCParticle1->GetVertex()), &mcVertex2(pMCParticle2->GetVertex());
    CartesianVector mcDirection1(pMCParticle1->GetMomentum()), mcDirection2(pMCParticle2->GetMomentum());

    if ((mcDirection1.GetMagnitude() < std::numeric_limits<float>::epsilon()) || (mcDirection2.GetMagnitude() < std::numeric_limits<float>::epsilon()))
        return false;

    mcDirection1 = mcDirection1 * (1.f / mcDirection1.GetMagnitude());

    mcDirection2 = mcDirection2 * (1.f / mcDirection2.GetMagnitude());

    const float openingAngle(mcDirection1.GetOpeningAngle(mcDirection2) / 3.14 * 180.f);

    if (openingAngle < 10.f)
    {
        std::cout << "opening angle too low" << std::endl;
        return false;
    }

    const CartesianVector displacement1(mcVertex2 - mcVertex1), displacement2(mcVertex1 - mcVertex2);
    const float impactT1((mcDirection1.GetCrossProduct(displacement1)).GetMagnitude());
    const float impactT2((mcDirection2.GetCrossProduct(displacement2)).GetMagnitude());

    if ((impactT1 < 1.f) || (impactT2 < 1.f))
    {
        std::cout << "transverse impact parameter is too small" << std::endl;
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingGammaRefinementAlgorithm::RemoveContaminantHits(const ParticleFlowObject *const pGammaPfo, const MCParticleList &mcContaminantList) const
{
    const MCParticle *pMCGamma(LArMCParticleHelper::GetMainMCParticle(pGammaPfo));
    const CartesianVector mcVertex(pMCGamma->GetVertex());
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


    int hitsRemoved(0);
    for (const Cluster *const pCluster : twoDClusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        const HitType hitType(caloHitList.front()->GetHitType());

        MCParticleToHitListMap &viewMCParticleToHitListMap(hitType == TPC_VIEW_U ? mcParticleToHitListMapU : hitType == TPC_VIEW_V ? mcParticleToHitListMapV : mcParticleToHitListMapW);
        const CartesianVector &viewMCVertex(hitType == TPC_VIEW_U ? mcVertexU : hitType == TPC_VIEW_V ? mcVertexV : mcVertexW);
        const CartesianVector &viewMCDirection(hitType == TPC_VIEW_U ? mcDirectionU : hitType == TPC_VIEW_V ? mcDirectionV : mcDirectionW);
        const std::string &clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            try
            {
                const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

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

                const StatusCode statusCodeCluster(PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));

                if (statusCodeCluster != STATUS_CODE_SUCCESS)
                {
                    //std::cout << "CCCCCCCCCCCC" << std::endl;
                    if (statusCodeCluster != STATUS_CODE_NOT_ALLOWED)
                    {
                        std::cout << "jam" << std::endl;
                        throw StatusCodeException(statusCodeCluster);
                    }

                    const StatusCode statusCodePfo(PandoraContentApi::RemoveFromPfo(*this, pGammaPfo, pCluster));

                    unsigned int nHits(LArPfoHelper::GetNumberOfTwoDHits(pGammaPfo));

                    if (nHits == 0)
                        std::cout << "CheatingCCElectronRefinementAlgorithm: ISOBEL - PFO HAS ZERO HITS" << std::endl;

                    if (statusCodePfo != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "pfo jam" << std::endl;
                        throw;
                    }

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pCluster));
                }

                if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                    std::cout << "CALO HIT IS NOT AVAILABLE!!" << std::endl;
    
                ++hitsRemoved;
                    
                viewMCParticleToHitListMap[pMCParticle].push_back(pCaloHit);
            }
            catch (const StatusCodeException &)
            {
            }
        }
    }

    MCParticleVector mcContaminantVector(mcContaminantList.begin(), mcContaminantList.end());
    std::sort(mcContaminantVector.begin(), mcContaminantVector.end(), LArMCParticleHelper::SortByMomentum);

    PfoList created;

    int hitsCreated(0);
    for (const MCParticle *const pMCParticle : mcContaminantVector)
    {
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = (this->IsShower(pMCParticle) ? 11 : 13);
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = pMCParticle->GetEnergy();
        pfoParameters.m_momentum = pMCParticle->GetMomentum();

        if (mcParticleToHitListMapU.find(pMCParticle) != mcParticleToHitListMapU.end())
            pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapU.at(pMCParticle), TPC_VIEW_U));

        if (mcParticleToHitListMapV.find(pMCParticle) != mcParticleToHitListMapV.end())
            pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapV.at(pMCParticle), TPC_VIEW_V));

        if (mcParticleToHitListMapW.find(pMCParticle) != mcParticleToHitListMapW.end())
            pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapW.at(pMCParticle), TPC_VIEW_W));

        for (const Cluster *const jam : pfoParameters.m_clusterList)
        {
            CaloHitList jamList;
            jam->GetOrderedCaloHitList().FillCaloHitList(jamList);
            std::cout << "cluster created size: " << jamList.size() << std::endl;
            hitsCreated += jamList.size();
        }

        const PfoList *pTemporaryList(nullptr);
        //PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
        //PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, this->IsShower(pMCParticle) ? "ShowerParticles3D" : "TrackParticles3D"));

        std::string temporaryListName, currentListName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<ParticleFlowObject>(*this, currentListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::CreateTemporaryListAndSetCurrent<PfoList>(*this, pTemporaryList, temporaryListName));

        const ParticleFlowObject *pPfo(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
            PandoraContentApi::SaveList<ParticleFlowObject>(*this, temporaryListName, this->IsShower(pMCParticle) ? "ShowerParticles3D" : "TrackParticles3D"));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, currentListName));

        created.push_back(pPfo);
    }

    std::cout << "hitsRemoved: " << hitsRemoved << std::endl;
    std::cout << "hitsCreated: " << hitsCreated << std::endl;

    if (hitsCreated > 0)
    {
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &created, "CREATED", BLUE, true, true);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *CheatingGammaRefinementAlgorithm::CreateCluster(const MCParticle *const pMCParticle, const CaloHitList &caloHitList, const HitType hitType) const
{
    const std::string &clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

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

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, temporaryListName, currentListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentListName));

    return pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingGammaRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TruncateMode", m_truncateMode));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content


