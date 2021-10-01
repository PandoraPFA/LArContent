/**
 *  @file   larpandoracontent/LArCheating/CheatingCCLeptonClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingCCLeptonClusterCreationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingCCLeptonClusterCreationAlgorithm::CheatingCCLeptonClusterCreationAlgorithm() : m_collapseToPrimaryMCParticles(false), m_collectAllElectronHits(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCCLeptonClusterCreationAlgorithm::Run()
{
    //PandoraMonitoringApi::Create(this->GetPandora());
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    MCParticleToHitListMap mcParticleToHitListMap;
    this->GetMCParticleToHitListMap(mcParticleToHitListMap);
    this->CreateClusters(mcParticleToHitListMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCLeptonClusterCreationAlgorithm::GetMCParticleToHitListMap(MCParticleToHitListMap &mcParticleToHitListMap) const
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;

    if (m_collapseToPrimaryMCParticles)
    {
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    }

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "Input", pMCParticleList));

    const bool isNueCC(LArMCParticleHelper::IsCCNuEvent(pMCParticleList, 12));
    const bool isNumuCC(LArMCParticleHelper::IsCCNuEvent(pMCParticleList, 14));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    if (isNueCC)
    {
        MCParticleVector mcParticleVector(pMCParticleList->begin(), pMCParticleList->end());
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            const bool isLeadingElectron((std::abs(pMCParticle->GetParticleId()) == 11) && (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle));

            if (isLeadingElectron)
                this->ClusterLeadingElectrons(pMCParticle, pCaloHitList, mcParticleToHitListMap);
        }
    }

    CaloHitList hitsToCluster;

    for (auto &entry : mcParticleToHitListMap)
        hitsToCluster.insert(hitsToCluster.begin(), entry.second.begin(), entry.second.end());

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (std::find(hitsToCluster.begin(), hitsToCluster.end(), pCaloHit) != hitsToCluster.end())
            continue;

        try
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            this->SimpleMCParticleCollection(pCaloHit, mcPrimaryMap, isNueCC, isNumuCC, mcParticleToHitListMap);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCLeptonClusterCreationAlgorithm::ClusterLeadingElectrons(const MCParticle *const pMCLeadingElectron, const CaloHitList *pCaloHitList, 
    MCParticleToHitListMap &mcParticleToHitListMap) const
{
    CaloHitList hitsToCluster;

    for (auto &entry : mcParticleToHitListMap)
        hitsToCluster.insert(hitsToCluster.begin(), entry.second.begin(), entry.second.end());

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (std::find(hitsToCluster.begin(), hitsToCluster.end(), pCaloHit) != hitsToCluster.end())
            continue;

        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            if (pMCParticle != pMCLeadingElectron)
                continue;

            hitsToCluster.push_back(pCaloHit);

            mcParticleToHitListMap[pMCLeadingElectron].push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }

    if (!m_collectAllElectronHits)
        return;

    if (mcParticleToHitListMap.find(pMCLeadingElectron) == mcParticleToHitListMap.end())
        return;

    const CartesianVector mcVertexPosition(pMCLeadingElectron->GetVertex());
    const CartesianVector mcVertexPositionU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertexPosition, TPC_VIEW_U));
    const CartesianVector mcVertexPositionV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertexPosition, TPC_VIEW_V));
    const CartesianVector mcVertexPositionW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertexPosition, TPC_VIEW_W));

    CartesianVector closestPositionU(0.f, 0.f, 0.f), closestPositionV(0.f, 0.f, 0.f), closestPositionW(0.f, 0.f, 0.f);
    float closestDistanceU(std::numeric_limits<float>::max()), closestDistanceV(std::numeric_limits<float>::max()), closestDistanceW(std::numeric_limits<float>::max());
    bool foundU(false), foundV(false), foundW(false);

    const CaloHitList leadingElectronHits(mcParticleToHitListMap.at(pMCLeadingElectron));

    for (const CaloHit *const pCaloHit : leadingElectronHits)
    {
        const HitType hitType(pCaloHit->GetHitType());
        const CartesianVector projectedMCParticleVertexPosition(hitType == TPC_VIEW_U ? mcVertexPositionU : hitType == TPC_VIEW_V ? mcVertexPositionV : mcVertexPositionW);
        float &closestDistance(hitType == TPC_VIEW_U ? closestDistanceU : hitType == TPC_VIEW_V ? closestDistanceV : closestDistanceW);

        const float separation((pCaloHit->GetPositionVector() - projectedMCParticleVertexPosition).GetMagnitudeSquared());

        if (separation < closestDistance)
        {
            closestDistance = separation;

            CartesianVector &closestPosition(hitType == TPC_VIEW_U ? closestPositionU : hitType == TPC_VIEW_V ? closestPositionV : closestPositionW);
            closestPosition = pCaloHit->GetPositionVector();

            bool &found(hitType == TPC_VIEW_U ? foundU : hitType == TPC_VIEW_V ? foundV : foundW);
            found = true;
        }
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (std::find(hitsToCluster.begin(), hitsToCluster.end(), pCaloHit) != hitsToCluster.end())
            continue;

        MCParticleVector contributingMCParticleVector;

        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap())
            contributingMCParticleVector.push_back(mapEntry.first);

        std::sort(contributingMCParticleVector.begin(), contributingMCParticleVector.end(), PointerLessThan<MCParticle>());

        for (const MCParticle *const pMCParticle : contributingMCParticleVector)
        {
            if (pMCParticle == pMCLeadingElectron)
            {
                const HitType hitType(pCaloHit->GetHitType());

                const bool found(hitType == TPC_VIEW_U ? foundU : hitType == TPC_VIEW_V ? foundV : foundW);

                if (!found)
                    break;

                const CartesianVector projectedMCParticleVertexPosition(hitType == TPC_VIEW_U ? mcVertexPositionU : hitType == TPC_VIEW_V ? mcVertexPositionV : mcVertexPositionW);
                const CartesianVector closestPosition(hitType == TPC_VIEW_U ? closestPositionU : hitType == TPC_VIEW_V ? closestPositionV : closestPositionW);

                const float minZ(std::min(projectedMCParticleVertexPosition.GetZ(), closestPosition.GetZ()));
                const float maxZ(std::max(projectedMCParticleVertexPosition.GetZ(), closestPosition.GetZ()));

                if (pCaloHit->GetPositionVector().GetZ() < minZ)
                    break;

                if (pCaloHit->GetPositionVector().GetZ() > maxZ)
                    break;

                //CartesianVector hitPosition(pCaloHit->GetPositionVector());
                //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "ADDED", GREEN, 2);

                mcParticleToHitListMap[pMCLeadingElectron].push_back(pCaloHit);

                break;
            }
        }
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCLeptonClusterCreationAlgorithm::SimpleMCParticleCollection(const CaloHit *const pCaloHit,
    const LArMCParticleHelper::MCRelationMap &mcPrimaryMap, const bool isNueCC, const bool isNumuCC, MCParticleToHitListMap &mcParticleToHitListMap) const
{
    /*    
    if (isNueCC && m_collectAllElectronHits)
    {
        MCParticleVector contributingMCParticleVector;
        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap())
            contributingMCParticleVector.push_back(mapEntry.first);

        std::sort(contributingMCParticleVector.begin(), contributingMCParticleVector.end(), PointerLessThan<MCParticle>());

        for (const MCParticle *const pMCParticle : contributingMCParticleVector)
        {
            const bool isLeadingElectron((std::abs(pMCParticle->GetParticleId()) == 11) && (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle));

            if (isLeadingElectron)
            {
                mcParticleToHitListMap[pMCParticle].push_back(pCaloHit);
                return;
            }
        }
    }
    */

    const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

    if (!this->SelectMCParticlesForClustering(pMCParticle, isNueCC, isNumuCC))
        return;

    if (m_collapseToPrimaryMCParticles)
    {
        LArMCParticleHelper::MCRelationMap::const_iterator primaryIter = mcPrimaryMap.find(pMCParticle);

        if (mcPrimaryMap.end() == primaryIter)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        pMCParticle = primaryIter->second;
    }

    mcParticleToHitListMap[pMCParticle].push_back(pCaloHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

    bool CheatingCCLeptonClusterCreationAlgorithm::SelectMCParticlesForClustering(const MCParticle *const pMCParticle, const bool isNueCC, const bool isNumuCC) const
{
    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle));
    const int pdg(std::abs(pParentMCParticle->GetParticleId()));
    const bool selected(std::find(m_particleIdList.begin(), m_particleIdList.end(), pdg) != m_particleIdList.end());

    if (!selected)
        return false;

    const bool isLeadingLepton((pdg == E_MINUS) || (pdg == MU_MINUS));

    if (isLeadingLepton)
    {
        if ((pdg == E_MINUS) && isNueCC)
            return true;

        if ((pdg == MU_MINUS) && isNumuCC)
            return true;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCLeptonClusterCreationAlgorithm::CreateClusters(const MCParticleToHitListMap &mcParticleToHitListMap) const
{
    MCParticleVector mcParticleVector;
    for (const auto &mapEntry : mcParticleToHitListMap)
        mcParticleVector.push_back(mapEntry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleVector)
    {
        const CaloHitList &caloHitList(mcParticleToHitListMap.at(pMCParticle));

        if (caloHitList.empty())
            continue;

        const Cluster *pCluster(nullptr);
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList = caloHitList;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

        PandoraContentApi::Cluster::Metadata metadata;

        switch (pMCParticle->GetParticleId())
        {
            case PHOTON:
            case E_PLUS:
            case E_MINUS:
            case MU_PLUS:
            case MU_MINUS:
                metadata.m_particleId = pMCParticle->GetParticleId();
                break;
            default:
                break;
        }

        if (metadata.m_particleId.IsInitialized())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCCLeptonClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CollapseToPrimaryMCParticles", m_collapseToPrimaryMCParticles));

    if (m_collapseToPrimaryMCParticles)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CollectAllElectronHits", m_collectAllElectronHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ParticleIdList", m_particleIdList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
