/**
 *  @file   LArContent/src/LArMonitoring/ParticleMonitoringAlgorithm.cc
 * 
 *  @brief  Implementation of the particle monitoring algorithm.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"

#include "LArMonitoring/ParticleMonitoringAlgorithm.h"

using namespace pandora;

namespace lar
{

ParticleMonitoringAlgorithm::ParticleMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(Create());
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleMonitoringAlgorithm::~ParticleMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleMonitoringAlgorithm::Run()
{
    // Tree elements
    int nMCParticles(0), nPfos(0);
    IntVector pidVector, neutrinoVector, nMCHitsVector, nPfoHitsVector, nMatchedHitsVector;
    FloatVector completenessVector, purityVector, pxVector, pyVector, pzVector, thetaVector, energyVector, pVector;
    FloatVector vtxXPosVector, vtxYPosVector, vtxZPosVector, endXPosVector, endYPosVector, endZPosVector, mcVtxXPosVector, mcVtxYPosVector, mcVtxZPosVector, mcEndXPosVector, mcEndYPosVector, mcEndZPosVector;

    try
    {
        // Input lists
        const MCParticleList *pMCParticleList = NULL;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        const CaloHitList *pCaloHitList = NULL;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

        const PfoList *pPfoList = NULL;
        (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);
        const PfoList pfoList((NULL != pPfoList) ? PfoList(*pPfoList) : PfoList());
        nPfos = pfoList.size();

        // Obtain required information
        UidRelationMap uidToPrimaryMap;
        UidToMCParticleMap uidToMCParticleMap;
        this->GetMCParticleMaps(pMCParticleList, uidToPrimaryMap, uidToMCParticleMap);

        UidToPfoMap uidToPfoMap;
        ContributionMap mcPfoContributionMap;
        this->GetMCParticleToPfoMatches(pCaloHitList, pfoList, uidToPrimaryMap, uidToPfoMap, mcPfoContributionMap);

        ContributionMap mcHitContributionMap;
        this->GetMCParticleToCaloHitMatches(pCaloHitList, uidToPrimaryMap, mcHitContributionMap);

        // Populate tree elements
        for (ContributionMap::const_iterator iter = mcHitContributionMap.begin(), iterEnd = mcHitContributionMap.end(); iter != iterEnd; ++iter)
        {
            const Uid uid(iter->first);
            const float nMCHits(iter->second);

            UidToMCParticleMap::const_iterator mcIter = uidToMCParticleMap.find(uid);

            if (uidToMCParticleMap.end() == mcIter)
                return STATUS_CODE_FAILURE;

            const MCParticle *pMCParticle3D(mcIter->second);
            ++nMCParticles;

            neutrinoVector.push_back(LArMCParticleHelper::GetPrimaryNeutrino(pMCParticle3D));

            pidVector.push_back(pMCParticle3D->GetParticleId());
            pxVector.push_back(pMCParticle3D->GetMomentum().GetX());
            pyVector.push_back(pMCParticle3D->GetMomentum().GetY());
            pzVector.push_back(pMCParticle3D->GetMomentum().GetZ());
            pVector.push_back(pMCParticle3D->GetMomentum().GetMagnitude());
            energyVector.push_back(pMCParticle3D->GetEnergy());
            thetaVector.push_back(pMCParticle3D->GetMomentum().GetOpeningAngle(CartesianVector(0., 0., 1.)));

            mcVtxXPosVector.push_back(pMCParticle3D->GetVertex().GetX());
            mcVtxYPosVector.push_back(pMCParticle3D->GetVertex().GetY());
            mcVtxZPosVector.push_back(pMCParticle3D->GetVertex().GetZ());
            mcEndXPosVector.push_back(pMCParticle3D->GetEndpoint().GetX());
            mcEndYPosVector.push_back(pMCParticle3D->GetEndpoint().GetY());
            mcEndZPosVector.push_back(pMCParticle3D->GetEndpoint().GetZ());

            ContributionMap::const_iterator cIter = mcPfoContributionMap.find(uid);
            const int nMatchedHits((mcPfoContributionMap.end() == cIter) ? 0 : cIter->second);

            UidToPfoMap::const_iterator pfoIter = uidToPfoMap.find(uid);
            const ParticleFlowObject *pPfo((uidToPfoMap.end() == pfoIter) ? NULL : pfoIter->second);
            int nPfoHits(0);

            if (NULL != pPfo)
            {
                const ClusterList &clusterList(pPfo->GetClusterList());

                for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
                    nPfoHits += (*cIter)->GetNCaloHits() + (*cIter)->GetNIsolatedCaloHits();

                const CartesianVector vertex(this->GetSpacePoint(clusterList, pMCParticle3D->GetVertex(), true));
                const CartesianVector endpoint(this->GetSpacePoint(clusterList, pMCParticle3D->GetEndpoint(), false));
                vtxXPosVector.push_back(vertex.GetX());
                vtxYPosVector.push_back(vertex.GetY());
                vtxZPosVector.push_back(vertex.GetZ());
                endXPosVector.push_back(endpoint.GetX());
                endYPosVector.push_back(endpoint.GetY());
                endZPosVector.push_back(endpoint.GetZ());
            }
            else
            {
                vtxXPosVector.push_back(-std::numeric_limits<float>::max());
                vtxYPosVector.push_back(-std::numeric_limits<float>::max());
                vtxZPosVector.push_back(-std::numeric_limits<float>::max());
                endXPosVector.push_back(-std::numeric_limits<float>::max());
                endYPosVector.push_back(-std::numeric_limits<float>::max());
                endZPosVector.push_back(-std::numeric_limits<float>::max());
            }

            nMCHitsVector.push_back(nMCHits);
            nPfoHitsVector.push_back(nPfoHits);
            nMatchedHitsVector.push_back(nMatchedHits);

            const float purity((nPfoHits == 0) ? 0 : static_cast<float>(nMatchedHits) / static_cast<float>(nPfoHits));
            purityVector.push_back(purity);
            const float completeness((nPfoHits == 0) ? 0 : static_cast<float>(nMatchedHits) / static_cast<float>(nMCHits));
            completenessVector.push_back(completeness);
        }
    }
    catch (const StatusCodeException &except)
    {
        std::cout << "ParticleMonitoringAlgorithm: exception on loading event objects: " << except.what() << std::endl;
    }

    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "nMCParticles", nMCParticles));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "nPfos", nPfos));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "pdg", &pidVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "nupdg", &neutrinoVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "px", &pxVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "py", &pyVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "pz", &pzVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "pTot", &pVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "energy", &energyVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "theta", &thetaVector));

    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "mcVtxXPos", &mcVtxXPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "mcVtxYPos", &mcVtxYPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "mcVtxZPos", &mcVtxZPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "mcEndXPos", &mcEndXPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "mcEndYPos", &mcEndYPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "mcEndZPos", &mcEndZPosVector));

    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "vtxXPos", &vtxXPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "vtxYPos", &vtxYPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "vtxZPos", &vtxZPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "endXPos", &endXPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "endYPos", &endYPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "endZPos", &endZPosVector));

    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "completeness", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "purity", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "nMCHits", &nMCHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "nPfoHits", &nPfoHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "nMatchedHits", &nMatchedHitsVector));

    PANDORA_MONITORING_API(FillTree(m_treeName.c_str()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetMCParticleMaps(const MCParticleList *const pMCParticleList, UidRelationMap &uidToPrimaryMap,
    UidToMCParticleMap &uidToMCParticleMap) const
{
    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        const MCParticle *pMCParticle = *iter;
        uidToMCParticleMap[pMCParticle->GetUid()] = pMCParticle;

        if (pMCParticle->GetParentList().empty() || LArMCParticleHelper::IsNeutrinoInduced(pMCParticle))
        {
            if (!LArMCParticleHelper::IsNeutrino(pMCParticle))
                uidToPrimaryMap[pMCParticle->GetUid()] = pMCParticle->GetUid();
        }

        else
        {
            const MCParticle *pParentMCParticle = pMCParticle;

            while (true)
            {
                pParentMCParticle = *(pParentMCParticle->GetParentList().begin());

                if (uidToPrimaryMap.find(pParentMCParticle->GetUid()) != uidToPrimaryMap.end())
                {
                    uidToPrimaryMap[pMCParticle->GetUid()] = uidToPrimaryMap[pParentMCParticle->GetUid()];
                    break;
                }
                else if (pParentMCParticle->GetParentList().empty() || LArMCParticleHelper::IsNeutrinoInduced(pParentMCParticle))
                {
                    uidToPrimaryMap[pMCParticle->GetUid()] = pParentMCParticle->GetUid();
                    break;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetMCParticleToPfoMatches(const CaloHitList *const pCaloHitList, const PfoList &pfoList, 
    const UidRelationMap &uidToPrimaryMap, UidToPfoMap &uidToPfoMap, ContributionMap &contributionMap) const
{
    for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *pPfo = *iter;
        const ClusterList &clusterList(pPfo->GetClusterList());
        ContributionMap pfoContributionMap;

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;

            CaloHitList clusterHits;
            pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterHits);
            clusterHits.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

            for (CaloHitList::const_iterator hIter = clusterHits.begin(), hIterEnd = clusterHits.end(); hIter != hIterEnd; ++hIter)
            {
                try
                {
                    CaloHit *pCaloHit = *hIter;

                    if (pCaloHitList->count(pCaloHit) == 0)
                        continue;

                    const MCParticle *pHitMCParticle(pCaloHit->GetMainMCParticle());
                    UidRelationMap::const_iterator rIter = uidToPrimaryMap.find(pHitMCParticle->GetUid());

                    if (uidToPrimaryMap.end() == rIter)
                        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                    pfoContributionMap[rIter->second] += 1;
                }
                catch (const StatusCodeException &)
                {
                }
            }
        }

        float maxHits = -1.f;
        Uid biggestContributor = NULL;

        for (ContributionMap::const_iterator cIter = pfoContributionMap.begin(), cIterEnd = pfoContributionMap.end(); cIter != cIterEnd; ++cIter)
        {
            if (cIter->second > maxHits)
            {
                maxHits = cIter->second;
                biggestContributor = cIter->first;
            }
        }

        if (biggestContributor && (contributionMap[biggestContributor] < maxHits))
        {
            contributionMap[biggestContributor] = maxHits;
            uidToPfoMap[biggestContributor] = pPfo;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetMCParticleToCaloHitMatches(const CaloHitList *const pCaloHitList, const UidRelationMap &uidToPrimaryMap,
    ContributionMap &contributionMap) const
{
    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            CaloHit *pCaloHit = *iter;
            const MCParticle *pMCParticle(pCaloHit->GetMainMCParticle());
            UidRelationMap::const_iterator rIter = uidToPrimaryMap.find(pMCParticle->GetUid());

            if (uidToPrimaryMap.end() == rIter)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            contributionMap[rIter->second] += 1;
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector ParticleMonitoringAlgorithm::GetSpacePoint(const ClusterList &clusterList, const CartesianVector &referencePoint, const bool useInnerLayer) const
{
    CartesianVector spacePoint(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    ClusterList clusterListU, clusterListV, clusterListW;

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        // TODO
    }

    return spacePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
