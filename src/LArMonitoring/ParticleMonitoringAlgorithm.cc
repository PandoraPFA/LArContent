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
        (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);
        const PfoList pfoList((NULL != pPfoList) ? PfoList(*pPfoList) : PfoList());
        nPfos = pfoList.size();

        // Match Pfos and MC Particles
        MCRelationMap mcPrimaryMap;      // [particles -> primary particle]
        this->GetMCParticleMaps(pMCParticleList, mcPrimaryMap);

        MCContributionMap trueHitMap;    // [particle -> number of true hits]
        CaloHitToMCMap mcHitMap;         // [hit -> parent primary]
        this->GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, mcHitMap, trueHitMap);

        PfoContributionMap recoHitMap;   // [pfo -> number of reco hits]
        CaloHitToPfoMap pfoHitMap;       // [hit -> parent pfo]
        this->GetPfoToCaloHitMatches(pCaloHitList, pfoList, pfoHitMap, recoHitMap);

        MCContributionMap matchedHitMap; // [particle -> number of matched hits with best pfo]
        MCToPfoMap matchedPfoMap;        // [particle -> best pfo]
        this->GetMCParticleToPfoMatches(pCaloHitList, pfoList, mcHitMap, matchedPfoMap, matchedHitMap);

        nMCParticles = trueHitMap.size();

        // Write out reco/true information
        for (MCContributionMap::const_iterator iter = trueHitMap.begin(), iterEnd = trueHitMap.end(); iter != iterEnd; ++iter)
        {
            const MCParticle *pMCParticle3D(iter->first);
            const int nMCHits(iter->second);

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

            int nPfoHits(0);

            MCToPfoMap::const_iterator pIter1 = matchedPfoMap.find(pMCParticle3D);
            MCContributionMap::const_iterator pIter2 = matchedHitMap.find(pMCParticle3D);

            const ParticleFlowObject *pPfo((matchedPfoMap.end() == pIter1) ? NULL : pIter1->second);
            const int nMatchedHits((matchedHitMap.end() == pIter2) ? 0 : pIter2->second);

            if (NULL != pPfo)
            {
                PfoContributionMap::const_iterator pIter3 = recoHitMap.find(pPfo);
                if (recoHitMap.end() == pIter3)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                nPfoHits = pIter3->second;

                const CartesianVector vertex(pMCParticle3D->GetVertex());
                const CartesianVector endpoint(pMCParticle3D->GetEndpoint());
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
            const float completeness((nPfoHits == 0) ? 0 : static_cast<float>(nMatchedHits) / static_cast<float>(nMCHits));

            purityVector.push_back(purity);
            completenessVector.push_back(completeness);
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if ((STATUS_CODE_NOT_INITIALIZED != statusCodeException.GetStatusCode()) &&
            (STATUS_CODE_NOT_FOUND != statusCodeException.GetStatusCode()))
            throw statusCodeException;
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

void ParticleMonitoringAlgorithm::GetMCParticleMaps(const MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap) const
{
    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        const MCParticle *pMCParticle = *iter;

        if (pMCParticle->GetParentList().empty() || LArMCParticleHelper::IsNeutrinoInduced(pMCParticle))
        {
            if (!LArMCParticleHelper::IsNeutrino(pMCParticle))
                mcPrimaryMap[pMCParticle] = pMCParticle;
        }

        else
        {
            const MCParticle *pParentMCParticle = pMCParticle;

            while (true)
            {
                pParentMCParticle = *(pParentMCParticle->GetParentList().begin());

                if (mcPrimaryMap.find(pParentMCParticle) != mcPrimaryMap.end())
                {
                    mcPrimaryMap[pMCParticle] = mcPrimaryMap[pParentMCParticle];
                    break;
                }
                else if (pParentMCParticle->GetParentList().empty() || LArMCParticleHelper::IsNeutrinoInduced(pParentMCParticle))
                {
                    mcPrimaryMap[pMCParticle] = pParentMCParticle;
                    break;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::GetMCParticleToPfoMatches(const CaloHitList *const pCaloHitList, const PfoList &pfoList,
    const CaloHitToMCMap &mcHitMap, MCToPfoMap &outputPrimaryMap, MCContributionMap &outputContributionMap) const
{
    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pPfo = *pIter;

        CaloHitList clusterHits;
        this->CollectCaloHits(pPfo, clusterHits);

        MCContributionMap truthContributionMap;

        for (CaloHitList::const_iterator hIter = clusterHits.begin(), hIterEnd = clusterHits.end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;

            if (pCaloHitList->count(pCaloHit) == 0)
                continue;

            CaloHitToMCMap::const_iterator mcIter = mcHitMap.find(pCaloHit);
            if (mcIter == mcHitMap.end())
                continue;

            const MCParticle *pPrimaryParticle = mcIter->second;
            truthContributionMap[pPrimaryParticle] += 1;
        }

        int maxHits(0);
        MCParticle *biggestContributor = NULL;

        for (MCContributionMap::const_iterator cIter = truthContributionMap.begin(), cIterEnd = truthContributionMap.end(); cIter != cIterEnd; ++cIter)
        {
            if (cIter->second > maxHits)
            {
                maxHits = cIter->second;
                biggestContributor = const_cast<MCParticle*>(cIter->first);
            }
        }

        if (biggestContributor)
        {
            int matchedHits(0.f);

            MCContributionMap::const_iterator cIter = outputContributionMap.find(biggestContributor);
            if (outputContributionMap.end() != cIter)
                matchedHits = cIter->second;

            if (maxHits > matchedHits)
            {
                outputPrimaryMap[biggestContributor] = pPfo;
                outputContributionMap[biggestContributor] = maxHits;
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
            CaloHit *pCaloHit = *iter;
            const MCParticle *pHitParticle(pCaloHit->GetMainMCParticle());

            MCRelationMap::const_iterator mcIter = mcPrimaryMap.find(pHitParticle);
            if (mcIter == mcPrimaryMap.end())
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const MCParticle *pPrimaryParticle = mcIter->second;
            mcContributionMap[pPrimaryParticle] += 1;
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
        const ParticleFlowObject *pPfo = *pIter;
        const ClusterList &pfoClusterList(pPfo->GetClusterList());
        ClusterList clusterList(pfoClusterList.begin(), pfoClusterList.end());

        if (m_useDaughterPfos)
            this->CollectDaughterClusters(pPfo, clusterList);

        int nPfoHits(0);

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;

            CaloHitList clusterHits;
            pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterHits);
            clusterHits.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

            for (CaloHitList::const_iterator hIter = clusterHits.begin(), hIterEnd = clusterHits.end(); hIter != hIterEnd; ++hIter)
            {
                CaloHit *pCaloHit = *hIter;

                if (pCaloHitList->count(pCaloHit) == 0)
                    continue;

                pfoHitMap[pCaloHit] = pPfo;
                ++nPfoHits;
            }
        }

        pfoContributionMap[pPfo] = nPfoHits;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::CollectCaloHits(const ParticleFlowObject *const pParentPfo, CaloHitList &caloHitList) const
{
    ClusterList clusterList;
    const ClusterList &pfoClusterList(pParentPfo->GetClusterList());
    clusterList.insert(pfoClusterList.begin(), pfoClusterList.end());

    if (m_useDaughterPfos)
    {
        if (!pParentPfo->GetParentPfoList().empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        this->CollectDaughterClusters(pParentPfo, clusterList);
    }

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        caloHitList.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleMonitoringAlgorithm::CollectDaughterClusters(const ParticleFlowObject *const pParentPfo, ClusterList &clusterList) const
{
    const PfoList &daughterPfoList(pParentPfo->GetDaughterPfoList());
    for (PfoList::const_iterator dIter = daughterPfoList.begin(), dIterEnd = daughterPfoList.end(); dIter != dIterEnd; ++dIter)
    {
        const ParticleFlowObject *pDaughterPfo = *dIter;
        const ClusterList &pfoClusterList(pDaughterPfo->GetClusterList());
        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;

            if (clusterList.count(pCluster))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            clusterList.insert(pCluster);
        }

        this->CollectDaughterClusters(pDaughterPfo, clusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    m_useDaughterPfos = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CountDaughters", m_useDaughterPfos));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
