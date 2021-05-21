/**
 *  @file   larpandoracontent/LArMonitoring/DeltaRayTensorVisualizationTool.cc
 *
 *  @brief  Implementation of the delta ray tensor visualization tool class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMonitoring/DeltaRayTensorVisualizationTool.h"

using namespace pandora;

namespace lar_content
{

DeltaRayTensorVisualizationTool::DeltaRayTensorVisualizationTool() :
    m_minClusterConnections(1),
    m_ignoreUnavailableClusters(true),
    m_showEachIndividualElement(false),
    m_showContext(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayTensorVisualizationTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterSet usedKeyClusters;
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (m_ignoreUnavailableClusters && !pKeyCluster->IsAvailable())
            continue;

        if (usedKeyClusters.count(pKeyCluster))
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, m_ignoreUnavailableClusters, elementList, nU, nV, nW);

        if ((nU < m_minClusterConnections) && (nV < m_minClusterConnections) && (nW < m_minClusterConnections))
            continue;

        if (nU * nV * nW == 0)
            continue;

        int counter(0);
        ClusterList allClusterListU, allClusterListV, allClusterListW;
        std::cout << " Connections: nU " << nU << ", nV " << nV << ", nW " << nW << ", nElements " << elementList.size() << std::endl;

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (allClusterListU.end() == std::find(allClusterListU.begin(), allClusterListU.end(), eIter->GetClusterU())) allClusterListU.push_back(eIter->GetClusterU());
            if (allClusterListV.end() == std::find(allClusterListV.begin(), allClusterListV.end(), eIter->GetClusterV())) allClusterListV.push_back(eIter->GetClusterV());
            if (allClusterListW.end() == std::find(allClusterListW.begin(), allClusterListW.end(), eIter->GetClusterW())) allClusterListW.push_back(eIter->GetClusterW());
            usedKeyClusters.insert(eIter->GetClusterU());

            std::cout << " Element " << counter++ << ": MatchedFraction " << eIter->GetOverlapResult().GetMatchedFraction()
                      << ", MatchedSamplingPoints " << eIter->GetOverlapResult().GetNMatchedSamplingPoints()
                      << ", xSpanU " << eIter->GetOverlapResult().GetXOverlap().GetXSpanU()
                      << ", xSpanV " << eIter->GetOverlapResult().GetXOverlap().GetXSpanV()
                      << ", xSpanW " << eIter->GetOverlapResult().GetXOverlap().GetXSpanW()
                      << ", xOverlapSpan " << eIter->GetOverlapResult().GetXOverlap().GetXOverlapSpan()
                      << ", chiSquared: " << eIter->GetOverlapResult().GetReducedChi2() << std::endl;

            if (m_showEachIndividualElement)
            {
                if (this->GetPandora().GetGeometry()->GetLArTPC().GetCenterX() > (-370.))
                {
                const ClusterList clusterListU(1, eIter->GetClusterU()), clusterListV(1, eIter->GetClusterV()), clusterListW(1, eIter->GetClusterW());
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListU, "UCluster", RED));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListV, "VCluster", VIOLET));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListW, "WCluster", BLUE));

                MCParticleToIDMap mcParticleToIDMap;
                IDToHitMap idToUHitMap, idToVHitMap, idToWHitMap;
                this->FillMCParticleIDMap(clusterListU.front(), mcParticleToIDMap, idToUHitMap);
                this->FillMCParticleIDMap(clusterListV.front(), mcParticleToIDMap, idToVHitMap);
                this->FillMCParticleIDMap(clusterListW.front(), mcParticleToIDMap, idToWHitMap);                

                std::cout << "UCluster: " << std::endl;
                this->PrintClusterHitOwnershipMap(idToUHitMap);
                std::cout << "VCluster: " << std::endl;
                this->PrintClusterHitOwnershipMap(idToVHitMap);
                std::cout << "WCluster: " << std::endl;
                this->PrintClusterHitOwnershipMap(idToWHitMap);

                std::cout << "COMMON MUON SIZE: " <<  eIter->GetOverlapResult().GetCommonMuonPfoList().size() << std::endl;
                
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                }
            }
        }

        std::cout << " All Connected Clusters " << std::endl;
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterListU, "AllUClusters", RED));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterListV, "AllVClusters", VIOLET));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterListW, "AllWClusters", BLUE));

        if (m_showContext)
        {
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &(pAlgorithm->GetInputClusterList(TPC_VIEW_U)), "InputClusterListU", GRAY));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &(pAlgorithm->GetInputClusterList(TPC_VIEW_V)), "InputClusterListV", GRAY));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &(pAlgorithm->GetInputClusterList(TPC_VIEW_W)), "InputClusterListW", GRAY));
        }

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
   
void DeltaRayTensorVisualizationTool::FillMCParticleIDMap(const Cluster *const pCluster, MCParticleToIDMap &mcParticleToIDMap, IDToHitMap &idToHitMap)
{
    unsigned int particleIDCounter(0);
    for (auto &entry : mcParticleToIDMap)
    {
        if (entry.second > particleIDCounter)
            particleIDCounter = entry.second;
    }

    ++particleIDCounter;

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            try
            {
                const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                const MCParticle *const pLeadingParticle(this->GetLeadingParticle(pHitParticle));

                unsigned int particleID(0);
                MCParticleToIDMap::iterator iter(mcParticleToIDMap.find(pLeadingParticle));

                if (iter == mcParticleToIDMap.end())
                {
                    mcParticleToIDMap[pLeadingParticle] = particleIDCounter;
                    particleID = particleIDCounter;
                    ++particleIDCounter;
                }
                else
                {
                    particleID = mcParticleToIDMap.at(pLeadingParticle);
                }

                idToHitMap[particleID].push_back(pCaloHit);
            }
            catch (StatusCodeException &)
            {
                continue;
            }
        }
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayTensorVisualizationTool::PrintClusterHitOwnershipMap(IDToHitMap &idToHitMap)
{
    std::vector<unsigned int> particleIDVector;
    unsigned int maxHits(0), maxID(0);

    while(particleIDVector.size() != idToHitMap.size())
    {
        for (auto &entry : idToHitMap)
        {
            if (std::find(particleIDVector.begin(), particleIDVector.end(), entry.first) != particleIDVector.end())
                continue;

            if (entry.second.size() > maxHits)
            {
                maxHits = entry.second.size();
                maxID = entry.first;
            }
        }

        particleIDVector.push_back(maxID);
        maxHits = 0;
        maxID = 0;
    }

    unsigned int totalHits(0);
    for (const auto &entry : idToHitMap)
        totalHits += entry.second.size();

    for (const unsigned int id : particleIDVector)
    {
        std::cout << "(" << std::to_string(id) << ", " << static_cast<float>(idToHitMap.at(id).size()) * 100 /static_cast<float>(totalHits) << "%) ";
    }
    std::cout << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
const MCParticle *DeltaRayTensorVisualizationTool::GetLeadingParticle(const MCParticle *const pMCParticle)
{
    if (pMCParticle == LArMCParticleHelper::GetParentMCParticle(pMCParticle))
        return pMCParticle;

    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcParticleVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcParticleVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcParticleVector.push_back(pParentMCParticle);
    }

    return *(mcParticleVector.end() - 2);
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
StatusCode DeltaRayTensorVisualizationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterConnections", m_minClusterConnections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IgnoreUnavailableClusters", m_ignoreUnavailableClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowEachIndividualElement", m_showEachIndividualElement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowContext", m_showContext));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
