/**
 *  @file   larpandoracontent/LArMonitoring/DeltaRayMatrixVisualizationTool.cc
 *
 *  @brief  Implementation of the delta ray tensor visualization tool class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMonitoring/DeltaRayMatrixVisualizationTool.h"

using namespace pandora;

namespace lar_content
{

DeltaRayMatrixVisualizationTool::DeltaRayMatrixVisualizationTool() :
    m_minClusterConnections(1),
    m_ignoreUnavailableClusters(true),
    m_showEachIndividualElement(false),
    m_showContext(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMatrixVisualizationTool::Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterSet usedKeyClusters;
    ClusterVector sortedKeyClusters;
    overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (m_ignoreUnavailableClusters && !pKeyCluster->IsAvailable())
            continue;

        if (usedKeyClusters.count(pKeyCluster))
            continue;

        unsigned int n1(0), n2(0);
        MatrixType::ElementList elementList;
        overlapMatrix.GetConnectedElements(pKeyCluster, m_ignoreUnavailableClusters, elementList, n1, n2);

        if ((n1 < m_minClusterConnections) && (n1 < m_minClusterConnections))
            continue;

        if (n1 * n2 == 0)
            continue;

        int counter(0);
        ClusterList allClusterList1, allClusterList2;
        std::cout << " Connections: n1 " << n1 << ", n2 " << n2 << ", nElements " << elementList.size() << std::endl;

        for (MatrixType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (allClusterList1.end() == std::find(allClusterList1.begin(), allClusterList1.end(), eIter->GetCluster1())) allClusterList1.push_back(eIter->GetCluster1());
            if (allClusterList2.end() == std::find(allClusterList2.begin(), allClusterList2.end(), eIter->GetCluster2())) allClusterList2.push_back(eIter->GetCluster2());
            usedKeyClusters.insert(eIter->GetCluster1());

            std::cout << " Element " << counter++ << ": MatchedFraction " << eIter->GetOverlapResult().GetMatchedFraction()
                      << ", MatchedSamplingPoints " << eIter->GetOverlapResult().GetNMatchedSamplingPoints()
                      << ", xSpan1 " << eIter->GetOverlapResult().GetXOverlap().GetXSpan0()
                      << ", xSpan2 " << eIter->GetOverlapResult().GetXOverlap().GetXSpan1()
                      << ", xOverlapSpan " << eIter->GetOverlapResult().GetXOverlap().GetTwoViewXOverlapSpan()
                      << ", xOverlapFraction1 " << eIter->GetOverlapResult().GetXOverlap().GetTwoViewXOverlapSpan() / eIter->GetOverlapResult().GetXOverlap().GetXSpan0()
                      << ", xOverlapFraction2 " << eIter->GetOverlapResult().GetXOverlap().GetTwoViewXOverlapSpan() / eIter->GetOverlapResult().GetXOverlap().GetXSpan1()
                      << ", chiSquared: " << eIter->GetOverlapResult().GetReducedChi2() << std::endl;

            if (m_showEachIndividualElement)
            {
                if (this->GetPandora().GetGeometry()->GetLArTPC().GetCenterX() > -170)//(-370.))
                {
                    const ClusterList clusterList1(1, eIter->GetCluster1()), clusterList2(1, eIter->GetCluster2());
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList1, "1Cluster", VIOLET));
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList2, "2Cluster", BLUE));

                std::cout << "HIT COLLECTION VECTOR SIZE: " << eIter->GetOverlapResult().GetProjectedCaloHits().size() << std::endl;
                for (const CaloHit *const pCaloHit : eIter->GetOverlapResult().GetProjectedCaloHits())
                {
                    const CartesianVector &position(pCaloHit->GetPositionVector());
                    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "PROJECTED HIT", RED, 2));
                }

                MCParticleToIDMap mcParticleToIDMap;
                IDToHitMap idTo1HitMap, idTo2HitMap;
                this->FillMCParticleIDMap(clusterList1.front(), mcParticleToIDMap, idTo1HitMap);
                this->FillMCParticleIDMap(clusterList2.front(), mcParticleToIDMap, idTo2HitMap);

                std::cout << "1Cluster: " << std::endl;
                this->PrintClusterHitOwnershipMap(idTo1HitMap);
                std::cout << "2Cluster: " << std::endl;
                this->PrintClusterHitOwnershipMap(idTo2HitMap);
                
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                }
            }
        }

        std::cout << " All Connected Clusters " << std::endl;
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterList1, "All1Clusters", VIOLET));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &allClusterList2, "All2Clusters", BLUE));

        if (m_showContext)
        {
            std::cout << "ISOBEL: I COULDN'T BE BOTHERED TO DO THIS - SORRY FUTURE SELF" << std::endl; 
            //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &(pAlgorithm->GetInputClusterList(TPC_VIEW_U)), "InputClusterList1", GRAY));
            //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &(pAlgorithm->GetInputClusterList(TPC_VIEW_V)), "InputClusterList2", GRAY));
        }

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
   
void DeltaRayMatrixVisualizationTool::FillMCParticleIDMap(const Cluster *const pCluster, MCParticleToIDMap &mcParticleToIDMap, IDToHitMap &idToHitMap)
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

void DeltaRayMatrixVisualizationTool::PrintClusterHitOwnershipMap(IDToHitMap &idToHitMap)
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
    
const MCParticle *DeltaRayMatrixVisualizationTool::GetLeadingParticle(const MCParticle *const pMCParticle)
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
    
StatusCode DeltaRayMatrixVisualizationTool::ReadSettings(const TiXmlHandle xmlHandle)
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
