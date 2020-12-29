/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMergeTool.cc
 *
 *  @brief  Implementation of the delta ray merge tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayMergeTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

DeltaRayMergeTool::DeltaRayMergeTool() :
    m_maxUnambiguousClusterSeparation(1.f),
    m_maxDRSeparationFromTrack(1.f),
    m_maxVertexSeparation(3.f),
    m_maxClusterSeparation(4.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool mergesMade(false);

    this->MakeMerges(pAlgorithm, overlapTensor, mergesMade);

    return mergesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::MakeMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor, bool &mergesMade) const
{
    bool mergeMade(true);

    while (mergeMade)
    {
        mergeMade = false;

        ClusterVector sortedKeyClusters;
        overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

        ClusterSet usedKeyClusters, modifiedClusters;
        for (const Cluster *const pKeyCluster : sortedKeyClusters)
        {
            if (usedKeyClusters.count(pKeyCluster))
                continue;

            unsigned int nU(0), nV(0), nW(0);
            TensorType::ElementList elementList;
            overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

            if (elementList.size() < 2)
                continue;

            // need to reloop if this happens
            //this->MakeClearMerges(pAlgorithm, elementList, modifiedClusters);

            for (const TensorType::Element &element : elementList)
                usedKeyClusters.insert(element.GetCluster(TPC_VIEW_U));

            const Cluster *pClusterToEnlarge(nullptr), *pClusterToDelete(nullptr);
            while(this->SearchForMerge(elementList, modifiedClusters, pClusterToEnlarge, pClusterToDelete))
            {
                mergeMade = true, mergesMade = true;
                modifiedClusters.insert(pClusterToEnlarge), modifiedClusters.insert(pClusterToDelete);

                ////////////////////////////
                /*
                ClusterList enlarge1({pClusterToEnlarge}), destroy({pClusterToDelete});
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &enlarge1, "pClusterToEnlarge", RED);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &destroy, "pClusterToDelete", BLUE);
                */
                ////////////////////////////
                
                pAlgorithm->UpdateUponDeletion(pClusterToEnlarge); pAlgorithm->UpdateUponDeletion(pClusterToDelete);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                    pAlgorithm->GetClusterListName(LArClusterHelper::GetClusterHitType(pClusterToEnlarge))));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge, pClusterToDelete));
                pAlgorithm->UpdateForNewClusters({pClusterToEnlarge}, {nullptr});

                ////////////////////////////
                /*
                ClusterList enlarge({pClusterToEnlarge});
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &enlarge, "pClusterToEnlargeAFTER", VIOLET);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                */
                ////////////////////////////                
            }
            
            if (this->Jam(pAlgorithm, elementList, modifiedClusters))
            {
                mergeMade = true; mergesMade = true;
            }
            
            // need to take into account simple merges <--- have that return a boolean
            if (modifiedClusters.empty())
                this->PickOutGoodMatches(pAlgorithm, elementList, mergesMade);

            modifiedClusters.clear();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::MakeClearMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, ClusterSet &modifiedClusters) const
{
    ClusterVector uClusters, vClusters, wClusters;
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    
    for (const TensorType::Element &element : elementList)
    {
        for (const HitType &hitType : hitTypeVector)
        {
            ClusterVector &clusterVector(hitType == TPC_VIEW_U ? uClusters : hitType == TPC_VIEW_V ? vClusters : wClusters);

            if (std::find(clusterVector.begin(), clusterVector.end(), element.GetCluster(hitType)) != clusterVector.end())
                continue;
            
            clusterVector.push_back(element.GetCluster(hitType));
        }
    }

    std::sort(uClusters.begin(), uClusters.end(), LArClusterHelper::SortByNHits);
    std::sort(vClusters.begin(), vClusters.end(), LArClusterHelper::SortByNHits);
    std::sort(wClusters.begin(), wClusters.end(), LArClusterHelper::SortByNHits);
    
    this->MakeClearMerges(pAlgorithm, uClusters, modifiedClusters);
    this->MakeClearMerges(pAlgorithm, vClusters, modifiedClusters);
    this->MakeClearMerges(pAlgorithm, wClusters, modifiedClusters);    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::MakeClearMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const ClusterVector &clusterVector, ClusterSet &modifiedClusters) const
{
    for (const Cluster *const pCluster1 : clusterVector)
    {
        for (const Cluster *const pCluster2 : clusterVector)
        {
            if (pCluster1 == pCluster2)
                continue;

            if (modifiedClusters.count(pCluster1) || modifiedClusters.count(pCluster2))
                continue;

            if (LArClusterHelper::GetClosestDistance(pCluster1, pCluster2) < m_maxUnambiguousClusterSeparation)
            {
                modifiedClusters.insert(pCluster1); modifiedClusters.insert(pCluster2);

                pAlgorithm->UpdateUponDeletion(pCluster1); pAlgorithm->UpdateUponDeletion(pCluster2);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                    pAlgorithm->GetClusterListName(LArClusterHelper::GetClusterHitType(pCluster1))));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pCluster1, pCluster2));
                pAlgorithm->UpdateForNewClusters({pCluster1}, {nullptr});
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::SearchForMerge(const TensorType::ElementList &elementList, const ClusterSet &modifiedClusters, const Cluster *&pClusterToEnlarge, const Cluster *&pClusterToDelete) const
{
    HitTypeVector hitTypeVector1({TPC_VIEW_U, TPC_VIEW_V});
    HitTypeVector hitTypeVector2({TPC_VIEW_V, TPC_VIEW_W});
    
    for (const TensorType::Element &element1 : elementList)
    {
        for (const TensorType::Element &element2 : elementList)
        {
            if ((element1.GetCluster(TPC_VIEW_U) == element2.GetCluster(TPC_VIEW_U)) && (element1.GetCluster(TPC_VIEW_V) == element2.GetCluster(TPC_VIEW_V)) &&
                (element1.GetCluster(TPC_VIEW_W) == element2.GetCluster(TPC_VIEW_W)))
            {
                continue;
            }

            for(const HitType &hitType1 : hitTypeVector1)
            {
                if ((element1.GetCluster(hitType1) == element2.GetCluster(hitType1)) && (!modifiedClusters.count(element1.GetCluster(hitType1))))
                {
                    for(const HitType &hitType2 : hitTypeVector2)
                    {
                        if (hitType1 == hitType2)
                            continue;
                        
                        if ((element1.GetCluster(hitType2) == element2.GetCluster(hitType2)) && (!modifiedClusters.count(element1.GetCluster(hitType2))))
                        {
                            const HitType mergeHitType(hitType1 == TPC_VIEW_U ? (hitType2 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_V) : TPC_VIEW_U);

                            pClusterToEnlarge = element1.GetCluster(mergeHitType);
                            pClusterToDelete = element2.GetCluster(mergeHitType);

                            if (modifiedClusters.count(element1.GetCluster(mergeHitType)) || modifiedClusters.count(element2.GetCluster(mergeHitType)))
                                continue;

                            PfoList commonMuonPfoList;
                            this->CombineCommonMuonPfoLists(element1.GetOverlapResult().GetCommonMuonPfoList(), element2.GetOverlapResult().GetCommonMuonPfoList(), commonMuonPfoList);

                            if (commonMuonPfoList.empty())
                                continue;

                            if (this->AreAssociated(commonMuonPfoList, pClusterToEnlarge, pClusterToDelete))
                                return true;

                        }
                    }
                }
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::Jam(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, ClusterSet &modifiedClusters) const
{
    bool mergeMade(true), mergesMade(false);
    HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    
    while(mergeMade)
    {
        mergeMade = false;
    
        for (const TensorType::Element &element1 : elementList)
        {
            for (const TensorType::Element &element2 : elementList)
            {
                if ((element1.GetCluster(TPC_VIEW_U) == element2.GetCluster(TPC_VIEW_U)) && (element1.GetCluster(TPC_VIEW_V) == element2.GetCluster(TPC_VIEW_V)) &&
                    (element1.GetCluster(TPC_VIEW_W) == element2.GetCluster(TPC_VIEW_W)))
                {
                    continue;
                }

                for(const HitType &hitType : hitTypeVector)
                {
                    if ((element1.GetCluster(hitType) == element2.GetCluster(hitType)) && (!modifiedClusters.count(element1.GetCluster(hitType))))
                    {
                        const HitType mergeHitType1(hitType == TPC_VIEW_U ? TPC_VIEW_V : hitType == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
                        const HitType mergeHitType2(mergeHitType1 == TPC_VIEW_U ? TPC_VIEW_V : mergeHitType1 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
                        
                        const Cluster *const pClusterToEnlarge1 = element1.GetCluster(mergeHitType1),  *const pClusterToEnlarge2 = element1.GetCluster(mergeHitType2);
                        const Cluster *const pClusterToDelete1 = element2.GetCluster(mergeHitType1),  *const pClusterToDelete2 = element2.GetCluster(mergeHitType2);

                        if ((pClusterToEnlarge1 == pClusterToDelete1) || (pClusterToEnlarge2 == pClusterToDelete2))
                        {
                            continue;
                        }
                        
                        if (modifiedClusters.count(pClusterToEnlarge1) || modifiedClusters.count(pClusterToEnlarge2) || modifiedClusters.count(pClusterToDelete1) ||
                            modifiedClusters.count(pClusterToDelete2))
                        {
                            continue;
                        }
                        
                        CaloHitList caloHitList1, caloHitList2, caloHitList3;
                        pClusterToEnlarge1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
                        pClusterToDelete1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
                        pClusterToEnlarge2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
                        pClusterToDelete2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
                        element1.GetCluster(hitType)->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);

                        float reducedChiSquared(pAlgorithm->CalculateChiSquared(caloHitList1, caloHitList2, caloHitList3));

                        /*
                        ////////////////////////////
                        ClusterList enlarge({pClusterToEnlarge1, pClusterToEnlarge2}), destroy({pClusterToDelete1, pClusterToDelete2});
                        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &enlarge, "pClusterToEnlarge", RED);
                        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &destroy, "pClusterToDelete", BLUE);
                        std::cout << "element1 chi: " << element1.GetOverlapResult().GetReducedChi2() << std::endl;
                        std::cout << "element2 chi: " << element2.GetOverlapResult().GetReducedChi2() << std::endl;
                        std::cout << "overall chi: " << reducedChiSquared << std::endl;
                        PandoraMonitoringApi::ViewEvent(this->GetPandora());
                        ////////////////////////////
                        */
                        //float reducedChiSquared(pAlgorithm->CalculateChiSquared(caloHitList1, caloHitList2, caloHitList3));

                        if (reducedChiSquared < 1.f)
                        {
                            mergeMade = true; mergesMade = true;
                            
                            modifiedClusters.insert(pClusterToEnlarge1); modifiedClusters.insert(pClusterToEnlarge2);
                            modifiedClusters.insert(pClusterToDelete1); modifiedClusters.insert(pClusterToDelete2);

                            ////////////////////////////
                            //ClusterList enlarge({pClusterToEnlarge1, pClusterToEnlarge2}), destroy({pClusterToDelete1, pClusterToDelete2});
                            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &enlarge, "pClusterToEnlarge", RED);
                            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &destroy, "pClusterToDelete", BLUE);
                            ////////////////////////////

                            pAlgorithm->UpdateUponDeletion(pClusterToEnlarge1); pAlgorithm->UpdateUponDeletion(pClusterToEnlarge2);
                            pAlgorithm->UpdateUponDeletion(pClusterToDelete1); pAlgorithm->UpdateUponDeletion(pClusterToDelete2);

                            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                                pAlgorithm->GetClusterListName(mergeHitType1)));
                        
                            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge1, pClusterToDelete1));

                            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                                pAlgorithm->GetClusterListName(mergeHitType2)));
                        
                            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge2, pClusterToDelete2));

                            ////////////////////////////
                            /*
                            ClusterList enlarge1({pClusterToEnlarge1, pClusterToEnlarge2});
                            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &enlarge1, "pClusterToEnlarge", BLACK);
                            PandoraMonitoringApi::ViewEvent(this->GetPandora());
                            */
                            ////////////////////////////

                            pAlgorithm->UpdateForNewClusters({pClusterToEnlarge1}, {nullptr}); pAlgorithm->UpdateForNewClusters({pClusterToEnlarge2}, {nullptr});
                        }
                    }
                }
            }
        }
    }

    return mergesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::CombineCommonMuonPfoLists(const PfoList &commonMuonPfoList1, const PfoList &commonMuonPfoList2, PfoList &commonMuonPfoList) const
{
    for (const ParticleFlowObject *const pCommonMuonPfo1 : commonMuonPfoList1)
    {
        for (const ParticleFlowObject *const pCommonMuonPfo2 : commonMuonPfoList2)
        {
            if (pCommonMuonPfo1 == pCommonMuonPfo2)
                commonMuonPfoList.push_back(pCommonMuonPfo1);
        }
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::AreAssociated(const PfoList &commonMuonPfoList, const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete) const
{
    /*
    ClusterList enlarge({pClusterToEnlarge}), destroy({pClusterToDelete});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &enlarge, "pClusterToEnlarge", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &destroy, "pClusterToDelete", BLUE);
    */
    
    bool areAttached(false);
    if (this->IsHiddenTrack(commonMuonPfoList, pClusterToEnlarge, pClusterToDelete, areAttached))
    {
        //std::cout << "vertices close enough" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return true;
    } 
    /*
    if (areAttached)
    {
        std::cout << "vertices not close enough" << std::endl;
    }
    */
    if (!areAttached)
    {
        //std::cout << "not attached" << std::endl;
        if (this->IsBrokenCluster(pClusterToEnlarge, pClusterToDelete))
        {
            //std::cout << "is broken cluster" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());            
            return true;
        }

        //std::cout << "clusters not close enough" << std::endl;
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    return false;
}
   
//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::IsHiddenTrack(const PfoList &commonMuonPfoList, const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete, bool &areAttached) const
{
    for (const ParticleFlowObject *const pCommonMuonPfo : commonMuonPfoList)
    {
        if (this->IsConnected(pCommonMuonPfo, pClusterToEnlarge) && this->IsConnected(pCommonMuonPfo, pClusterToDelete))
        {
            //std::cout << "are connected to track" << std::endl;

            CaloHitList enlargeVertices, deleteVertices;
            this->FindVertices(pCommonMuonPfo, pClusterToEnlarge, enlargeVertices);
            this->FindVertices(pCommonMuonPfo, pClusterToDelete, deleteVertices);            

            float closestDistance(std::numeric_limits<float>::max());
            for (const CaloHit *const pCaloHit : enlargeVertices)
            {
                const float separation(this->GetClosestDistance(pCaloHit, deleteVertices));
                if (separation < closestDistance)
                    closestDistance = separation;
            }

            //std::cout << "vertex separation: " << closestDistance  << std::endl;
            
            areAttached = true;
            
            if (closestDistance < m_maxVertexSeparation)
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::FindVertices(const Pfo *const pCommonMuonPfo, const Cluster *const pCluster, CaloHitList &vertexList) const
{
    HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pCommonMuonPfo, hitType, muonClusterList);    

    if (muonClusterList.size() != 1)
    {
        std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
        return;
    }
    
    CaloHitList caloHitList;
    muonClusterList.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pCluster) < m_maxDRSeparationFromTrack)
            vertexList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool DeltaRayMergeTool::IsConnected(const Pfo *const pCommonMuonPfo, const Cluster *const pCluster) const
{
    HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    
    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pCommonMuonPfo, hitType, muonClusterList);    

    if (muonClusterList.size() != 1)
    {
        std::cout << "ISOBEL SIZE DOES NOT EQUAL ONE" << std::endl;
        return false;
    }
    
    const float separation(LArClusterHelper::GetClosestDistance(pCluster, muonClusterList));

    return separation < m_maxDRSeparationFromTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::IsBrokenCluster(const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete) const
{
    const float clusterSeparation(LArClusterHelper::GetClosestDistance(pClusterToEnlarge, pClusterToDelete));

    // also demand clusters on same side or the delta ray?
    
    return clusterSeparation < m_maxClusterSeparation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMergeTool::GetClosestDistance(const CaloHit *const pCaloHit, const CaloHitList &caloHitList) const
{
    float shortestDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector referencePoint(pCaloHit->GetPositionVector());

    for (const CaloHit *const pTestCaloHit : caloHitList)
    {
        const CartesianVector &position(pTestCaloHit->GetPositionVector());
        float separationSquared((position - referencePoint).GetMagnitude());

        if (separationSquared < shortestDistanceSquared)
            shortestDistanceSquared = separationSquared;
    }

    return shortestDistanceSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::PickOutGoodMatches(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, bool &particlesCreated) const
{
    ProtoParticleVector protoParticleVector;

    bool found(true);
    ClusterSet usedClusters;
    
    while (found)
    {
        found = false;

        float highestHitCount(-std::numeric_limits<float>::max());//, bestChi(0);
        const Cluster *pBestClusterU(nullptr), *pBestClusterV(nullptr), *pBestClusterW(nullptr);

        for (const TensorType::Element &element : elementList)
        {
            if (element.GetOverlapResult().GetReducedChi2() > 1.f)
                continue;
            
            const Cluster *const pClusterU(element.GetCluster(TPC_VIEW_U)), *const pClusterV(element.GetCluster(TPC_VIEW_V)), *const pClusterW(element.GetCluster(TPC_VIEW_W));
            
            if (usedClusters.count(pClusterU) || usedClusters.count(pClusterV) || usedClusters.count(pClusterW))
                continue;
            
            if ((pClusterU->GetNCaloHits() + pClusterV->GetNCaloHits() + pClusterW->GetNCaloHits()) > highestHitCount)
            {
                //bestChi = element.GetOverlapResult().GetReducedChi2();
                highestHitCount = pClusterU->GetNCaloHits() + pClusterV->GetNCaloHits() + pClusterW->GetNCaloHits();
                pBestClusterU = pClusterU; pBestClusterV = pClusterV; pBestClusterW = pClusterW;
            }
        }

        if (pBestClusterU && pBestClusterV && pBestClusterW)
        {
            found = true;
            usedClusters.insert(pBestClusterU); usedClusters.insert(pBestClusterV); usedClusters.insert(pBestClusterW);

            /*
            std::cout << "WILL CREATE THIS PARTICLE" << std::endl;
            std::cout << "MATCH CHI SQUARED: " << bestChi << std::endl;
            ClusterList u({pBestClusterU}), v({pBestClusterV}), w({pBestClusterW});
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &u, "u", BLUE);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &v, "v", BLUE);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &w, "w", BLUE);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            */
            
            ProtoParticle protoParticle;
            protoParticle.m_clusterList.push_back(pBestClusterU);
            protoParticle.m_clusterList.push_back(pBestClusterV);
            protoParticle.m_clusterList.push_back(pBestClusterW);
            protoParticleVector.push_back(protoParticle);
        }
    }                

    if (!protoParticleVector.empty())
    {
        pAlgorithm->CreatePfos(protoParticleVector);
        particlesCreated = true;
    }
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayMergeTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDRSeparationFromTrack", m_maxDRSeparationFromTrack));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexSeparation", m_maxVertexSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterSeparation", m_maxClusterSeparation));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxUnambiguousClusterSeparation", m_maxUnambiguousClusterSeparation));            
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
