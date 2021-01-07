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
    m_maxDRSeparationFromTrack(1.5f),
    m_maxVertexSeparation(10.f),
    m_maxClusterSeparation(3.f)
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
        std::cout << "LOOP BEGINS" << std::endl;
        mergeMade = false;

        ClusterVector sortedKeyClusters;
        overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

        ClusterSet usedKeyClusters, modifiedClusters;
        for (const Cluster *const pKeyCluster : sortedKeyClusters)
        {
            if (usedKeyClusters.count(pKeyCluster))
                continue;

            if (modifiedClusters.count(pKeyCluster))
                continue;

            unsigned int nU(0), nV(0), nW(0);
            TensorType::ElementList elementList;
            overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

            for (const TensorType::Element &element : elementList)
	    {
		if (usedKeyClusters.count(element.GetClusterU()))
		    continue;

		usedKeyClusters.insert(element.GetClusterU());
	    }

            if (elementList.size() < 2)
                continue;

            if (this->MakeTwoCommonViewMerges(pAlgorithm, elementList, modifiedClusters))
	    {
	        std::cout << "MADE COMMON TWO VIEW MERGES REPEATING" << std::endl;
                mergeMade = true; mergesMade = true;
		break;
	    }

	    if (this->MakeOneCommonViewMerges(pAlgorithm, elementList, modifiedClusters))
	    {
	        std::cout << "MADE COMMON ONE VIEW MERGES REPEATING" << std::endl;
                mergeMade = true; mergesMade = true;
		break;
	    }

	    this->PickOutGoodMatches(pAlgorithm, elementList, modifiedClusters);

	    if (modifiedClusters.size())
	      break;
	}

	modifiedClusters.clear();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::MakeTwoCommonViewMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, ClusterSet &modifiedClusters) const
{
    std::cout << "making two view common merges..." << std::endl;
    const HitTypeVector hitTypeVector1({TPC_VIEW_U, TPC_VIEW_V});
    const HitTypeVector hitTypeVector2({TPC_VIEW_V, TPC_VIEW_W});

    bool mergeMade(true), mergesMade(false);
    
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

                for (const HitType &hitType1 : hitTypeVector1)
                {
		    if ((element1.GetCluster(hitType1) == element2.GetCluster(hitType1))) //&& (!modifiedClusters.count(element1.GetCluster(hitType1))))
                    {
		        if (modifiedClusters.count(element1.GetCluster(hitType1)))
			  continue;

                        for(const HitType &hitType2 : hitTypeVector2)
                        {
                            if (hitType1 == hitType2)
                                continue;
                        
                            if ((element1.GetCluster(hitType2) == element2.GetCluster(hitType2))) //&& (!modifiedClusters.count(element1.GetCluster(hitType2))))
                            {
			        if (modifiedClusters.count(element1.GetCluster(hitType2)))
			            continue;

                                const HitType mergeHitType(hitType1 == TPC_VIEW_U ? (hitType2 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_V) : TPC_VIEW_U);

                                const Cluster *pClusterToEnlarge(element1.GetCluster(mergeHitType)), *pClusterToDelete(element2.GetCluster(mergeHitType));

                                if (modifiedClusters.count(pClusterToEnlarge) || modifiedClusters.count(pClusterToDelete))
                                    continue;
				
				ClusterList e({pClusterToEnlarge}), d({pClusterToDelete});
				PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &e, "enlarge", RED);
				PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &d, "delete", BLUE);
				

                                if (this->AreAssociated(element1, element2, mergeHitType))
                                {
				  std::cout << "making merge" << std::endl;
				  PandoraMonitoringApi::ViewEvent(this->GetPandora());

				    mergesMade = true;
                                    
                                    modifiedClusters.insert(pClusterToEnlarge), modifiedClusters.insert(pClusterToDelete);
                                    pAlgorithm->UpdateUponDeletion(pClusterToEnlarge); pAlgorithm->UpdateUponDeletion(pClusterToDelete);
                                    
                                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                                        pAlgorithm->GetClusterListName(LArClusterHelper::GetClusterHitType(pClusterToEnlarge))));
                                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge, pClusterToDelete));
                                    
                                    pAlgorithm->UpdateForNewClusters({pClusterToEnlarge}, {nullptr});
                                }
				else
				{
				 std::cout << "not associated" << std::endl;
				  PandoraMonitoringApi::ViewEvent(this->GetPandora());
				  }
                            }
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

bool DeltaRayMergeTool::AreAssociated(const TensorType::Element &element1, const TensorType::Element &element2, const HitType &mergeHitType) const
{
  // Demand the elements to have a shared common muon
  PfoList commonMuonPfoList;
  this->CombineCommonMuonPfoLists(element1.GetOverlapResult().GetCommonMuonPfoList(), element2.GetOverlapResult().GetCommonMuonPfoList(), commonMuonPfoList);

  if (commonMuonPfoList.empty())
  {
    std::cout << "no common muons" << std::endl;
    return false;
  }

  const Cluster *const pCluster1(element1.GetCluster(mergeHitType)), *const pCluster2(element2.GetCluster(mergeHitType));

  PfoList connectedMuonPfoList1, connectedMuonPfoList2;
  this->GetConnectedMuons(element1.GetOverlapResult().GetCommonMuonPfoList(), pCluster1, connectedMuonPfoList1);
  this->GetConnectedMuons(element2.GetOverlapResult().GetCommonMuonPfoList(), pCluster2, connectedMuonPfoList2);

  if (connectedMuonPfoList1.empty() || connectedMuonPfoList2.empty())
  {
    if (this->IsBrokenCluster(pCluster1, pCluster2))                                                                                                                                   
    {
      std::cout << "broken cluster"<< std::endl;
      return true;       
    }
  }

  for (const ParticleFlowObject *const pConnectedMuon1 : connectedMuonPfoList1)
  {
    for (const ParticleFlowObject *const pConnectedMuon2 : connectedMuonPfoList2)
    {
      if (pConnectedMuon1 == pConnectedMuon2)
      {
	if (this->IsHiddenTrack(pConnectedMuon1, pCluster1, pCluster2))
	{
	  std::cout << "is hidden track" << std::endl;
	  return true;
	}
      }
    }
  }

  std::cout << "connected to different muons" << std::endl;

  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// of the element not the combined one!
void DeltaRayMergeTool::GetConnectedMuons(const PfoList &commonMuonPfoList, const Cluster *const pClusterToEnlarge, PfoList &connectedMuonPfoList) const 
{
  for (const ParticleFlowObject *const pCommonMuonPfo : commonMuonPfoList)
  {
    if (this->IsConnected(pCommonMuonPfo, pClusterToEnlarge))
      connectedMuonPfoList.push_back(pCommonMuonPfo);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::IsHiddenTrack(const ParticleFlowObject *const pMuonPfo, const Cluster *const pCluster1, const Cluster *const pCluster2) const 
{
  CaloHitList vertices1, vertices2;
  this->FindVertices(pMuonPfo, pCluster1, vertices1);
  this->FindVertices(pMuonPfo, pCluster2, vertices2);   

  float closestDistance(std::numeric_limits<float>::max());
  for (const CaloHit *const pCaloHit : vertices1)
  {
    const float separation(this->GetClosestDistance(pCaloHit, vertices2));
    if (separation < closestDistance)
      closestDistance = separation;
  }

  return (closestDistance < m_maxVertexSeparation);
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


bool DeltaRayMergeTool::IsBrokenCluster(const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete) const
{
    const float clusterSeparation(LArClusterHelper::GetClosestDistance(pClusterToEnlarge, pClusterToDelete));

    // also demand clusters on same side or the delta ray?
    
    return clusterSeparation < m_maxClusterSeparation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMergeTool::MakeOneCommonViewMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, ClusterSet &modifiedClusters) const
{
    std::cout << "making one view common merges..." << std::endl;
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
    
    bool mergeMade(true), mergesMade(false);
    
    while (mergeMade)
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

                for (const HitType &hitType : hitTypeVector)
                {
                    if ((element1.GetCluster(hitType) == element2.GetCluster(hitType)))// && (!modifiedClusters.count(element1.GetCluster(hitType))))
                    {
                        if (modifiedClusters.count(element1.GetCluster(hitType)))
                            continue;


                        const HitType mergeHitType1(hitType == TPC_VIEW_U ? TPC_VIEW_V : hitType == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
                        const HitType mergeHitType2(mergeHitType1 == TPC_VIEW_U ? TPC_VIEW_V : mergeHitType1 == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
                        
                        const Cluster *pClusterToEnlarge1 = element1.GetCluster(mergeHitType1),  *pClusterToEnlarge2 = element1.GetCluster(mergeHitType2);
                        const Cluster *pClusterToDelete1 = element2.GetCluster(mergeHitType1),  *pClusterToDelete2 = element2.GetCluster(mergeHitType2);

			ClusterList e1({pClusterToEnlarge1}), d1({pClusterToDelete1}), e2({pClusterToEnlarge2}), d2({pClusterToDelete2});
			PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &e1, "enlarge1", RED);
			PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &d1, "delete1", BLUE);
			PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &e2, "enlarge2", RED);
			PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &d2, "delete2", BLUE);


                        if ((pClusterToEnlarge1 == pClusterToDelete1) || (pClusterToEnlarge2 == pClusterToDelete2))
                        {
			  std::cout << "one of the clusters is the same " << std::endl;
			  PandoraMonitoringApi::ViewEvent(this->GetPandora());
                            continue;
                        }
                        
                        if (modifiedClusters.count(pClusterToEnlarge1) || modifiedClusters.count(pClusterToEnlarge2) || modifiedClusters.count(pClusterToDelete1) || modifiedClusters.count(pClusterToDelete2))
                        {
			  std::cout << "in modified clusters" << std::endl;
			  PandoraMonitoringApi::ViewEvent(this->GetPandora());
                            continue;
                        }

                        if (!this->AreAssociated(element1, element2, mergeHitType1))
			{
			  std::cout << "first clusters are not associated" << std::endl;
			  PandoraMonitoringApi::ViewEvent(this->GetPandora());
			  continue;
			}

                        if (!this->AreAssociated(element1, element2, mergeHitType2))
			{
			  std::cout << "second clusters are not associated" << std::endl;
			  PandoraMonitoringApi::ViewEvent(this->GetPandora());
			  continue;
			}
                        
                        CaloHitList caloHitList1, caloHitList2, caloHitList3;
                        pClusterToEnlarge1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
                        pClusterToDelete1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
                        pClusterToEnlarge2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
                        pClusterToDelete2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
                        element1.GetCluster(hitType)->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);

                        XOverlap xOverlapObject(0.f,0.f,0.f,0.f,0.f, 0.f,0.f);
                        float chiSquaredSum(0.f);
                        unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
                        
                        StatusCode status(pAlgorithm->PerformMatching(caloHitList1, caloHitList2, caloHitList3, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xOverlapObject));
                        
                        if (status == STATUS_CODE_NOT_FOUND)
			{
			  std::cout << "matching didn't work" << std::endl;
			  PandoraMonitoringApi::ViewEvent(this->GetPandora());
                            continue;
			}

                        if (status != STATUS_CODE_SUCCESS)
                        {
                            std::cout << "THIS IS WHERE" << std::endl;
                            throw status;
                        }

                        float reducedChiSquared(chiSquaredSum / nSamplingPoints);
                        
                        if (reducedChiSquared < 1.f)
                        {
			  std::cout << "will merge" << std::endl;
			  PandoraMonitoringApi::ViewEvent(this->GetPandora());
			  //std::cout << "NOW" << std::endl;
			  //std::cout << "pClusterToEnlarge1: " << pClusterToEnlarge1->GetNCaloHits() << std::endl;
			  //std::cout << "pClusterToEnlarge2: " << pClusterToEnlarge2->GetNCaloHits() << std::endl;
			  //std::cout << "pClusterToDelete1: " << pClusterToDelete1->GetNCaloHits() << std::endl;
			  //std::cout << "pClusterToDelete2: " << pClusterToDelete2->GetNCaloHits() << std::endl;
			  //mergeMade = true; 
                            mergesMade = true;
                            
                            modifiedClusters.insert(pClusterToEnlarge1); modifiedClusters.insert(pClusterToEnlarge2);
                            modifiedClusters.insert(pClusterToDelete1); modifiedClusters.insert(pClusterToDelete2);

                            pAlgorithm->UpdateUponDeletion(pClusterToEnlarge1); pAlgorithm->UpdateUponDeletion(pClusterToDelete1);  

                            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                                pAlgorithm->GetClusterListName(mergeHitType1)));
                        
                            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge1, pClusterToDelete1));

                            pAlgorithm->UpdateForNewClusters({pClusterToEnlarge1}, {nullptr});

                            pAlgorithm->UpdateUponDeletion(pClusterToEnlarge2); pAlgorithm->UpdateUponDeletion(pClusterToDelete2);

                            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
                                pAlgorithm->GetClusterListName(mergeHitType2)));
                        
                            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pClusterToEnlarge2, pClusterToDelete2));
			   
			    //ClusterVector clusterVector; PfoVector pfoVector;
			    //clusterVector.push_back(pClusterToEnlarge1); clusterVector.push_back(pClusterToEnlarge2);
			    //pfoVector.push_back(nullptr); pfoVector.push_back(nullptr);

                            //pAlgorithm->UpdateForNewClusters(clusterVector, pfoVector);
			    

			    pAlgorithm->UpdateForNewClusters({pClusterToEnlarge2}, {nullptr});

			    //std::cout << "pClusterToEnlarge1AFTER: " << pClusterToEnlarge1->GetNCaloHits() << std::endl;
			    //std::cout << "pClusterToEnlarge2AFTER: " << pClusterToEnlarge2->GetNCaloHits() << std::endl;
                        }
			else
			{
			  std::cout << "chi squared not good enough" << std::endl;
			  PandoraMonitoringApi::ViewEvent(this->GetPandora());
			}
                    }
                }
            }
        }
    }
    return mergesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMergeTool::PickOutGoodMatches(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, ClusterSet &modifiedClusters) const
{
    ProtoParticleVector protoParticleVector;
    
    bool found(true);
   
    while (found)
    {
        found = false;

        float highestHitCount(-std::numeric_limits<float>::max()), bestChiSquared(0.f);
        const Cluster *pBestClusterU(nullptr), *pBestClusterV(nullptr), *pBestClusterW(nullptr);

        for (const TensorType::Element &element : elementList)
        {
            if (element.GetOverlapResult().GetReducedChi2() > 1.f)
                continue;
            
            const Cluster *const pClusterU(element.GetCluster(TPC_VIEW_U)), *const pClusterV(element.GetCluster(TPC_VIEW_V)), *const pClusterW(element.GetCluster(TPC_VIEW_W));
            
            if (modifiedClusters.count(pClusterU) || modifiedClusters.count(pClusterV) || modifiedClusters.count(pClusterW))
                continue;

            const float chiSquared = element.GetOverlapResult().GetReducedChi2();            
            const unsigned int hitSum(pClusterU->GetNCaloHits() + pClusterV->GetNCaloHits() + pClusterW->GetNCaloHits());

            if ((hitSum == highestHitCount) && (chiSquared < bestChiSquared))
            {
                bestChiSquared = chiSquared;
                highestHitCount = hitSum;
                pBestClusterU = pClusterU; pBestClusterV = pClusterV; pBestClusterW = pClusterW;
                
                continue;
            }
            
            if (hitSum > highestHitCount)
            {
                bestChiSquared = chiSquared;
                highestHitCount = hitSum;
                pBestClusterU = pClusterU; pBestClusterV = pClusterV; pBestClusterW = pClusterW;
            }
        }

        if (pBestClusterU && pBestClusterV && pBestClusterW)
        {
            found = true;
            modifiedClusters.insert(pBestClusterU); modifiedClusters.insert(pBestClusterV); modifiedClusters.insert(pBestClusterW);
            
            ProtoParticle protoParticle;
            protoParticle.m_clusterList.push_back(pBestClusterU);
            protoParticle.m_clusterList.push_back(pBestClusterV);
            protoParticle.m_clusterList.push_back(pBestClusterW);
            protoParticleVector.push_back(protoParticle);
        }
    }        

    if (!protoParticleVector.empty())
        pAlgorithm->CreatePfos(protoParticleVector);
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
