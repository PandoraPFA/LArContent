/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMergeTool.cc
 *
 *  @brief  Implementation of the delta ray merge tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMergeTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

TwoViewDeltaRayMergeTool::TwoViewDeltaRayMergeTool() :
    m_maxUnambiguousClusterSeparation(1.f),
    m_maxDRSeparationFromTrack(1.f),
    m_maxVertexSeparation(3.f),
    m_maxClusterSeparation(4.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMergeTool::Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool mergesMade(false);

    this->MakeMerges(pAlgorithm, overlapMatrix, mergesMade);

    return mergesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMergeTool::MakeMerges(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix, bool &mergesMade) const
{
    bool mergeMade(true);

    while (mergeMade)
    {
        mergeMade = false;

        ClusterVector sortedKeyClusters;
        overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

        ClusterSet usedKeyClusters;
        for (const Cluster *const pKeyCluster : sortedKeyClusters)
        {
            if (usedKeyClusters.count(pKeyCluster))
                continue;

            unsigned int n1(0), n2(0);
            MatrixType::ElementList elementList;
            overlapMatrix.GetConnectedElements(pKeyCluster, true, elementList, n1, n2);

            for (const MatrixType::Element &element : elementList)
	        {
                if (usedKeyClusters.count(element.GetCluster1()))
                    continue;

                usedKeyClusters.insert(element.GetCluster1());
            }

            if (this->PickOutGoodMatches(pAlgorithm, elementList))
            {
                mergeMade = true;
                mergesMade = true;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMergeTool::PickOutGoodMatches(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList) const
{
    bool found(false);   
        
    float highestHitCount(-std::numeric_limits<float>::max()), bestChiSquared(0.f);
    MatrixType::Element bestElement(nullptr, nullptr, TrackTwoViewTopologyOverlapResult(TwoViewXOverlap(0.f, 0.f, 0.f, 0.f), PfoList(), nullptr, ClusterList(), 0.f));
    
    for (const MatrixType::Element &element : elementList)
    {            
        const Cluster *const pCluster1(element.GetCluster1()), *const pCluster2(element.GetCluster2()), *const pCluster3(element.GetOverlapResult().GetBestMatchedCluster());

        //ATTN: Best matched cluster can be removed during pfo creation process and may not be replaceable
        if (!pCluster3)
            continue;

        const float chiSquared = element.GetOverlapResult().GetReducedChiSquared();
        const unsigned int hitSum(pCluster1->GetNCaloHits() + pCluster2->GetNCaloHits() + (pCluster3 ? pCluster3->GetNCaloHits() : 0));

        if ((hitSum == highestHitCount) && (chiSquared < bestChiSquared))
        {
            found = true;
            bestChiSquared = chiSquared;
            highestHitCount = hitSum;
            bestElement = element;
                
            continue;
        }
            
        if (hitSum > highestHitCount)
        {
            found = true;
            bestChiSquared = chiSquared;
            highestHitCount = hitSum;
            bestElement = element;
        }
    }
    
    if (found)
    {
        this->CreatePfo(pAlgorithm, bestElement);
	return true;
    }
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool TwoViewDeltaRayMergeTool::CreatePfo(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element) const
{
    //std::cout << "CREATING PFO" << std::endl;
    
    ProtoParticle protoParticle;
    protoParticle.m_clusterList.push_back(element.GetCluster1());
    protoParticle.m_clusterList.push_back(element.GetCluster2());

    const Cluster *const pBestMatchedCluster(element.GetOverlapResult().GetBestMatchedCluster());
    
    //std::cout << "best matched cluster hit number: " << pBestMatchedCluster->GetNCaloHits() << std::endl;

    if (pBestMatchedCluster)
    {
        this->GrowThirdView(pAlgorithm, element, protoParticle);
    }
    else
    {
        std::cout << "TwoViewCreation" << std::endl;
	/*  
        ClusterList c1({element.GetCluster1()}), c2({element.GetCluster2()});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &c1, "ClusterList1", BLUE);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &c2, "ClusterList2", RED);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
    }

    ProtoParticleVector protoParticleVector({protoParticle});

    return (pAlgorithm->CreateThreeDParticles(protoParticleVector));
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMergeTool::GrowThirdView(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, ProtoParticle &protoParticle) const
{
    const PfoList &commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());
    const Cluster *const pBestMatchedCluster(element.GetOverlapResult().GetBestMatchedCluster());
    const HitType &thirdViewHitType(LArClusterHelper::GetClusterHitType(pBestMatchedCluster));
    const ParticleFlowObject *pMatchedMuonPfo(nullptr);
    
    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
    {
        ClusterList muonClusterList;
        LArPfoHelper::GetClusters(pMuonPfo, thirdViewHitType, muonClusterList);

        if (std::find(muonClusterList.begin(), muonClusterList.end(), pBestMatchedCluster) != muonClusterList.end())
            pMatchedMuonPfo = pMuonPfo;
    }
    
    /*
    ClusterList c1({element.GetCluster1()}), c2({element.GetCluster2()}), c3({pBestMatchedCluster}), bestMatchedClusterList(element.GetOverlapResult().GetMatchedClusterList());
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &c1, "ClusterList1", BLUE);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &c2, "ClusterList2", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &c3, "BestMatchedCluster", BLACK);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &bestMatchedClusterList, "BestMatchedClusterList", VIOLET);

    std::cout << "bestMatchedClusterList size: " << bestMatchedClusterList.size() << std::endl;
    
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    
    if (pMatchedMuonPfo)
    {
      //std::cout << "BEST MATCHED IS A MUON" << std::endl;
        
        CartesianPointVector deltaRayProjectedPositions;
        this->ProjectPositions(element.GetCluster1(), element.GetCluster2(), deltaRayProjectedPositions);

        /*
        for (const CartesianVector &deltaRayPosition : deltaRayProjectedPositions)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &deltaRayPosition, "Projected Position", BLACK, 2);
        */
        CaloHitList deltaRayHits;
        if ((this->CollectDeltaRayHits(element, deltaRayProjectedPositions, pMatchedMuonPfo, deltaRayHits) != STATUS_CODE_SUCCESS) || (deltaRayHits.empty()))
        {
	  //std::cout << "couldn't pull oute a DR cluster" << std::endl;
            const Cluster *pSeedCluster(nullptr);
            this->GetBestMatchedAvailableCluster(element.GetOverlapResult().GetMatchedClusterList(), pSeedCluster);
            
            if (pSeedCluster)
            {
	      //std::cout << "NO SEED FOUND" << std::endl;
                //ClusterList seedBefore({pSeedCluster});
                //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &seedBefore, "seedBefore", BLUE);
                
                this->MergeThirdView(pAlgorithm, element, pSeedCluster);

                //ClusterList seedAfter({pSeedCluster});
                //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &seedAfter, "seedAfter", RED);

                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                
	        pAlgorithm->RemoveThirdViewCluster(pSeedCluster);                
                protoParticle.m_clusterList.push_back(pSeedCluster);
            }
            //else
	      //{
	      //std::cout << "no seed cluster found" << std::endl;
	      //}
        }
        else
        {
            /*
            for (const CaloHit *const pCaloHit : deltaRayHits)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "CollectedHits", GREEN, 2);
            }
            */

	  //std::cout << "oulled out DR cluster" << std::endl;

            const Cluster *const pSeedCluster(this->SplitCluster(pAlgorithm, pBestMatchedCluster, deltaRayHits));
            
            //ClusterList seedBefore({pSeedCluster});
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &seedBefore, "seedBefore", BLUE);
            
            this->MergeThirdView(pAlgorithm, element, pSeedCluster);

            //ClusterList seedAfter({pSeedCluster});
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &seedAfter, "seedAfter", RED);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
	    pAlgorithm->RemoveThirdViewCluster(pSeedCluster);            
            protoParticle.m_clusterList.push_back(pSeedCluster);
        }   
    }
    else
    {
      //std::cout << "BEST MATCHED IS AN AVAILBALE DELTA RAY" << std::endl;

        //ClusterList seedBefore({pBestMatchedCluster});
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &seedBefore, "seedBefore", BLUE);
        
        this->MergeThirdView(pAlgorithm, element, pBestMatchedCluster);

        //ClusterList seedAfter({pBestMatchedCluster});
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &seedAfter, "seedAfter", RED);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());        
	pAlgorithm->RemoveThirdViewCluster(pBestMatchedCluster);        
        protoParticle.m_clusterList.push_back(pBestMatchedCluster);
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void TwoViewDeltaRayMergeTool::MergeThirdView(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const Cluster *const pSeedCluster) const
{
    const Cluster *const pCluster1(element.GetCluster1()), *const pCluster2(element.GetCluster2());

    // the original copy of this will change throughout function... 
    ClusterList matchedClusters(element.GetOverlapResult().GetMatchedClusterList());
        
    //std::cout << "matched cluster list size: " << matchedClusters.size() << std::endl;

    ClusterSet checkedClusters;
    checkedClusters.insert(pSeedCluster);
        
    while (checkedClusters.size() != matchedClusters.size())
    {
        const Cluster *pClusterToDelete(nullptr);
            
        unsigned int highestHit(0);
        for (const Cluster *const pMatchedCluster : matchedClusters)
        {
            if (checkedClusters.count(pMatchedCluster))
                continue;

            if (pMatchedCluster->GetNCaloHits() > highestHit)
            {
                pClusterToDelete = pMatchedCluster;
                highestHit = pMatchedCluster->GetNCaloHits();
            }
        }

        checkedClusters.insert(pClusterToDelete);

        if (!pClusterToDelete->IsAvailable())
            continue;

        CaloHitList caloHitList1, caloHitList2, caloHitList3;
        pCluster1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
        pCluster2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
        pSeedCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);
        pClusterToDelete->GetOrderedCaloHitList().FillCaloHitList(caloHitList3);

        XOverlap xOverlapObject(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
        float chiSquaredSum(0.f);
        unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
                        
        StatusCode status(pAlgorithm->PerformMatching(caloHitList1, caloHitList2, caloHitList3, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xOverlapObject));
                        
        if (status == STATUS_CODE_NOT_FOUND)
            continue;

        if (status != STATUS_CODE_SUCCESS)
            throw StatusCodeException(status);

        const float reducedChiSquared(chiSquaredSum / nSamplingPoints);
                        
        if (reducedChiSquared < 1.f)
        {
	  //std::cout << "before removal... "<< std::endl;
            pAlgorithm->RemoveThirdViewCluster(pClusterToDelete);
	    //std::cout << "after removal... "<< std::endl;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm,
               pAlgorithm->GetThirdViewClusterListName()));
                        
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pSeedCluster, pClusterToDelete));
        }
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMergeTool::CollectDeltaRayHits(const MatrixType::Element &element, const CartesianPointVector &deltaRayProjectedPositions,
    const ParticleFlowObject *const pParentMuon, CaloHitList &collectedHits) const
{
    const HitType &thirdViewHitType(LArClusterHelper::GetClusterHitType(element.GetOverlapResult().GetBestMatchedCluster()));
    
    CartesianPointVector muonProjectedPositions;
    this->ProjectMuonPositions(thirdViewHitType, pParentMuon, muonProjectedPositions);    

    ClusterList muonClusterList;
    LArPfoHelper::GetClusters(pParentMuon, thirdViewHitType, muonClusterList);
            
    if (muonClusterList.size() != 1)
        return STATUS_CODE_NOT_FOUND;
    
    const Cluster *const pMuonCluster(muonClusterList.front());

    const float projectedHitsFraction(static_cast<float>(muonProjectedPositions.size()) / pMuonCluster->GetNCaloHits());    

    CartesianVector muonDirection(0.f, 0.f, 0.f), positionOnMuon(0.f, 0.f, 0.f);
    if (projectedHitsFraction < 0.8f)
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResult(pMuonCluster, 40, slidingFitPitch);

        CartesianVector deltaRayVertex(0.f,0.f,0.f), muonVertex(0.f,0.f,0.f);
        this->GetClosestPositions(deltaRayProjectedPositions, pMuonCluster, deltaRayVertex, muonVertex);

        positionOnMuon = this->GetClosestPosition(muonVertex, muonProjectedPositions, pMuonCluster);

	    if (positionOnMuon.GetMagnitude() < std::numeric_limits<float>::epsilon())
            return STATUS_CODE_NOT_FOUND;

        float rL(0.f), rT(0.f);
        slidingFitResult.GetLocalPosition(positionOnMuon, rL, rT);
        slidingFitResult.GetGlobalFitDirection(rL, muonDirection);
    }

    CaloHitList muonCaloHitList;
    pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);

    bool hitsAdded(true);
    while (hitsAdded)
    {
        hitsAdded = false;

        for (const CaloHit *const pCaloHit : muonCaloHitList)
        {
            if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            
            const float distanceToDeltaRayHits(std::min(this->GetClosestDistance(pCaloHit, deltaRayProjectedPositions), this->GetClosestDistance(pCaloHit, collectedHits)));
            const float distanceToMuonHits((projectedHitsFraction < 0.8f) ? muonDirection.GetCrossProduct(hitPosition - positionOnMuon).GetMagnitude() :
                this->GetClosestDistance(pCaloHit, muonProjectedPositions));

	        if ((std::fabs(distanceToMuonHits - distanceToDeltaRayHits) > std::numeric_limits<float>::epsilon()) && (distanceToDeltaRayHits < distanceToMuonHits)
                    && (distanceToMuonHits > 1.f) && (distanceToDeltaRayHits < 1.f))
	        {
                collectedHits.push_back(pCaloHit);
                hitsAdded = true;
            }
        }
    }

    // Catch if delta ray has travelled along muon
    if ((static_cast<float>(collectedHits.size()) / muonCaloHitList.size()) > 0.05)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TwoViewDeltaRayMergeTool::SplitCluster(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const Cluster *const pMuonCluster, CaloHitList &collectedHits) const
{
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, pAlgorithm->GetThirdViewClusterListName()));

    const Cluster *pDeltaRayCluster(nullptr);

    CaloHitList muonCaloHitList;
    pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);
    
    for (const CaloHit *const pCaloHit : muonCaloHitList)
    {
        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*pAlgorithm, pMuonCluster, pCaloHit));

            if (!pDeltaRayCluster)
            {
                const ClusterList *pTemporaryList(nullptr);
                std::string temporaryListName, currentListName;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*pAlgorithm, currentListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*pAlgorithm, pTemporaryList, temporaryListName));

                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, parameters, pDeltaRayCluster));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*pAlgorithm, temporaryListName, currentListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, currentListName));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*pAlgorithm, pDeltaRayCluster, pCaloHit));
            }
        }
    }

    return pDeltaRayCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMergeTool::ProjectMuonPositions(const HitType &thirdViewHitType, const ParticleFlowObject *const pParentMuon,
    CartesianPointVector &projectedPositions) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == thirdViewHitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == thirdViewHitType) || (hitType1 == hitType2))
                continue;

            ClusterList muonClusterList1, muonClusterList2;
            LArPfoHelper::GetClusters(pParentMuon, hitType1, muonClusterList1);
            LArPfoHelper::GetClusters(pParentMuon, hitType2, muonClusterList2);
            
            if ((muonClusterList1.size() != 1) || (muonClusterList1.size() != 1))
                return;

            this->ProjectPositions(muonClusterList1.front(), muonClusterList2.front(), projectedPositions);

            break;
        }
        
        break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMergeTool::ProjectPositions(const Cluster *const pCluster1, const Cluster *const pCluster2, CartesianPointVector &projectedPositions) const
{
    float m_xOverlapWindow(1.f);
    
    float xMin1(-std::numeric_limits<float>::max()), xMax1(+std::numeric_limits<float>::max());
    float xMin2(-std::numeric_limits<float>::max()), xMax2(+std::numeric_limits<float>::max());

    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);

    const float xPitch(0.5 * m_xOverlapWindow);

    const float xMin(std::max(xMin1, xMin2) - xPitch);
    const float xMax(std::min(xMax1, xMax2) + xPitch);
    const float xOverlap(xMax - xMin);

    // this has already been done...
    if (xOverlap < std::numeric_limits<float>::epsilon())
      return;
    
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

    if (hitType1 == hitType2)
        return;

    const unsigned int nPoints(1 + static_cast<unsigned int>(xOverlap / xPitch));

    // Projection into third view
    for (unsigned int n = 0; n < nPoints; ++n)
    {
        const float x(xMin + (xMax - xMin) * (static_cast<float>(n) + 0.5f) / static_cast<float>(nPoints));
        const float xmin(x - xPitch);
        const float xmax(x + xPitch);

        try
        {
            float zMin1(0.f), zMin2(0.f), zMax1(0.f), zMax2(0.f);
            pCluster1->GetClusterSpanZ(xmin, xmax, zMin1, zMax1);
            pCluster2->GetClusterSpanZ(xmin, xmax, zMin2, zMax2);

            const float z1(0.5f * (zMin1 + zMax1));
            const float z2(0.5f * (zMin2 + zMax2));

            // could make use of the chi-squared?
            float chi2(0.f);
            CartesianVector projection(0.f, 0.f, 0.f);
            LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, CartesianVector(x, 0.f, z1), CartesianVector(x, 0.f, z2), projection, chi2);

            projectedPositions.push_back(projection);
        }
        catch(StatusCodeException &statusCodeException)
        {
            if (statusCodeException.GetStatusCode() != STATUS_CODE_NOT_FOUND)
                throw statusCodeException.GetStatusCode();

            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMergeTool::ReadSettings(const TiXmlHandle xmlHandle)
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

void TwoViewDeltaRayMergeTool::GetBestMatchedAvailableCluster(const ClusterList &matchedClusters, const Cluster *&pBestMatchedCluster) const
{
    unsigned int highestNHits(0);
    for (const Cluster *const pMatchedCluster : matchedClusters)
    {
        if (!pMatchedCluster->IsAvailable())
            continue;
        
        if (pMatchedCluster->GetNCaloHits() > highestNHits)
        {
            highestNHits = pMatchedCluster->GetNCaloHits();
            pBestMatchedCluster = pMatchedCluster;
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoViewDeltaRayMergeTool::GetClosestPosition(const CartesianVector &referencePoint, const CartesianPointVector &cartesianPointVector, const Cluster *const pCluster) const
{
    CartesianVector closestPoint(0.f,0.f,0.f);
    float shortestDistanceSquared(std::numeric_limits<float>::max());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        if (LArClusterHelper::GetClosestDistance(testPosition, pCluster) > 0.5f)
            continue;

        const float separationSquared((testPosition - referencePoint).GetMagnitude());

        if (separationSquared > 5.f)
            continue;

        if (separationSquared < shortestDistanceSquared)
        {
            shortestDistanceSquared = separationSquared;
            closestPoint = testPosition;
        }
    }

    return closestPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoViewDeltaRayMergeTool::GetClosestDistance(const CaloHit *const pCaloHit, const CartesianPointVector &cartesianPointVector) const
{
    float shortestDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector referencePoint(pCaloHit->GetPositionVector());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        const float separationSquared((testPosition - referencePoint).GetMagnitudeSquared());

        if (separationSquared < shortestDistanceSquared)
            shortestDistanceSquared = separationSquared;
    }

    return std::sqrt(shortestDistanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------    

float TwoViewDeltaRayMergeTool::GetClosestDistance(const CaloHit *const pCaloHit, const CaloHitList &caloHitList) const
{
    float shortestDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector referencePoint(pCaloHit->GetPositionVector());

    for (const CaloHit *const pTestCaloHit : caloHitList)
    {
        const CartesianVector &position(pTestCaloHit->GetPositionVector());
        float separationSquared((position - referencePoint).GetMagnitudeSquared());

        if (separationSquared < shortestDistanceSquared)
            shortestDistanceSquared = separationSquared;
    }

    return std::sqrt(shortestDistanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMergeTool::GetClosestPositions(const CartesianPointVector &pCluster1, const Cluster *const pCluster2, CartesianVector &outputPosition1,
    CartesianVector &outputPosition2) const
{
    bool distanceFound(false);
    float minDistanceSquared(std::numeric_limits<float>::max());

    CartesianVector closestPosition1(0.f, 0.f, 0.f);
    CartesianVector closestPosition2(0.f, 0.f, 0.f);

    CaloHitList caloHitList2;
    pCluster2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);

    for (const CartesianVector &positionVector1 : pCluster1)
    {
        for (const CaloHit *const pCaloHit : caloHitList2)
        {
            const CartesianVector &positionVector2(pCaloHit->GetPositionVector());

            const float distanceSquared((positionVector1 - positionVector2).GetMagnitudeSquared());

            if (distanceSquared < minDistanceSquared)
            {
                minDistanceSquared = distanceSquared;
                closestPosition1 = positionVector1;
                closestPosition2 = positionVector2;
                distanceFound = true;
            }
        }
    }

    if (!distanceFound)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    outputPosition1 = closestPosition1;
    outputPosition2 = closestPosition2;
}


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
