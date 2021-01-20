/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the three view delta ray matching class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

TwoViewDeltaRayMatchingAlgorithm::TwoViewDeltaRayMatchingAlgorithm() :
    m_nMaxMatrixToolRepeats(10),
    m_minClusterCaloHits(3),
    m_maxDisplacementSquared(1.0f),
    m_minProjectedPositions(3),
    m_maxDistanceFromPrediction(2.f),
    m_maxGoodMatchReducedChiSquared(1.f)
{
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::GetConnectedElements(const Cluster *const pClusterA, const bool hasAssociatedMuon, MatrixType::ElementList &elementList, ClusterSet &checkedClusters)
{
    if (checkedClusters.count(pClusterA))
        return;

    if (!pClusterA->IsAvailable())
        return;

    auto &theMatrix(this->GetMatchingControl().GetOverlapMatrix());
    
    const HitType &hitType(LArClusterHelper::GetClusterHitType(pClusterA));
    const unsigned int hitTypeIndex(m_matchingControl.m_hitTypeToIndexMap.at(hitType));
    
    if (hitTypeIndex == 1)
        checkedClusters.insert(pClusterA);
    
    auto &navigationMap((hitTypeIndex == 1) ? theMatrix.GetClusterNavigationMap12() : theMatrix.GetClusterNavigationMap21());

    auto iter = navigationMap.find(pClusterA);

    if (iter == navigationMap.end())
        throw StatusCodeException(STATUS_CODE_FAILURE);
    
    for (const Cluster *const pClusterB : iter->second)
    {
        if (checkedClusters.count(pClusterB))
            continue;

        if (!pClusterB->IsAvailable())
            continue; 

        const Cluster *const pCluster1((hitTypeIndex == 1) ? pClusterA : pClusterB);
        const Cluster *const pCluster2((hitTypeIndex == 1) ? pClusterB : pClusterA);

        if (pCluster1 == pCluster2)
            throw StatusCodeException(STATUS_CODE_FAILURE);      
            
        // now find in the tensor
        try
        {
            auto &overlapResult(theMatrix.GetOverlapResult(pCluster1, pCluster2));
            
            const PfoList &commonMuonPfoList(overlapResult.GetCommonMuonPfoList());

            if (!hasAssociatedMuon && commonMuonPfoList.size())
                continue;

            if (hasAssociatedMuon && commonMuonPfoList.empty())
                continue;
                
            bool found = false;
            for (const MatrixType::Element &t : elementList)
            {
                if ((t.GetCluster1() == pCluster1) && (t.GetCluster2() == pCluster2))
                    found = true;
            }

            if (!found)
            {
                MatrixType::Element element(pCluster1, pCluster2, overlapResult);
                elementList.push_back(element);
            }
            
            this->GetConnectedElements(pClusterB, hasAssociatedMuon, elementList, checkedClusters);
        }
        catch (StatusCodeException &)
        {
            continue;
        }
    }

    std::sort(elementList.begin(), elementList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void TwoViewDeltaRayMatchingAlgorithm::GetUnambiguousElements(const bool hasAssociatedMuon, MatrixType::ElementList &elementList)
{
    auto &theMatrix(this->GetMatchingControl().GetOverlapMatrix());

    for (auto iter1 = theMatrix.begin(); iter1 != theMatrix.end(); ++iter1)
    {
        ClusterSet checkedClusters;
        MatrixType::ElementList tempElementList;
        this->GetConnectedElements(iter1->first, hasAssociatedMuon, tempElementList, checkedClusters);

        if (tempElementList.size() != 1)
            continue;

        MatrixType::Element &unambiguousElement(tempElementList.front());
        const Cluster *const pCluster1(unambiguousElement.GetCluster1()), *const pCluster2(unambiguousElement.GetCluster2());
        
        // ATTN With HIT_CUSTOM definitions, it is possible to navigate from different U clusters to same combination
        if (iter1->first != pCluster1)
            continue;

        if (!pCluster1 || !pCluster2)
            continue;

        auto iter2 = iter1->second.find(pCluster2);
        if (iter1->second.end() == iter2)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        MatrixType::Element element(pCluster1, pCluster2, iter2->second);
        elementList.push_back(element);
    }

    std::sort(elementList.begin(), elementList.end());
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMatchingAlgorithm::DoesClusterPassTesorThreshold(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
        return false;

    return true;
}    

//------------------------------------------------------------------------------------------------------------------------------------------
    
void TwoViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const)
{
    TrackTwoViewTopologyOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pCluster1, pCluster2, overlapResult));
    
    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapMatrix().SetOverlapResult(pCluster1, pCluster2, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------     

StatusCode TwoViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, TrackTwoViewTopologyOverlapResult &overlapResult) const
{
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1); 
    pCluster2->GetClusterSpanX(xMin2, xMax2); 

    const float overlapMinX(std::max(xMin1, xMin2));
    const float overlapMaxX(std::min(xMax1, xMax2));
    const float xOverlap(overlapMaxX - overlapMinX);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;

    PfoList commonMuonPfoList;
    this->FindCommonMuonParents(pCluster1, pCluster2, commonMuonPfoList);
    
    if (commonMuonPfoList.size())
        return STATUS_CODE_NOT_FOUND;

    // Project delta ray clusters into the theird view
    CartesianPointVector projectedPositions;
    StatusCode status(this->GetProjectedPositions(pCluster1, pCluster2, projectedPositions));
    
    if (status != STATUS_CODE_SUCCESS)
        return status;

    // Find all matched clusters (including unavailable)
    ClusterList matchedClusterList;
    this->CollectThirdViewClusters(pCluster1, pCluster2, projectedPositions, matchedClusterList);

    if (matchedClusterList.empty())
        return STATUS_CODE_NOT_FOUND;
    
    const Cluster *pBestMatchedCluster(nullptr);
    float reducedChiSquared(std::numeric_limits<float>::max());    
    this->GetBestMatchedCluster(pCluster1, pCluster2, commonMuonPfoList, matchedClusterList, pBestMatchedCluster, reducedChiSquared);

    //ATTN: Ignore if other clusters matches have more hits
    if (pBestMatchedCluster && (pBestMatchedCluster->IsAvailable()))
    {
        unsigned int hitSum12(pCluster1->GetNCaloHits() + pCluster2->GetNCaloHits());
        unsigned int hitSum13(pCluster1->GetNCaloHits() + pBestMatchedCluster->GetNCaloHits());
        unsigned int hitSum23(pCluster2->GetNCaloHits() + pBestMatchedCluster->GetNCaloHits());        

        if (hitSum23 > hitSum12)
            return STATUS_CODE_NOT_FOUND;

        if (hitSum13 > hitSum12)
            return STATUS_CODE_NOT_FOUND;
    }

    TwoViewXOverlap xOverlapObject(xMin1, xMax1, xMin2, xMax2);

    overlapResult = TrackTwoViewTopologyOverlapResult(xOverlapObject, commonMuonPfoList, pBestMatchedCluster, matchedClusterList, reducedChiSquared);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::FindCommonMuonParents(const Cluster *const pCluster1, const Cluster *const pCluster2, PfoList &commonMuonPfoList) const
{
    ClusterList consideredClusters1, consideredClusters2;
    PfoList nearbyMuonPfos1, nearbyMuonPfos2;
    this->GetNearbyMuonPfos(pCluster1, consideredClusters1, nearbyMuonPfos1);
    this->GetNearbyMuonPfos(pCluster2, consideredClusters2, nearbyMuonPfos2);

    for (const ParticleFlowObject *const pNearbyMuon1 : nearbyMuonPfos1)
    {
        for (const ParticleFlowObject *const pNearbyMuon2 : nearbyMuonPfos2)
        {
            if (pNearbyMuon1 == pNearbyMuon2)
            {
                commonMuonPfoList.push_back(pNearbyMuon1);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::ProjectMuonPositions(const HitType &thirdViewHitType, const ParticleFlowObject *const pParentMuon,
    CartesianPointVector &projectedPositions) const
{
    const HitTypeVector hitTypeVector({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});

    ClusterList muonClusterList1, muonClusterList2;    
    for (const HitType &hitType1 : hitTypeVector)
    {
        if (hitType1 == thirdViewHitType)
            continue;
        
        for (const HitType &hitType2 : hitTypeVector)
        {
            if ((hitType2 == thirdViewHitType) || (hitType1 == hitType2))
                continue;

            LArPfoHelper::GetClusters(pParentMuon, hitType1, muonClusterList1);
            LArPfoHelper::GetClusters(pParentMuon, hitType2, muonClusterList2);
            
            if ((muonClusterList1.size() != 1) || (muonClusterList1.size() != 1))
                return STATUS_CODE_NOT_FOUND;
            
            break;
        }
    }

    return (this->GetProjectedPositions(muonClusterList1.front(), muonClusterList2.front(), projectedPositions));
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::GetProjectedPositions(const Cluster *const pCluster1, const Cluster *const pCluster2,
    CartesianPointVector &projectedPositions) const
{
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1); 
    pCluster2->GetClusterSpanX(xMin2, xMax2); 

    const float xPitch(0.5 * m_xOverlapWindow);
    const float xMin(std::max(xMin1, xMin2) - xPitch);
    const float xMax(std::min(xMax1, xMax2) + xPitch);
    const float xOverlap(xMax - xMin);

    if (xOverlap < std::numeric_limits<float>::epsilon())
         return STATUS_CODE_NOT_FOUND;
    
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

    if (hitType1 == hitType2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

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

            float chi2;
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

    // Reject if projection is not good
    if (projectedPositions.size() < m_minProjectedPositions)
        return STATUS_CODE_NOT_FOUND;
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void TwoViewDeltaRayMatchingAlgorithm::CollectThirdViewClusters(const Cluster *const pCluster1, const Cluster *const pCluster2, const CartesianPointVector &projectedPositions,
    ClusterList &matchedClusters) const
{
    const ClusterList *pInputClusterList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pInputClusterList));
    
    for (const Cluster *const pCluster : *pInputClusterList)
    {        
        const float separation(this->GetClosestDistance(pCluster, projectedPositions));
            
        if (separation > m_maxDistanceFromPrediction)
            continue;

        float reducedChiSquared(0.f);
        if (this->PerformThreeViewMatching(pCluster1, pCluster2, pCluster, reducedChiSquared) == STATUS_CODE_NOT_FOUND)
            continue;

        if (reducedChiSquared > m_maxGoodMatchReducedChiSquared)
            continue;
            
        matchedClusters.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::GetBestMatchedCluster(const Cluster *const pCluster1, const Cluster *const pCluster2, const PfoList &commonMuonPfoList,
    const ClusterList &matchedClusters, const Cluster *&pBestMatchedCluster, float &reducedChiSquared) const
{
    const ClusterList *pInputClusterList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pInputClusterList));
    
    ClusterList muonClusterList;
    
    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
        LArPfoHelper::GetClusters(pMuonPfo, LArClusterHelper::GetClusterHitType(pInputClusterList->front()), muonClusterList);

    ClusterSet checkedClusters;
    while (true)
    {    
        pBestMatchedCluster = nullptr; reducedChiSquared = std::numeric_limits<float>::max();

        unsigned int highestNHits(0);
        for (const Cluster *const pMatchedCluster : matchedClusters)
        {
            if (checkedClusters.count(pMatchedCluster))
                continue;

            if (!pMatchedCluster->IsAvailable())
            {
                if (std::find(muonClusterList.begin(), muonClusterList.end(), pMatchedCluster) == muonClusterList.end())
                    continue;
            }
                
            if (pMatchedCluster->GetNCaloHits() > highestNHits)
            {
                highestNHits = pMatchedCluster->GetNCaloHits();
                pBestMatchedCluster = pMatchedCluster;
            }
        }

        if (!pBestMatchedCluster)
            return;

        checkedClusters.insert(pBestMatchedCluster);

        if (this->PerformThreeViewMatching(pCluster1, pCluster2, pBestMatchedCluster, reducedChiSquared) == STATUS_CODE_NOT_FOUND)
        {
            if (std::find(muonClusterList.begin(), muonClusterList.end(), pBestMatchedCluster) != muonClusterList.end())
            {
                continue;
            }
            else
            {
                std::cout << "THIS SHOULD NEVER HAPPEN" << std::endl;
                throw STATUS_CODE_NOT_ALLOWED;
            }
        }

        if (reducedChiSquared > m_maxGoodMatchReducedChiSquared)
            continue;

        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (MatrixToolVector::const_iterator toolIter = m_algorithmToolVector.begin(); toolIter != m_algorithmToolVector.end(); )
    {
        if ((*toolIter)->Run(this, this->GetMatchingControl().GetOverlapMatrix()))
        {
            toolIter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxMatrixToolRepeats)
                break;
        }
        else
        {
            ++toolIter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName)); 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DeltaRayPfoListName", m_deltaRayPfoListName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "DeltaRayTools", algorithmToolVector));
    
    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        DeltaRayMatrixTool *const pDeltaRayMatrixTool(dynamic_cast<DeltaRayMatrixTool*>(*iter));

        if (!pDeltaRayMatrixTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pDeltaRayMatrixTool);
    }
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxMatrixToolRepeats", m_nMaxMatrixToolRepeats));    
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SearchRegion1D", m_searchRegion1D));      

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverlapWindow", m_xOverlapWindow));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDisplacementSquared", m_maxDisplacementSquared));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedPoints", m_minMatchedPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinProjectedPositions", m_minProjectedPositions));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromPrediction", m_maxDistanceFromPrediction)); 

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGoodMatchReducedChiSquared", m_maxGoodMatchReducedChiSquared));

                                                                        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));    

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------     
    
void TwoViewDeltaRayMatchingAlgorithm::RemoveThirdViewCluster(const Cluster *const pCluster)
{
    auto &theMatrix(this->GetMatchingControl().GetOverlapMatrix());

    for (auto iter1 = theMatrix.begin(); iter1 != theMatrix.end(); ++iter1)
    {
        auto iter2(iter1->second);
        
        for (auto &entry : iter2)
        {
            TrackTwoViewTopologyOverlapResult &overlapResult(entry.second);

            ClusterList matchedClusters(overlapResult.GetMatchedClusterList());

            auto matchedClustersIter(std::find(matchedClusters.begin(), matchedClusters.end(), pCluster));

            if (matchedClustersIter == matchedClusters.end())
                continue;
            
            matchedClusters.erase(matchedClustersIter);

            const Cluster *pBestMatchedCluster(nullptr); float reducedChiSquared(std::numeric_limits<float>::max());
            this->GetBestMatchedCluster(iter1->first, entry.first, overlapResult.GetCommonMuonPfoList(), matchedClusters, pBestMatchedCluster, reducedChiSquared);

            overlapResult = TrackTwoViewTopologyOverlapResult(overlapResult.GetXOverlap(), overlapResult.GetCommonMuonPfoList(), pBestMatchedCluster, matchedClusters, reducedChiSquared);
            theMatrix.ReplaceOverlapResult(iter1->first, entry.first, overlapResult);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMatchingAlgorithm::CreatePfo(const MatrixType::Element &element)
{
    ProtoParticle protoParticle;
    protoParticle.m_clusterList.push_back(element.GetCluster1());
    protoParticle.m_clusterList.push_back(element.GetCluster2());

    const Cluster *const pBestMatchedCluster(element.GetOverlapResult().GetBestMatchedCluster());

    if (pBestMatchedCluster)
    {
        this->GrowThirdView(element, protoParticle);
    }
    else
    {
        std::cout << "TwoViewCreation" << std::endl;
    }

    ProtoParticleVector protoParticleVector({protoParticle});

    return (this->CreateThreeDParticles(protoParticleVector));
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::GrowThirdView(const MatrixType::Element &element, ProtoParticle &protoParticle)
{
    const PfoList &commonMuonPfoList(element.GetOverlapResult().GetCommonMuonPfoList());
    const Cluster *const pBestMatchedCluster(element.GetOverlapResult().GetBestMatchedCluster());
    const HitType &thirdViewHitType(LArClusterHelper::GetClusterHitType(pBestMatchedCluster));

    // Determine whether best matched cluster is a muon
    const ParticleFlowObject *pMatchedMuonPfo(nullptr);
    for (const ParticleFlowObject *const pMuonPfo : commonMuonPfoList)
    {
        ClusterList muonClusterList;
        LArPfoHelper::GetClusters(pMuonPfo, thirdViewHitType, muonClusterList);

        if (std::find(muonClusterList.begin(), muonClusterList.end(), pBestMatchedCluster) != muonClusterList.end())
            pMatchedMuonPfo = pMuonPfo;
    }
    
    if (pMatchedMuonPfo)
    {
        CaloHitList deltaRayHits;
        if ((this->PullOutDeltaRayHits(element, pMatchedMuonPfo, deltaRayHits) != STATUS_CODE_SUCCESS) || (deltaRayHits.empty()))
        {
            const Cluster *pSeedCluster(nullptr);
            this->GetBestMatchedAvailableCluster(element.GetOverlapResult().GetMatchedClusterList(), pSeedCluster);
            
            if (pSeedCluster)
            {
                this->MergeThirdView(element, pSeedCluster);
                this->RemoveThirdViewCluster(pSeedCluster);
                
                protoParticle.m_clusterList.push_back(pSeedCluster);
            }
        }
        else
        {
            const Cluster *const pSeedCluster(this->SplitCluster(pBestMatchedCluster, deltaRayHits));

            this->MergeThirdView(element, pSeedCluster);
            this->RemoveThirdViewCluster(pSeedCluster);

            protoParticle.m_clusterList.push_back(pSeedCluster);
        }   
    }
    else
    {
        this->MergeThirdView(element, pBestMatchedCluster);
        this->RemoveThirdViewCluster(pBestMatchedCluster);

        protoParticle.m_clusterList.push_back(pBestMatchedCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMatchingAlgorithm::PullOutDeltaRayHits(const MatrixType::Element &element, const ParticleFlowObject *const pParentMuon, CaloHitList &collectedHits) const
{
    CartesianPointVector deltaRayProjectedPositions;
    if (this->GetProjectedPositions(element.GetCluster1(), element.GetCluster2(), deltaRayProjectedPositions) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;
    
    const HitType &thirdViewHitType(LArClusterHelper::GetClusterHitType(element.GetOverlapResult().GetBestMatchedCluster()));
    
    CartesianPointVector muonProjectedPositions;
    if (this->ProjectMuonPositions(thirdViewHitType, pParentMuon, muonProjectedPositions) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

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

void TwoViewDeltaRayMatchingAlgorithm::MergeThirdView(const MatrixType::Element &element, const Cluster *const pSeedCluster)
{
    const Cluster *const pCluster1(element.GetCluster1()), *const pCluster2(element.GetCluster2());

    // the original copy of this will change throughout function... 
    ClusterList matchedClusters(element.GetOverlapResult().GetMatchedClusterList());
        
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

        float reducedChiSquared(0.f);
        if (this->PerformThreeViewMatching(pCluster1, pCluster2, pClusterToDelete, reducedChiSquared) == STATUS_CODE_NOT_FOUND)
            continue;

        if (reducedChiSquared > m_maxGoodMatchReducedChiSquared)
            continue;

        this->RemoveThirdViewCluster(pClusterToDelete);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this,
            this->GetThirdViewClusterListName()));
                        
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pClusterToDelete));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TwoViewDeltaRayMatchingAlgorithm::SplitCluster(const Cluster *const pMuonCluster, CaloHitList &collectedHits) const
{
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, this->GetThirdViewClusterListName()));

    const Cluster *pDeltaRayCluster(nullptr);

    CaloHitList muonCaloHitList;
    pMuonCluster->GetOrderedCaloHitList().FillCaloHitList(muonCaloHitList);
    
    for (const CaloHit *const pCaloHit : muonCaloHitList)
    {
        if (std::find(collectedHits.begin(), collectedHits.end(), pCaloHit) != collectedHits.end())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pMuonCluster, pCaloHit));

            if (!pDeltaRayCluster)
            {
                const ClusterList *pTemporaryList(nullptr);
                std::string temporaryListName, currentListName;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*this, pTemporaryList, temporaryListName));

                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pDeltaRayCluster));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, temporaryListName, currentListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentListName));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pDeltaRayCluster, pCaloHit));
            }
        }
    }

    return pDeltaRayCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewDeltaRayMatchingAlgorithm::GetBestMatchedAvailableCluster(const ClusterList &matchedClusters, const Cluster *&pBestMatchedCluster) const
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

float TwoViewDeltaRayMatchingAlgorithm::GetClosestDistance(const Cluster *const pCluster, const CartesianPointVector &cartesianPointVector) const
{
    float closestDistance(std::numeric_limits<float>::max());
    
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float distance(this->GetClosestDistance(pCaloHit, cartesianPointVector));

        if (distance < closestDistance)
            closestDistance = distance;
    }

    return closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoViewDeltaRayMatchingAlgorithm::GetClosestDistance(const CaloHit *const pCaloHit, const CartesianPointVector &cartesianPointVector) const
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

float TwoViewDeltaRayMatchingAlgorithm::GetClosestDistance(const CaloHit *const pCaloHit, const CaloHitList &caloHitList) const
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

void TwoViewDeltaRayMatchingAlgorithm::GetClosestPositions(const CartesianPointVector &pCluster1, const Cluster *const pCluster2, CartesianVector &outputPosition1,
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

CartesianVector TwoViewDeltaRayMatchingAlgorithm::GetClosestPosition(const CartesianVector &referencePoint, const CartesianPointVector &cartesianPointVector, const Cluster *const pCluster) const
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

} // namespace lar_content
