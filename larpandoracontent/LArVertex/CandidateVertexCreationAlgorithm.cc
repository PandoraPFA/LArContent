/**
 *  @file   larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the candidate vertex creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.h"

#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TProfile.h"

#include <utility> 

using namespace pandora;

namespace lar_content
{

CandidateVertexCreationAlgorithm::CandidateVertexCreationAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_slidingFitWindow(20),
    m_minClusterCaloHits(0),
    m_minClusterLengthSquared(0.f * 0.f),
    m_maxClusterXDiscrepancy(4.f),
    m_chiSquaredCut(2.f),
    m_enableCrossingCandidates(false),
    m_enableEnergyCandidates(false),
    m_energyPlot(false),
    m_minCrossingClusterSize(10),
    m_extrapolationLength(20.0f),
    m_extrapolationStepSize(0.1f),
    m_minClusterCrossingApproach(0.25f),
    m_postCrossingSkipDistance(5.0f),
    m_minEnergyVertexClusterSize(60)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::Run()
{
    try
    {
        ClusterVector clusterVectorU, clusterVectorV, clusterVectorW;
        this->SelectClusters(clusterVectorU, clusterVectorV, clusterVectorW);

        const VertexList *pVertexList(NULL); std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

        this->CreateClusterEndPointComparisonVertices(clusterVectorU, clusterVectorV);
        this->CreateClusterEndPointComparisonVertices(clusterVectorU, clusterVectorW);
        this->CreateClusterEndPointComparisonVertices(clusterVectorV, clusterVectorW);
        
        if (m_enableCrossingCandidates)
            this->CreateCrossingVertices(clusterVectorU, clusterVectorV, clusterVectorW);
            
        if (m_enableEnergyCandidates)
        {
            this->CreateEnergyVertices(clusterVectorU, clusterVectorV, clusterVectorW);
            this->CreateEnergyVertices(clusterVectorV, clusterVectorW, clusterVectorU);
            this->CreateEnergyVertices(clusterVectorW, clusterVectorU, clusterVectorV);
        }

        if (!pVertexList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_vertexListName));

            if (m_replaceCurrentVertexList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_vertexListName));
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }

    this->TidyUp();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::SelectClusters(ClusterVector &clusterVectorU, ClusterVector &clusterVectorV, ClusterVector &clusterVectorW)
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CandidateVertexCreationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        ClusterVector &selectedClusterVector((TPC_VIEW_U == hitType) ? clusterVectorU : (TPC_VIEW_V == hitType) ? clusterVectorV : clusterVectorW);

        if (!selectedClusterVector.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterVector sortedClusters(pClusterList->begin(), pClusterList->end());
        std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pCluster : sortedClusters)
        {
            if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
                continue;

            if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
                continue;

            try
            {
                this->AddToSlidingFitCache(pCluster);
                selectedClusterVector.push_back(pCluster);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateClusterEndPointComparisonVertices(const ClusterVector &clusterVector1, const ClusterVector &clusterVector2) const
{
    for (const Cluster *const pCluster1 : clusterVector1)
    {
        const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));

        const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(pCluster1));
        const CartesianVector minLayerPosition1(fitResult1.GetGlobalMinLayerPosition());
        const CartesianVector maxLayerPosition1(fitResult1.GetGlobalMaxLayerPosition());

        for (const Cluster *const pCluster2 : clusterVector2)
        {
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

            const TwoDSlidingFitResult &fitResult2(this->GetCachedSlidingFitResult(pCluster2));
            const CartesianVector minLayerPosition2(fitResult2.GetGlobalMinLayerPosition());
            const CartesianVector maxLayerPosition2(fitResult2.GetGlobalMaxLayerPosition());

            this->CreateVertex(maxLayerPosition1, hitType1, fitResult2);
            this->CreateVertex(minLayerPosition1, hitType1, fitResult2);
            this->CreateVertex(maxLayerPosition2, hitType2, fitResult1);
            this->CreateVertex(minLayerPosition2, hitType2, fitResult1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateCrossingVertices(ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW)
{
    std::vector<CartesianVector> crossingsU;
    std::vector<CartesianVector> crossingsV;
    std::vector<CartesianVector> crossingsW;
    
    this->Find2DClusterCrossings(clusterListU, crossingsU);
    this->Find2DClusterCrossings(clusterListV, crossingsV);
    this->Find2DClusterCrossings(clusterListW, crossingsW);
    
    this->CreateMatchedVertices(crossingsU, crossingsV, TPC_VIEW_U, TPC_VIEW_V);
    this->CreateMatchedVertices(crossingsU, crossingsW, TPC_VIEW_U, TPC_VIEW_W);
    this->CreateMatchedVertices(crossingsV, crossingsW, TPC_VIEW_V, TPC_VIEW_W);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateEnergyVertices(const ClusterList &clusterList1, const ClusterList &clusterList2, const ClusterList &clusterList3)
{
    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iter1End = clusterList1.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pCluster = *iter1;
        
        if (pCluster->GetNCaloHits() < m_minEnergyVertexClusterSize)
            continue;
        
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.GetCaloHitList(caloHitList);
        
        std::vector<CartesianVector> energyAlongRLvector;
        this->CreateEnergyAlongRLVector(pCluster, energyAlongRLvector);
        
        std::vector<CartesianVector> filteredEnergyAlongRLvector;
        this->FilterEnergyVector(energyAlongRLvector, filteredEnergyAlongRLvector);
        
        std::vector<float> energySpikeRLvector;
        this->FindEnergySpike(filteredEnergyAlongRLvector, energySpikeRLvector);
        
        this->CreateVerticesFromSpikes(energySpikeRLvector, filteredEnergyAlongRLvector, caloHitList, clusterList1, clusterList2, clusterList3);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateVertex(const CartesianVector &position1, const HitType hitType1, const TwoDSlidingFitResult &fitResult2) const
{
    const CartesianVector minLayerPosition2(fitResult2.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition2(fitResult2.GetGlobalMaxLayerPosition());
    
    if ((((position1.GetX() < minLayerPosition2.GetX()) && (position1.GetX() < maxLayerPosition2.GetX())) ||
        ((position1.GetX() > minLayerPosition2.GetX()) && (position1.GetX() > maxLayerPosition2.GetX()))) &&
        (std::fabs(position1.GetX() - minLayerPosition2.GetX()) > m_maxClusterXDiscrepancy) &&
        (std::fabs(position1.GetX() - maxLayerPosition2.GetX()) > m_maxClusterXDiscrepancy))
    {
        return;
    }

    CartesianVector position2(0.f, 0.f, 0.f);
    if (STATUS_CODE_SUCCESS != fitResult2.GetExtrapolatedPositionAtX(position1.GetX(), position2))
        return;

    const HitType hitType2(LArClusterHelper::GetClusterHitType(fitResult2.GetCluster()));

    float chiSquared(0.f);
    CartesianVector position3D(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, position1, position2, position3D, chiSquared);
    const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_W));

    if (chiSquared > m_chiSquaredCut)
        return;
        
    //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionW, "Old Vertex", RED, 1));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = position3D;
    parameters.m_vertexLabel = VERTEX_INTERACTION;
    parameters.m_vertexType = VERTEX_3D;

    const Vertex *pVertex(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::Find2DClusterCrossings(const ClusterList &clusterList, std::vector<CartesianVector> &crossingsVector)
{
    for (ClusterList::const_iterator iter1 = clusterList.begin(), iter1End = clusterList.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pCluster1 = *iter1;
        
        if (pCluster1->GetNCaloHits() < m_minCrossingClusterSize)
            continue;
        
        const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(*iter1));
        
        const CartesianVector minLayerPosition1(fitResult1.GetGlobalMinLayerPosition());
        const CartesianVector maxLayerPosition1(fitResult1.GetGlobalMaxLayerPosition());
        
        CartesianVector minZHitPosition(0.f, 0.f, 0.f), maxZHitPosition(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(clusterList, minZHitPosition, maxZHitPosition);
        
        crossingsVector.push_back(minLayerPosition1);
        crossingsVector.push_back(maxLayerPosition1);
        
        crossingsVector.push_back(minZHitPosition);
        crossingsVector.push_back(maxZHitPosition);
        
        std::vector<CartesianVector> spacePointVector1;
        this->GetExtrapolatedClusterSpacepoints(spacePointVector1, pCluster1);
        
        for (ClusterList::const_iterator iter2 = clusterList.begin(), iter2End = clusterList.end(); iter2 != iter2End; ++iter2)
        {
            const Cluster *const pCluster2 = *iter2;
            
            if (pCluster1->GetNCaloHits() == pCluster2->GetNCaloHits() || pCluster2->GetNCaloHits() < m_minCrossingClusterSize)
                continue;

            std::vector<CartesianVector> spacePointVector2;
            this->GetExtrapolatedClusterSpacepoints(spacePointVector2, pCluster2);
            
            std::vector<float> distancesBetweenClusters;
            
            for (CartesianVector &position1: spacePointVector1)
            {
                for (CartesianVector &position2: spacePointVector2)
                    distancesBetweenClusters.push_back((position1-position2).GetMagnitude());
            }
            
            std::sort(spacePointVector1.begin(), spacePointVector1.end(), SortSpacePointsByZ);
            std::sort(spacePointVector2.begin(), spacePointVector2.end(), SortSpacePointsByZ);
            
            this->FindCrossingsFromSpacepoints(spacePointVector1, spacePointVector2, crossingsVector);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::GetExtrapolatedClusterSpacepoints(std::vector<CartesianVector> &spacePointVector, const Cluster *const pCluster)
{
    OrderedCaloHitList orderedCaloHitList1(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList1;
    orderedCaloHitList1.GetCaloHitList(caloHitList1);
    
    for (const CaloHit *const pCaloHit : caloHitList1)
    {
        CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        spacePointVector.push_back(caloHitPosition);
    }
    
    const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(pCluster));
    
    float minLayerRL(fitResult1.GetL(fitResult1.GetMinLayer()));
    float maxLayerRL(fitResult1.GetL(fitResult1.GetMaxLayer()));
    
    for (float i = m_extrapolationStepSize; i < m_extrapolationLength; i += m_extrapolationStepSize)
    {
        CartesianVector tempExtrapolatedPositionUnder(0.f, 0.f, 0.f);
        CartesianVector tempExtrapolatedPositionOver(0.f, 0.f, 0.f);
    
        fitResult1.GetExtrapolatedPosition(minLayerRL - i, tempExtrapolatedPositionUnder);
        fitResult1.GetExtrapolatedPosition(maxLayerRL + i, tempExtrapolatedPositionOver);
    
        spacePointVector.push_back(tempExtrapolatedPositionUnder);
        spacePointVector.push_back(tempExtrapolatedPositionOver);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindCrossingsFromSpacepoints(std::vector<CartesianVector> &spacePointVector1, std::vector<CartesianVector> &spacePointVector2, std::vector<CartesianVector> &crossingsVector)
{
    int skipCounter(0);
    bool shouldSkip(false);
    
    for (CartesianVector &position1: spacePointVector1)
    {
        if (shouldSkip)
        {
            if (skipCounter <= (m_postCrossingSkipDistance/m_extrapolationStepSize))
            {
                skipCounter++;
                continue;
            }
        }
        
        skipCounter = 0;
        shouldSkip = false;
        
        for (CartesianVector &position2: spacePointVector2)
        {
            if ((position1-position2).GetMagnitude() < m_minClusterCrossingApproach)
            {
                crossingsVector.push_back(position1);
                crossingsVector.push_back(position2);
                
                shouldSkip = true;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateMatchedVertices(std::vector<CartesianVector> &crossingsVector1, std::vector<CartesianVector> &crossingsVector2, HitType hitType1, HitType hitType2) const
{
    for (CartesianVector &position1: crossingsVector1)
    {
        std::vector<std::pair<CartesianVector*,float>> matched3DPositions;
        std::vector<float> chiSquareds;
        
        for (CartesianVector &position2: crossingsVector2)
        {
            if (std::fabs(position1.GetX() - position2.GetX()) > 0.1)
                continue;
                
            float chiSquared(0.f);
            CartesianVector position3D(0.f, 0.f, 0.f);
            LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, position1, position2, position3D, chiSquared);
            const CartesianVector vertexProjectionU(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_U));
            const CartesianVector vertexProjectionV(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_V));
            const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_W));
                
            if (chiSquared > 2.0)
                return;
            
            std::pair<CartesianVector*,float> positionChiSquaredPair;
            positionChiSquaredPair = std::make_pair(&position3D, chiSquared);
            
            matched3DPositions.push_back(positionChiSquaredPair);
            chiSquareds.push_back(chiSquared);
        }
        
        for (std::pair<CartesianVector*,float> &pair : matched3DPositions)
        {
            if (pair.second == (*std::min_element(chiSquareds.begin(), chiSquareds.end())))
            {
                PandoraContentApi::Vertex::Parameters parameters;
                parameters.m_position = *(pair.first);
                parameters.m_vertexLabel = VERTEX_INTERACTION;
                parameters.m_vertexType = VERTEX_3D;
                    
                const Vertex *pVertex(NULL);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateEnergyAlongRLVector(const Cluster *const pCluster, std::vector<CartesianVector> &energyAlongRLvector)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.GetCaloHitList(caloHitList);
    
    const TwoDSlidingFitResult &slidingFitResult(this->GetCachedSlidingFitResult(pCluster));
    
    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CaloHit *const pCaloHit(*hitIter);
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        
        float caloHitEnergy(pCaloHit->GetElectromagneticEnergy());
        float rL(0.f), rT(0.f);
        
        slidingFitResult.GetLocalPosition(caloHitPosition, rL, rT);
        CartesianVector energyAlongRL(rL, 0.f, 1000*caloHitEnergy); //units of MeV
        
        energyAlongRLvector.push_back(energyAlongRL);
    }
    
    std::sort(energyAlongRLvector.begin(), energyAlongRLvector.end(), SortEnergyVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FilterEnergyVector(const std::vector<CartesianVector> &unfilteredEnergyVector, std::vector<CartesianVector> &filteredEnergyVector)
{
    for (std::vector<CartesianVector>::const_iterator pairIter = std::next(unfilteredEnergyVector.begin(), 1), pairIterEnd = std::prev(unfilteredEnergyVector.end(), 1); pairIter != pairIterEnd; ++pairIter)
    {
        float chargeRatioNext(((*(std::next(pairIter, 1))).GetZ())/((*pairIter).GetZ()));
        float chargeRatioPrevious(((*(std::prev(pairIter, 1))).GetZ())/((*pairIter).GetZ()));
        
        if (!(chargeRatioPrevious < 0.8 && chargeRatioNext < 0.8))
            filteredEnergyVector.push_back(*pairIter);
    }
    
    std::sort(filteredEnergyVector.begin(), filteredEnergyVector.end(), SortEnergyVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::BinEnergyRLVector(const std::vector<CartesianVector> &energyAlongRLvector, std::vector<CartesianVector> &binnedEnergyAlongRLvector)
{
    int clusterLength(100);
    for (float i = 0.5; i != clusterLength; i += 0.5)
    {
        int nHitsBin(0);
        float averageBinPosition(0.f);
        float averageBinEnergy(0.f);
        
        for (const CartesianVector &energyRL : energyAlongRLvector)
        {
            if (energyRL.GetX() > i)
                break;
            
            if (energyRL.GetX() < i && energyRL.GetX() > (i - 1))
            {
                nHitsBin++;
                averageBinPosition += energyRL.GetX();
                averageBinEnergy += energyRL.GetZ();
            }
        }
        
        if (nHitsBin == 0)
            continue;
        
        averageBinPosition /= nHitsBin;
        averageBinEnergy /= nHitsBin;
        
        CartesianVector bin(averageBinPosition, 0.f, averageBinEnergy);
        binnedEnergyAlongRLvector.push_back(bin);

    }
    
    std::sort(binnedEnergyAlongRLvector.begin(), binnedEnergyAlongRLvector.end(), SortEnergyVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindEnergySpike(const std::vector<CartesianVector> &energyAlongRLvector, std::vector<float> &spikeRLvector)
{
    std::vector<CartesianVector> binnedEnergyAlongRLvector;
    this->BinEnergyRLVector(energyAlongRLvector, binnedEnergyAlongRLvector);
    
    std::vector<float> binRLvector;
    this->FindBinWithSpike(binnedEnergyAlongRLvector, binRLvector);
    
    this->ConvertBinRLToSpikeRL(binRLvector, spikeRLvector, energyAlongRLvector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindBinWithSpike(const std::vector<CartesianVector> &binnedEnergyAlongRLvector, std::vector<float> &binRLvector)
{
    float averageEnergyDifference(0.f);
    
    for (std::vector<CartesianVector>::const_iterator pairIter = std::next(binnedEnergyAlongRLvector.begin(), 1), pairIterEnd = std::prev(binnedEnergyAlongRLvector.end(), 1); pairIter != pairIterEnd; ++pairIter)
    {
        float thisBinAverageEnergy((*pairIter).GetZ());
        float nextBinAverageEnergy((*std::next(pairIter, 1)).GetZ());
        float energyDifference(std::abs(nextBinAverageEnergy / thisBinAverageEnergy));
        averageEnergyDifference += (energyDifference);
    }
    
    averageEnergyDifference /= binnedEnergyAlongRLvector.size();
    
    if (averageEnergyDifference > 1.5)
        return;
        
    static float bestScore(1.5f);
    
    for (std::vector<CartesianVector>::const_iterator pairIter = std::next(binnedEnergyAlongRLvector.begin(), 1), pairIterEnd = std::prev(binnedEnergyAlongRLvector.end(), 2); pairIter != pairIterEnd; ++pairIter)
    {
        float thisBinAveragePosition((*pairIter).GetX());
        float thisBinAverageEnergy((*pairIter).GetZ());
        
        float previousBinAverageEnergy((*std::prev(pairIter, 1)).GetZ());
        float nextBinAverageEnergy((*std::next(pairIter, 1)).GetZ());
        
        float previousPreviousBinAverageEnergy((*std::prev(pairIter, 2)).GetZ());
        float nextNextBinAverageEnergy((*std::next(pairIter, 2)).GetZ());
        
        float previousPreviousPreviousBinAverageEnergy((*std::prev(pairIter, 3)).GetZ());
        float nextNextNextBinAverageEnergy((*std::next(pairIter, 3)).GetZ());
        
        float previousPreviousPreviousPreviousBinAverageEnergy((*std::prev(pairIter, 4)).GetZ());
        float nextNextNextNextBinAverageEnergy((*std::next(pairIter, 4)).GetZ());
        

        //std::cout << "Jump position: " << thisBinAveragePosition << std::endl;
        //std::cout << "Jump ratio: " << std::abs(1 - std::abs(nextBinAverageEnergy / thisBinAverageEnergy)) << std::endl;
        //std::cout << "Next jump ratio: " << std::abs(1 - (std::abs(nextNextBinAverageEnergy / thisBinAverageEnergy))) << std::endl;
        //std::cout << "Previous jump ratio: " << std::abs(1 - std::abs(previousBinAverageEnergy / thisBinAverageEnergy)) << std::endl;
        //std::cout << "Previous previous jump ratio: " << std::abs(1 - std::abs(previousPreviousBinAverageEnergy / thisBinAverageEnergy)) << std::endl;
        
        //This if statement means: can we find a bin from which one one side the bin energies steadily rise and on the other side there are two bins close together in energy in either the + or - RL direction?
        if ((std::abs(1 - std::abs(nextBinAverageEnergy / thisBinAverageEnergy)) > 0.15 && std::abs(1 - (std::abs(nextNextBinAverageEnergy / thisBinAverageEnergy))) > 0.5 && std::abs(1 - (std::abs(nextNextNextBinAverageEnergy / thisBinAverageEnergy))) > 0.7
        && std::abs(1 - (std::abs(previousBinAverageEnergy / thisBinAverageEnergy))) < 0.15 && std::abs(1 - (std::abs(previousPreviousBinAverageEnergy / thisBinAverageEnergy))) < 0.3)
        || (std::abs(1 - std::abs(previousBinAverageEnergy / thisBinAverageEnergy)) > 0.15 && std::abs(1 - (std::abs(previousPreviousBinAverageEnergy / thisBinAverageEnergy))) > 0.5 && std::abs(1 - (std::abs(previousPreviousPreviousBinAverageEnergy / thisBinAverageEnergy))) > 0.7
        && std::abs(1 - (std::abs(nextBinAverageEnergy / thisBinAverageEnergy))) < 0.15 && std::abs(1 - (std::abs(nextNextBinAverageEnergy / thisBinAverageEnergy))) < 0.3))
        {
            if (thisBinAverageEnergy < 0.25 || thisBinAverageEnergy == 4.375)
                continue;
            
            float score((std::abs(std::abs(nextBinAverageEnergy / thisBinAverageEnergy)) + std::abs((std::abs(nextNextBinAverageEnergy / thisBinAverageEnergy)))
            + std::abs((std::abs(nextNextNextBinAverageEnergy / thisBinAverageEnergy))) + std::abs((std::abs(nextNextNextNextBinAverageEnergy / thisBinAverageEnergy)))) 
            / (std::abs((std::abs(previousBinAverageEnergy / thisBinAverageEnergy))) + std::abs((std::abs(previousPreviousBinAverageEnergy / thisBinAverageEnergy)))));
            
            float scoreTwo((std::abs(std::abs(previousBinAverageEnergy / thisBinAverageEnergy)) + std::abs((std::abs(previousPreviousBinAverageEnergy / thisBinAverageEnergy)))
            + std::abs((std::abs(previousPreviousPreviousBinAverageEnergy / thisBinAverageEnergy))) + std::abs((std::abs(previousPreviousPreviousPreviousBinAverageEnergy / thisBinAverageEnergy))))
            / (std::abs((std::abs(nextBinAverageEnergy / thisBinAverageEnergy))) + std::abs((std::abs(nextNextBinAverageEnergy / thisBinAverageEnergy)))));
            
            float workingScore(0.f);
            
            if (score > scoreTwo)
                workingScore = score;
            else
                workingScore = scoreTwo;
                
            //std::cout << "workingScore: " << workingScore << std::endl;
            
            if (workingScore > bestScore)
            {
                binRLvector.clear();
                binRLvector.push_back(thisBinAveragePosition);
                //binRLvector.push_back(nextBinAveragePosition);
                //binRLvector.push_back(previousBinAveragePosition);
                bestScore = workingScore;
            }
        }
    }
    
    //std::cout << "Number of spikes: " << binRLvector.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::ConvertBinRLToSpikeRL(const std::vector<float> &binRLvector, std::vector<float> &spikeRLvector, const std::vector<CartesianVector> &energyAlongRLvector)
{
    for (const float &binRL : binRLvector)
    {
        float closestApproach(100.f);
        float closestMatchingRL(0.f);
        
        for (const CartesianVector &energyRL : energyAlongRLvector)
        {
            if (std::abs(energyRL.GetX() - binRL) < closestApproach)
            {
                closestApproach = std::abs(energyRL.GetX() - binRL);
                closestMatchingRL = energyRL.GetX();
            }
        }
        
        for (const CartesianVector &energyRL : energyAlongRLvector)
        {
            if (energyRL.GetX() == closestMatchingRL)
            {
                spikeRLvector.push_back(energyRL.GetX());
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::ConvertRLtoCaloHit(const float &spikeRL, std::vector<CartesianVector> &filteredEnergyAlongRLvector, const CaloHitList &caloHitList, CartesianVector &hitPosition)
{
    std::vector<std::pair<float, float>> XtoRLvector;
    
    float targetEnergy(0.f);
    
    for (CartesianVector &energyRL : filteredEnergyAlongRLvector)
    {
        if (energyRL.GetX() == spikeRL)
            targetEnergy = energyRL.GetZ();
    }
    
    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CaloHit *const pCaloHit(*hitIter);
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        const float caloHitEnergy(1000* pCaloHit->GetElectromagneticEnergy()); //MeV
        
        if (caloHitEnergy == targetEnergy)
            hitPosition = caloHitPosition;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindMatchingHitsInDifferentView(const ClusterList &clusterList, CartesianVector &energySpikePosition, std::vector<CartesianVector> &matchedHits)
{
    for (ClusterList::const_iterator iter = clusterList.begin(), iter1End = clusterList.end(); iter != iter1End; ++iter)
    {
        const Cluster *const pCluster = *iter;
        
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.GetCaloHitList(caloHitList);
    
        for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit(*hitIter);
            const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
            
            if (caloHitPosition.GetX() + 0.25 > energySpikePosition.GetX() && caloHitPosition.GetX() - 0.25 < energySpikePosition.GetX())
                matchedHits.push_back(caloHitPosition);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateVerticesFromSpikes(const std::vector<float>energySpikeRLvector, std::vector<CartesianVector> filteredEnergyAlongRLvector, CaloHitList &caloHitList, const ClusterList &clusterList1, const ClusterList &clusterList2, const ClusterList &clusterList3)
{
    for (const float &energySpikeRL : energySpikeRLvector)
    {
        CartesianVector energySpikePosition(0.f, 0.f, 0.f);
        this->ConvertRLtoCaloHit(energySpikeRL, filteredEnergyAlongRLvector, caloHitList, energySpikePosition);
        
        std::vector<CartesianVector> energySpikes1, matchedHits2, matchedHits3;
        
        energySpikes1.push_back(energySpikePosition);
        
        this->FindMatchingHitsInDifferentView(clusterList2, energySpikePosition, matchedHits2);
        this->FindMatchingHitsInDifferentView(clusterList3, energySpikePosition, matchedHits3);
        
        if (!clusterList1.empty() && !clusterList2.empty())
        {
            const Cluster *const pCluster1(*(clusterList1.begin()));
            const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
        
            const Cluster *const pCluster2(*(clusterList2.begin()));
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
            
            this->CreateMatchedVertices(energySpikes1, matchedHits2, hitType1, hitType2);
        }
        
        if (!clusterList1.empty() && !clusterList3.empty())
        {
            const Cluster *const pCluster1(*(clusterList1.begin()));
            const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
            
            const Cluster *const pCluster3(*(clusterList3.begin()));
            const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));
            
            this->CreateMatchedVertices(energySpikes1, matchedHits3, hitType1, hitType3);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &CandidateVertexCreationAlgorithm::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CandidateVertexCreationAlgorithm::SortSpacePointsByZ(CartesianVector &vector1, CartesianVector &vector2)
{
    return vector1.GetZ() < vector2.GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CandidateVertexCreationAlgorithm::SortEnergyVectorByRL(CartesianVector &vector1, CartesianVector &vector2)
{
    return vector1.GetX() < vector2.GetX();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterXDiscrepancy", m_maxClusterXDiscrepancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableCrossingCandidates", m_enableCrossingCandidates));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableEnergyCandidates", m_enableEnergyCandidates));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyPlot", m_energyPlot));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCrossingClusterSize", m_minCrossingClusterSize));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtrapolationLength", m_extrapolationLength));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtrapolationStepSize", m_extrapolationStepSize));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCrossingApproach", m_minClusterCrossingApproach));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PostCrossingSkipDistance", m_postCrossingSkipDistance));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinEnergyVertexClusterSize", m_minEnergyVertexClusterSize));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content