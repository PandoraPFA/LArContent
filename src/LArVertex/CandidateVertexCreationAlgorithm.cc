/**
 *  @file   LArContent/src/LArVertex/CandidateVertexCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the candidate vertex creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArVertex/CandidateVertexCreationAlgorithm.h"

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
    m_minClusterCaloHits(5),
    m_minClusterLengthSquared(3.f * 3.f),
    m_maxClusterXDiscrepancy(4.f),
    m_chiSquaredCut(2.f),
    m_enableCrossingCandidates(false),
    m_enableEnergyCandidates(false),
    m_strictMatching(true),
    m_energyPlot(false),
    m_maxScatterRms(0.4f), //0.2
    m_maxScatterCosTheta(0.905f), //0.905f
    m_maxSlidingCosTheta(0.985f), //0.985f
    m_energyDifferenceThreshold(0.5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::Run()
{
    //std::cout << "Cluster list sizes: " << clusterListU.size() << " " << clusterListV.size() << " " << clusterListW.size() << std::endl;

    try
    {
        ClusterList clusterListU, clusterListV, clusterListW;
        this->SelectClusters(clusterListU, clusterListV, clusterListW);

        const VertexList *pVertexList(NULL); std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

        this->ClusterEndPointComparison(clusterListU, clusterListV);
        this->ClusterEndPointComparison(clusterListU, clusterListW);
        this->ClusterEndPointComparison(clusterListV, clusterListW);
        
        if (m_enableCrossingCandidates)
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

        if (!pVertexList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_topologyVertexListName));

            if (m_replaceCurrentVertexList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_topologyVertexListName));
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }

    this->TidyUp();

    if (m_enableEnergyCandidates)
    {
        try
        {   
            ClusterList energyClusterListU, energyClusterListV, energyClusterListW;
            this->SelectClusters(energyClusterListU, energyClusterListV, energyClusterListW);
            
            VertexList energyVertices;
            std::vector<CartesianVector> energySpikesU, energySpikesV, energySpikesW;

            const VertexList *pEnergyVerticesTemporaryList(NULL);
            std::string energyVerticesTemporaryList;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pEnergyVerticesTemporaryList, energyVerticesTemporaryList));

            //this->FindEnergySpikes(energyClusterListU, energySpikesU);
            //this->FindEnergySpikes(energyClusterListV, energySpikesV);
            //this->FindEnergySpikes(energyClusterListW, energySpikesW);
            
            //LONGITUDINAL
            //this->CreateMatchedVertices(energySpikesU, energySpikesV, TPC_VIEW_U, TPC_VIEW_V);
            //this->CreateMatchedVertices(energySpikesU, energySpikesW, TPC_VIEW_U, TPC_VIEW_W);
            //this->CreateMatchedVertices(energySpikesV, energySpikesW, TPC_VIEW_V, TPC_VIEW_W);
            
            //TRANSVERSE
            this->CreateXBasedVertices(energyClusterListU, energyClusterListV, energyClusterListW);
            this->CreateXBasedVertices(energyClusterListV, energyClusterListW, energyClusterListU);
            this->CreateXBasedVertices(energyClusterListW, energyClusterListU, energyClusterListV);
            
            for (const Vertex *const pVertex : (*pEnergyVerticesTemporaryList))
                energyVertices.insert(pVertex);
                
            if (!pEnergyVerticesTemporaryList->empty())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_energyVertexListName, energyVertices));
        }
        catch (StatusCodeException &statusCodeException)
        {
            this->TidyUp();
        }

        this->TidyUp();
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::SelectClusters(ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW)
{
    this->SelectClusters(m_inputClusterListNameU, clusterListU);
    this->SelectClusters(m_inputClusterListNameV, clusterListV);
    this->SelectClusters(m_inputClusterListNameW, clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::ClusterEndPointComparison(const ClusterList &clusterList1, const ClusterList &clusterList2) const
{
    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iter1End = clusterList1.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pCluster1(*iter1);
        const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));

        const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(pCluster1));
        const CartesianVector minLayerPosition1(fitResult1.GetGlobalMinLayerPosition());
        const CartesianVector maxLayerPosition1(fitResult1.GetGlobalMaxLayerPosition());

        for (ClusterList::const_iterator iter2 = clusterList2.begin(), iter2End = clusterList2.end(); iter2 != iter2End; ++iter2)
        {
            const Cluster *const pCluster2(*iter2);
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

            const TwoDSlidingFitResult &fitResult2(this->GetCachedSlidingFitResult(*iter2));
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

void CandidateVertexCreationAlgorithm::SelectClusters(const std::string &clusterListName, ClusterList &selectedClusterList)
{
    const ClusterList *pInputClusterList(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pInputClusterList));

    if (!pInputClusterList || pInputClusterList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CandidateVertexCreationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

        return;
    }

    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            const Cluster *const pCluster = *iter;

            if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
                continue;

            if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
                continue;

            this->AddToSlidingFitCache(pCluster);
            selectedClusterList.insert(pCluster);
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::Find2DClusterCrossings(const ClusterList &clusterList, std::vector<CartesianVector> &crossingsVector) const
{
for (ClusterList::const_iterator iter1 = clusterList.begin(), iter1End = clusterList.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pCluster1 = *iter1;
        
        if (pCluster1->GetNCaloHits() < 10)
            continue;
        
        const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(*iter1));
        const CartesianVector minLayerPosition1(fitResult1.GetGlobalMinLayerPosition());
        const CartesianVector maxLayerPosition1(fitResult1.GetGlobalMaxLayerPosition());
        
        crossingsVector.push_back(minLayerPosition1);
        crossingsVector.push_back(maxLayerPosition1);

        std::vector<CartesianVector> spacePointVector1;
        
        OrderedCaloHitList orderedCaloHitList1(pCluster1->GetOrderedCaloHitList());
        CaloHitList caloHitList1;
        orderedCaloHitList1.GetCaloHitList(caloHitList1);
        
        for (CaloHitList::const_iterator caloIter1 = caloHitList1.begin(), caloIterEnd1 = caloHitList1.end(); caloIter1 != caloIterEnd1; ++caloIter1)
        {
            const CaloHit *const pCaloHit = (*caloIter1);
            CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
            spacePointVector1.push_back(caloHitPosition);
        }
        
        for (float i = 0.1; i < 20.0; i += 0.1)
        {
            CartesianVector tempExtrapolatedPositionUnder(0.f, 0.f, 0.f);
            CartesianVector tempExtrapolatedPositionOver(0.f, 0.f, 0.f);

            float minLayerRL(fitResult1.GetL(fitResult1.GetMinLayer()));
            float maxLayerRL(fitResult1.GetL(fitResult1.GetMaxLayer()));

            fitResult1.GetExtrapolatedPosition(minLayerRL - i, tempExtrapolatedPositionUnder);
            fitResult1.GetExtrapolatedPosition(maxLayerRL + i, tempExtrapolatedPositionOver);
            
            //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tempExtrapolatedPositionUnder, "Target Vertex", BLUE, 1));
            //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tempExtrapolatedPositionOver, "Target Vertex", BLUE, 1));

            spacePointVector1.push_back(tempExtrapolatedPositionUnder);
            spacePointVector1.push_back(tempExtrapolatedPositionOver);
        }

        for (ClusterList::const_iterator iter2 = clusterList.begin(), iter2End = clusterList.end(); iter2 != iter2End; ++iter2)
        {
            if ((*iter1)->GetNCaloHits() == (*iter2)->GetNCaloHits())
                continue;
                
            const Cluster *const pCluster2 = *iter2;
            
            if (pCluster2->GetNCaloHits() < 10)
            continue;
            
            const TwoDSlidingFitResult &fitResult2(this->GetCachedSlidingFitResult(*iter2));
            const CartesianVector minLayerPosition2(fitResult2.GetGlobalMinLayerPosition());
            const CartesianVector maxLayerPosition2(fitResult2.GetGlobalMaxLayerPosition());

            std::vector<CartesianVector> spacePointVector2;
            
            OrderedCaloHitList orderedCaloHitList2(pCluster2->GetOrderedCaloHitList());
            CaloHitList caloHitList2;
            orderedCaloHitList2.GetCaloHitList(caloHitList2);
            
            for (CaloHitList::const_iterator caloIter2 = caloHitList2.begin(), caloIterEnd2 = caloHitList2.end(); caloIter2 != caloIterEnd2; ++caloIter2)
            {
                const CaloHit *const pCaloHit = (*caloIter2);
                CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
                spacePointVector2.push_back(caloHitPosition);
            }

            for (float i = 0.1; i < 20.0; i += 0.1)
            {
                CartesianVector tempExtrapolatedPositionUnder(0.f, 0.f, 0.f);
                CartesianVector tempExtrapolatedPositionOver(0.f, 0.f, 0.f);

                float minLayerRL(fitResult2.GetL(fitResult2.GetMinLayer()));
                float maxLayerRL(fitResult2.GetL(fitResult2.GetMaxLayer()));

                fitResult2.GetExtrapolatedPosition(minLayerRL - i, tempExtrapolatedPositionUnder);
                fitResult2.GetExtrapolatedPosition(maxLayerRL + i, tempExtrapolatedPositionOver);
                
                //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tempExtrapolatedPositionUnder, "Extrapolated Point", BLUE, 1));
                //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tempExtrapolatedPositionOver, "Extrapolated Point", BLUE, 1));

                spacePointVector2.push_back(tempExtrapolatedPositionUnder);
                spacePointVector2.push_back(tempExtrapolatedPositionOver);
            }
            
            std::vector<float> distancesBetweenClusters;
            
            for (CartesianVector &position1: spacePointVector1)
            {
                for (CartesianVector &position2: spacePointVector2)
                    distancesBetweenClusters.push_back((position1-position2).GetMagnitude());
            }
            
            std::sort(spacePointVector1.begin(), spacePointVector1.end(), SortSpacePointsByZ);
            std::sort(spacePointVector2.begin(), spacePointVector2.end(), SortSpacePointsByZ);
            
            int skipCounter(0);
            bool shouldSkip(false);
            
            for (CartesianVector &position1: spacePointVector1)
            {
                if (shouldSkip)
                {
                    if (skipCounter <= 50)
                    {
                        skipCounter++;
                        continue;
                    }
                }
                
                skipCounter = 0;
                shouldSkip = false;
                
                for (CartesianVector &position2: spacePointVector2)
                {
                    if ((position1-position2).GetMagnitude() < 0.25)
                    {
                        crossingsVector.push_back(position1);
                        crossingsVector.push_back(position2);
                        
                        shouldSkip = true;
                        break;
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindEnergySpikes(const ClusterList &clusterList, std::vector<CartesianVector> &energySpikesVector)
{
    for (ClusterList::const_iterator iter1 = clusterList.begin(), iter1End = clusterList.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pCluster = *iter1;
        
        if (pCluster->GetNCaloHits() < 60)
            continue;
    
        const TwoDSlidingFitResult &slidingFitResult(this->GetCachedSlidingFitResult(*iter1));
        
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.GetCaloHitList(caloHitList);
        
        std::vector<CartesianVector> energyAlongRLvector;
        this->CreateEnergyAlongRLVector(slidingFitResult, caloHitList, energyAlongRLvector);
        
        std::sort(energyAlongRLvector.begin(), energyAlongRLvector.end(), SortEnergyVectorByRL);
        
        std::vector<CartesianVector> filteredEnergyAlongRLvector;
        this->FilterEnergyVector(energyAlongRLvector, filteredEnergyAlongRLvector);
        this->FilterEnergyVector(filteredEnergyAlongRLvector, filteredEnergyAlongRLvector);
        this->FilterEnergyVector(filteredEnergyAlongRLvector, filteredEnergyAlongRLvector);
        
        std::sort(filteredEnergyAlongRLvector.begin(), filteredEnergyAlongRLvector.end(), SortEnergyVectorByRL);
        
        //const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        //const TwoDSlidingFitResult energySlidingFitResult(filteredEnergyAlongRLvector, 2*m_slidingFitWindow, slidingFitPitch); 
        
        bool spikeFound(false);
        std::vector<float> energySpikeRLvector;
        
        this->FindEnergySpike(filteredEnergyAlongRLvector, energySpikeRLvector);
        
        if (!spikeFound)
            continue;
        
        if (m_energyPlot)
            this->DrawEnergyVector(filteredEnergyAlongRLvector, pCluster);
            
        for (const float &energySpikeRL : energySpikeRLvector)
        {
            CartesianVector energySpikePosition(0.f, 0.f, 0.f);
            this->ConvertRLtoCaloHit(energySpikeRL, filteredEnergyAlongRLvector, caloHitList, energySpikePosition);
            energySpikesVector.push_back(energySpikePosition);
        }
        
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateXBasedVertices(const ClusterList &clusterList1, const ClusterList &clusterList2, const ClusterList &clusterList3)
{
    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iter1End = clusterList1.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pCluster = *iter1;
        
        if (pCluster->GetNCaloHits() < 60)
            continue;
    
        const TwoDSlidingFitResult &slidingFitResult(this->GetCachedSlidingFitResult(*iter1));
        
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.GetCaloHitList(caloHitList);
        
        std::vector<CartesianVector> energyAlongRLvector;
        this->CreateEnergyAlongRLVector(slidingFitResult, caloHitList, energyAlongRLvector);
        
        std::vector<CartesianVector> filteredEnergyAlongRLvector;
        this->FilterEnergyVector(energyAlongRLvector, filteredEnergyAlongRLvector);
        
        std::vector<float> energySpikeRLvector;
        this->FindEnergySpike(filteredEnergyAlongRLvector, energySpikeRLvector);
        
        if (energySpikeRLvector.empty())
            std::cout << "No spike found." << std::endl;
        
        this->CreateVerticesFromSpikes(energySpikeRLvector, filteredEnergyAlongRLvector, caloHitList, clusterList1, clusterList2, clusterList3);
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

bool CandidateVertexCreationAlgorithm::SortEnergyVectorByRL(CartesianVector &position1, CartesianVector &position2)
{
    return position1.GetX() < position2.GetX();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateEnergyAlongRLVector(const TwoDSlidingFitResult &slidingFitResult, const CaloHitList &caloHitList, std::vector<CartesianVector> &energyAlongRLvector)
{
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

void CandidateVertexCreationAlgorithm::DrawEnergyVector(std::vector<CartesianVector> &energyAlongRLvector, const Cluster* pCluster)
{
    std::sort(energyAlongRLvector.begin(), energyAlongRLvector.end(), SortEnergyVectorByRL);
    
    TGraph *HitEnergy_vs_rL = new TGraph(energyAlongRLvector.size());
    ClusterList tempClusterList;
    tempClusterList.insert(pCluster);
    
    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &tempClusterList, "Visualised Cluster", BLUE, 1));
    
    int n(0);
    
    for (std::vector<CartesianVector>::const_iterator pairIter = energyAlongRLvector.begin(), pairIterEnd = std::prev(energyAlongRLvector.end(), 4); pairIter != pairIterEnd; ++pairIter)
    {
        HitEnergy_vs_rL->SetPoint(n, (*pairIter).GetX(), (*pairIter).GetZ());
        n++;
    }
    
    TCanvas *canvas1 = new TCanvas("HitEnergy_vs_rL", "HitEnergy_vs_rL", 900, 600);
    canvas1->cd();
    HitEnergy_vs_rL->SetMarkerStyle(6);
    HitEnergy_vs_rL->Draw("AP");
    PANDORA_MONITORING_API(Pause(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FilterEnergyVector(const std::vector<CartesianVector> &unfilteredEnergyVector, std::vector<CartesianVector> &filteredEnergyVector)
{
    for (std::vector<CartesianVector>::const_iterator pairIter = std::next(unfilteredEnergyVector.begin(), 1), pairIterEnd = std::prev(unfilteredEnergyVector.end(), 1); pairIter != pairIterEnd; ++pairIter)
    {
        float chargeRatioNext(((*(std::next(pairIter, 1))).GetZ())/((*pairIter).GetZ()));
        float chargeRatioPrevious(((*(std::prev(pairIter, 1))).GetZ())/((*pairIter).GetZ()));
        //float chargeDifference(std::abs(((*(std::prev(pairIter, 1))).GetZ()) - ((*pairIter).GetZ())));
        
        //std::cout << "(X, Y, Z): " << (*pairIter).GetX() << ", " << (*pairIter).GetY() << ", " << (*pairIter).GetZ() << std::endl;
        
        if (!(chargeRatioPrevious < 0.8 && chargeRatioNext < 0.8))
            filteredEnergyVector.push_back(*pairIter);
    }
    
    std::sort(filteredEnergyVector.begin(), filteredEnergyVector.end(), SortEnergyVectorByRL);
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

void CandidateVertexCreationAlgorithm::FindEnergySpike(const std::vector<CartesianVector> &energyAlongRLvector, std::vector<float> &spikeRLvector)
{
    std::vector<CartesianVector> binnedEnergyAlongRLvector;
    this->BinEnergyRLVector(energyAlongRLvector, binnedEnergyAlongRLvector);
    
    std::vector<float> binRLvector;
    this->FindBinRLFromBinnedVector(binnedEnergyAlongRLvector, binRLvector);
    
    //DRAW BINNED DISTRIBUTION
    //TGraph *HitEnergy_vs_rL = new TGraph(binnedEnergyAlongRLvector.size());
    //int n(0);    
    //for (std::vector<CartesianVector>::const_iterator pairIter = binnedEnergyAlongRLvector.begin(), pairIterEnd = binnedEnergyAlongRLvector.end(); pairIter != pairIterEnd; ++pairIter)
    //{
    //    HitEnergy_vs_rL->SetPoint(n, (*pairIter).GetX(), (*pairIter).GetZ());
    //    n++;
    //}    
    //TCanvas *canvas1 = new TCanvas("HitEnergy_vs_rL", "HitEnergy_vs_rL", 900, 600);
    //canvas1->cd();
    //HitEnergy_vs_rL->SetMarkerStyle(6);
    //HitEnergy_vs_rL->Draw("AP");
    //PANDORA_MONITORING_API(Pause(this->GetPandora()));
    
    this->ConvertBinRLToSpikeRL(binRLvector, spikeRLvector, energyAlongRLvector);
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

void CandidateVertexCreationAlgorithm::FindBinRLFromBinnedVector(const std::vector<CartesianVector> &binnedEnergyAlongRLvector, std::vector<float> &binRLvector)
{
    for (std::vector<CartesianVector>::const_iterator pairIter = std::next(binnedEnergyAlongRLvector.begin(), 1), pairIterEnd = std::prev(binnedEnergyAlongRLvector.end(), 2); pairIter != pairIterEnd; ++pairIter)
    {
        float thisBinAveragePosition((*pairIter).GetX());
        float thisBinAverageEnergy((*pairIter).GetZ());
        
        float previousBinAverageEnergy((*std::prev(pairIter, 1)).GetZ());
        float nextBinAverageEnergy((*std::next(pairIter, 1)).GetZ());
        
        float previousPreviousBinAverageEnergy((*std::prev(pairIter, 2)).GetZ());
        float nextNextBinAverageEnergy((*std::next(pairIter, 2)).GetZ());
        
        //std::cout << thisBinAveragePosition << std::endl;
        //std::cout << thisBinAverageEnergy << std::endl;
        //std::cout << std::abs(nextBinAverageEnergy - thisBinAverageEnergy) << std::endl;
        //std::cout << "******************" << std::endl;
        

            //std::cout << "Jump position: " << thisBinAveragePosition << std::endl;
            //std::cout << "Jump ratio: " << std::abs(1 - std::abs(nextBinAverageEnergy / thisBinAverageEnergy)) << std::endl;
            //std::cout << "Next jump ratio: " << std::abs(1 - (std::abs(nextNextBinAverageEnergy / nextBinAverageEnergy))) << std::endl;
            //std::cout << "Previous jump ratio: " << std::abs(1 - std::abs(previousBinAverageEnergy / thisBinAverageEnergy)) << std::endl;
            //std::cout << "Previous previous jump ratio: " << std::abs(1 - std::abs(previousPreviousBinAverageEnergy / previousBinAverageEnergy)) << std::endl;
            //std::cout << "*******************" << std::endl;
            
            //This if statement means: can we find two bins close to one another, then a jump and then again two bins that are close to one another, in either the + or - RL direction?
            if ((std::abs(1 - std::abs(nextBinAverageEnergy / thisBinAverageEnergy)) > 0.3 && std::abs(1 - (std::abs(nextNextBinAverageEnergy / nextBinAverageEnergy))) > 0.3
            && std::abs(1 - (std::abs(previousBinAverageEnergy / thisBinAverageEnergy))) < 0.15 && std::abs(1 - (std::abs(previousPreviousBinAverageEnergy / previousBinAverageEnergy))) < 0.15)
            || (std::abs(1 - std::abs(previousBinAverageEnergy / thisBinAverageEnergy)) > 0.3 && std::abs(1 - (std::abs(previousPreviousBinAverageEnergy / previousBinAverageEnergy))) > 0.3
            && std::abs(1 - (std::abs(nextBinAverageEnergy / thisBinAverageEnergy))) < 0.15 && std::abs(1 - (std::abs(nextNextBinAverageEnergy / nextBinAverageEnergy))) < 0.15))
            {
                if (thisBinAverageEnergy < 0.25 || thisBinAverageEnergy == 4.375)
                    continue;
                
                binRLvector.push_back(thisBinAveragePosition);
                //std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>JUMP" << std::endl;

                //binRLvector.push_back(nextBinAveragePosition);
                //binRLvector.push_back(previousBinAveragePosition);
            }
    }
    
    //std::cout << "Number of spikes: " << binRLvector.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::ConvertBinRLToSpikeRL(const std::vector<float> &binRLvector, std::vector<float> &spikeRLvector, const std::vector<CartesianVector> &energyAlongRLvector)
{
    for (const float &binRL : binRLvector)
    {
        //std::cout << "bin RL: " << binRL << std::endl;
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
                //std::cout << "spike position: " << energyRL.GetX() << std::endl;
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

void CandidateVertexCreationAlgorithm::CreateMatchedVertices(std::vector<CartesianVector> &crossingsVector1, std::vector<CartesianVector> &crossingsVector2, HitType hitType1, HitType hitType2) const
{
     if (crossingsVector1.empty() || crossingsVector2.empty())
        return;
    
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
            
            if (m_strictMatching)
            {
                std::pair<CartesianVector*,float> positionChiSquaredPair;
                positionChiSquaredPair = std::make_pair(&position3D, chiSquared);
                
                matched3DPositions.push_back(positionChiSquaredPair);
                chiSquareds.push_back(chiSquared);
            }
            else
            {
                PandoraContentApi::Vertex::Parameters parameters;
                parameters.m_position = position3D;
                parameters.m_vertexLabel = VERTEX_INTERACTION;
                parameters.m_vertexType = VERTEX_3D;
                        
                const Vertex *pVertex(NULL);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
                
                //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionU, "New Vertex", GREEN, 1));
                //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionV, "New Vertex", GREEN, 1));
                //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionW, "New Vertex", GREEN, 1));
                //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }
            
            
        }
        
        if (m_strictMatching)
        {
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
                    //std::cout << ">>>>>> Created vertex." << std::endl;
                    break;
                }
            }
        }
    }
    
    //std::cout << "***********************" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TopologyVertexListName", m_topologyVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "EnergyVertexListName", m_energyVertexListName));

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
        "StrictMatching", m_strictMatching));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyPlot", m_energyPlot));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyDifferenceThreshold", m_energyDifferenceThreshold));
        

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
