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

using namespace pandora;

namespace lar_content
{

CandidateVertexCreationAlgorithm::CandidateVertexCreationAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_slidingFitWindow(20),
    m_minClusterCaloHits(5),
    m_minClusterLengthSquared(3.f * 3.f),
    m_maxClusterXDiscrepancy(4.f),
    m_chiSquaredCut(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::Run()
{
    try
    {
        ClusterList clusterListU, clusterListV, clusterListW;
        this->SelectClusters(clusterListU, clusterListV, clusterListW);

        const VertexList *pVertexList(NULL); std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

        this->ClusterEndPointComparison(clusterListU, clusterListV);
        this->ClusterEndPointComparison(clusterListU, clusterListW);
        this->ClusterEndPointComparison(clusterListV, clusterListW);
        
        //std::vector<CartesianVector> crossingsU;
        //std::vector<CartesianVector> crossingsV;
        //std::vector<CartesianVector> crossingsW;
        //
        //this->Find2DClusterCrossings(clusterListU, crossingsU);
        //this->Find2DClusterCrossings(clusterListV, crossingsV);
        //this->Find2DClusterCrossings(clusterListW, crossingsW);
        //
        //this->CreateCrossingVertex(crossingsU, crossingsV, TPC_VIEW_U, TPC_VIEW_V);
        //this->CreateCrossingVertex(crossingsU, crossingsW, TPC_VIEW_U, TPC_VIEW_W);
        //this->CreateCrossingVertex(crossingsV, crossingsW, TPC_VIEW_V, TPC_VIEW_W);

        if (!pVertexList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

            if (m_replaceCurrentVertexList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
        }

        this->TidyUp();
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
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
            
            bool thisClusterCrossingHasBeenFound(false);
            
            std::vector<float> distancesBetweenClusters;
            
            for (CartesianVector &position1: spacePointVector1)
            {
                for (CartesianVector &position2: spacePointVector2)
                    distancesBetweenClusters.push_back((position1-position2).GetMagnitude());
            }

            for (CartesianVector &position1: spacePointVector1)
            {
                for (CartesianVector &position2: spacePointVector2)
                {
                    //if ((position1-position2).GetMagnitude() + 0.1 > (*(std::min_element(distancesBetweenClusters.begin(), distancesBetweenClusters.end())))
                    //&& (position1-position2).GetMagnitude() - 0.1 < (*(std::min_element(distancesBetweenClusters.begin(), distancesBetweenClusters.end()))))
                    if ((position1-position2).GetMagnitude() < 0.5)
                    {
                        crossingsVector.push_back(position1);
                        crossingsVector.push_back(position2);
                        
                        //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position1, "Crossover Point", RED, 1));
                        //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position2, "Crossover Point", RED, 1));
                        
                        thisClusterCrossingHasBeenFound = true;
                        break;
                    }
                }
                
                if (thisClusterCrossingHasBeenFound)
                    break;
            }
            
            //CartesianVector minimalApproachPointOne(0.f, 0.f, 0.f);
            //CartesianVector minimalApproachPointTwo(0.f, 0.f, 0.f);
            //float minimalApproach(100.f);
            //
            //for (CartesianVector &position1: spacePointVector1)
            //{
            //    for (CartesianVector &position2: spacePointVector2)
            //    {
            //        if ((position1-position2).GetMagnitude() > 0.5)
            //            continue;
            //        
            //        if ((position1-position2).GetMagnitude() < minimalApproach)
            //        {
            //            minimalApproachPointOne  = position1;
            //            minimalApproachPointTwo = position2;
            //            minimalApproach = (position1-position2).GetMagnitude();
            //        }
            //    }
            //    
            //    crossingsVector.push_back(minimalApproachPointOne);
            //    crossingsVector.push_back(minimalApproachPointTwo);
            //    //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &minimalApproachPointOne, "Target Vertex", RED, 1));
            //    //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &minimalApproachPointTwo, "Target Vertex", RED, 1));
            //}
            
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateCrossingVertex(std::vector<CartesianVector> &crossingsVector1, std::vector<CartesianVector> &crossingsVector2, HitType hitType1, HitType hitType2) const
{
    std::vector<CartesianVector> spacePointVector;
    
    for (CartesianVector &position1: crossingsVector1)
    {
        for (CartesianVector &position2: crossingsVector2)
        {
            if (std::fabs(position1.GetX() - position2.GetX()) > 0.1)
                continue;
                
            float chiSquared(0.f);
            CartesianVector position3D(0.f, 0.f, 0.f);
            LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, position1, position2, position3D, chiSquared);
            const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_W));
                
            if (chiSquared > 2.0)
                return;
                
            //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionW, "New Vertex", GREEN, 1));
                
            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = position3D;
            parameters.m_vertexLabel = VERTEX_INTERACTION;
            parameters.m_vertexType = VERTEX_3D;
                
            const Vertex *pVertex(NULL);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
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

StatusCode CandidateVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

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

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
