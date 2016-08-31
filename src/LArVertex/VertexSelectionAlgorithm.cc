/**
 *  @file   LArContent/src/LArVertex/VertexSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex selection algorithm class.
 * 
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArUtility/KDTreeLinkerAlgoT.h"

#include "LArVertex/VertexSelectionAlgorithm.h"
#include "LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

VertexSelectionAlgorithm::VertexSelectionAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_nDecayLengthsInZSpan(2.f),
    m_kappa(0.42f),
    m_selectSingleVertex(true),
    m_maxTopScoreSelections(3),
    m_minCandidateDisplacement(2.f),
    m_minCandidateScoreFraction(0.5f),
    m_useDetectorGaps(true),
    m_gapTolerance(0.f),
    m_isEmptyViewAcceptable(true),
    m_enableClustering(false),
    m_directionFilter(false),
    m_beamWeightFilter(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::Run()
{
    const VertexList *pTopologyVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_topologyVertexListName, pTopologyVertexList));
    
    if (!pTopologyVertexList || pTopologyVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexSelectionAlgorithm: unable to find current vertex list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }
    
    std::vector<const VertexList*> vertexListVector = m_pVertexClusteringTool->ClusterVertices(pTopologyVertexList);
    
    int nVertices(0);
    
    for (const VertexList* pVertexList : vertexListVector)
    {
        nVertices += pVertexList->size();
        
        //for (const Vertex *const pVertex : (*pVertexList))
        //{
        //    CartesianVector vertexPosition(pVertex->GetPosition());
        //    
        //    const CartesianVector vertexProjectionU(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_U));
        //    const CartesianVector vertexProjectionV(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_V));
        //    const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_W));
        //    
        //    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionU, "Top 5 Vertex U", RED, 1));
        //    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionV, "Top 5 Vertex V", RED, 1));
        //    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionW, "Top 5 Vertex W", RED, 1));
        //}        
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
    
    std::vector<VertexScoringTool::VertexScoreList> scoredClusterCollection;
    m_pVertexScoringTool->ScoreVertices(this, pTopologyVertexList, vertexListVector, scoredClusterCollection);
    
    VertexList selectedVertexList;
    this->SelectTopScoreVertices(scoredClusterCollection, selectedVertexList);
    
    if (!selectedVertexList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, selectedVertexList));
    
    //---------------------------------------------------------------------------------------------------------------------------------------
    
    const VertexList *pEnergyVertexList(NULL);
    VertexScoringTool::VertexScoreList energyVertexScoreList;
    
    try
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_energyVertexListName, pEnergyVertexList));
        m_pVertexScoringTool->ScoreEnergyVertices(this, pEnergyVertexList, energyVertexScoreList);
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "There are no energy vertices present." << std::endl;
    }

    //---------------------------------------------------------------------------------------------------------------------------------------
    
    this->StoreTopAllInformation(pTopologyVertexList, selectedVertexList, pEnergyVertexList);
    
    VertexScoringTool::VertexScoreList top5VertexScoreList;
    this->FindTop5Vertices(scoredClusterCollection, energyVertexScoreList, top5VertexScoreList);
    
    std::cout << "Before filtering: " << top5VertexScoreList.size() << std::endl;
    
    VertexScoringTool::VertexScoreList filteredTop5VertexScoreList;
    this->FilterTop5Vertices(top5VertexScoreList, filteredTop5VertexScoreList);
    
    std::cout << "After filtering: " << filteredTop5VertexScoreList.size() << std::endl;
    
    this->StoreTop5Information(filteredTop5VertexScoreList);
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    if (!selectedVertexList.empty())
    {
        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }
    
    return STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectTopScoreVertices(std::vector<VertexScoringTool::VertexScoreList> &scoredClusterCollection, VertexList &selectedVertexList) const
{
    float bestScore(0.f);
    
    VertexScoringTool::VertexScoreList vertexScoreList;
    for (VertexScoringTool::VertexScoreList &tempVertexScoreList : scoredClusterCollection)
        vertexScoreList.insert(vertexScoreList.begin(), tempVertexScoreList.begin(), tempVertexScoreList.end());
    
    std::sort(vertexScoreList.begin(), vertexScoreList.end());
    
    for (const VertexScoringTool::VertexScore &vertexScore : vertexScoreList)
    {
        if (selectedVertexList.size() >= m_maxTopScoreSelections)
            break;

        if (!selectedVertexList.empty() && !this->AcceptVertexLocation(vertexScore.GetVertex(), selectedVertexList))
            continue;

        if (!selectedVertexList.empty() && (vertexScore.GetScore() < m_minCandidateScoreFraction * bestScore))
            continue;

        selectedVertexList.insert(vertexScore.GetVertex());

        if (m_selectSingleVertex)
            return;

        if (vertexScore.GetScore() > bestScore)
            bestScore = vertexScore.GetScore();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexLocation(const Vertex *const pVertex, const VertexList &selectedVertexList) const
{
    const CartesianVector &position(pVertex->GetPosition());
    const float minCandidateDisplacementSquared(m_minCandidateDisplacement * m_minCandidateDisplacement);

    for (const Vertex *const pSelectedVertex : selectedVertexList)
    {
        if (pVertex == pSelectedVertex)
            return false;

        if ((position - pSelectedVertex->GetPosition()).GetMagnitudeSquared() < minCandidateDisplacementSquared)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FindTop5Vertices(std::vector<VertexScoringTool::VertexScoreList> &scoredClusterCollection, VertexScoringTool::VertexScoreList &energyVertexScoreList, VertexScoringTool::VertexScoreList &top5VertexScoreList)
{
    std::sort(energyVertexScoreList.begin(), energyVertexScoreList.end());

    unsigned int clusterCounter(0);
    
    if (m_enableClustering)
    {
        if (!energyVertexScoreList.empty())
        {
            top5VertexScoreList.push_back(*energyVertexScoreList.begin());
            clusterCounter++;
        }
        
        for (VertexScoringTool::VertexScoreList &vertexScoreList : scoredClusterCollection)
        {
            if (clusterCounter == 5)
                break;
                
            if (vertexScoreList.size() == 0)
                continue;
            
            std::sort(vertexScoreList.begin(), vertexScoreList.end());
            
            top5VertexScoreList.push_back(*vertexScoreList.begin());
            clusterCounter++;
        }
    }
    else
    {
        VertexScoringTool::VertexScoreList temporaryVertexScoreList;
        
        for (VertexScoringTool::VertexScoreList &vertexScoreList : scoredClusterCollection)
        {
            if (vertexScoreList.size() == 0)
                continue;
            
            temporaryVertexScoreList.insert(temporaryVertexScoreList.begin(), vertexScoreList.begin(), vertexScoreList.end());
        }
        
        std::sort(temporaryVertexScoreList.begin(), temporaryVertexScoreList.end());
        
        if (!energyVertexScoreList.empty())
        {
            top5VertexScoreList.push_back(*energyVertexScoreList.begin());
            clusterCounter++;
        }
        
        for (VertexScoringTool::VertexScore &vertexScore : temporaryVertexScoreList)
        {
            if (clusterCounter == 5)
                break;
            
            top5VertexScoreList.push_back(vertexScore);
            clusterCounter++;
        }
        
    }
    
    std::sort(top5VertexScoreList.begin(), top5VertexScoreList.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::StoreTop5Information(VertexScoringTool::VertexScoreList &top5VertexScoreList)
{
    //for (VertexScoringTool::VertexScore &vertexScore : top5VertexScoreList)
    //{
    //    CartesianVector vertexPosition(vertexScore.GetVertex()->GetPosition());
    //    
    //    const CartesianVector vertexProjectionU(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_U));
    //    const CartesianVector vertexProjectionV(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_V));
    //    const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_W));
    //    
    //    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionU, "Top 5 Vertex U", RED, 1));
    //    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionV, "Top 5 Vertex V", RED, 1));
    //    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertexProjectionW, "Top 5 Vertex W", RED, 1));
    //}
    //
    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    
    const VertexList *pTop5TemporaryList(NULL);
    std::string top5TemporaryListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTop5TemporaryList, top5TemporaryListName));
    
    VertexList top5VerticesList;
    
    for (const VertexScoringTool::VertexScore &vertexScore: top5VertexScoreList)
    {
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = vertexScore.GetVertex()->GetPosition();
        parameters.m_vertexLabel = vertexScore.GetVertex()->GetVertexLabel();
        parameters.m_vertexType = vertexScore.GetVertex()->GetVertexType();

        const Vertex *pVertexClone(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertexClone));
        top5VerticesList.insert(pVertexClone);
    }
    
    if (!top5VerticesList.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_top5VertexListName, top5VerticesList));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::StoreTopAllInformation(const VertexList* pTopologyVertexList, VertexList selectedVertexList, const VertexList* pEnergyVertexList)
{
    const VertexList *pAllVerticesTemporaryList(NULL);
    std::string allVerticesTemporaryList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pAllVerticesTemporaryList, allVerticesTemporaryList));
    
    VertexList allVerticesList;
    
    if (pTopologyVertexList != NULL)
    {
        for (const Vertex *const pVertex: (*pTopologyVertexList))
        {
            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = pVertex->GetPosition();
            parameters.m_vertexLabel = pVertex->GetVertexLabel();
            parameters.m_vertexType = pVertex->GetVertexType();
        
            const Vertex *pVertexClone(NULL);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertexClone));
            
            allVerticesList.insert(pVertexClone);
        }
    }
        
    if (!selectedVertexList.empty())
    {
        for (const Vertex *const pVertex : selectedVertexList)
        {
            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = pVertex->GetPosition();
            parameters.m_vertexLabel = pVertex->GetVertexLabel();
            parameters.m_vertexType = pVertex->GetVertexType();
    
            const Vertex *pVertexClone(NULL);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertexClone));
            
            allVerticesList.insert(pVertexClone);
        }
    }
    
    if (pEnergyVertexList != NULL)
    {
        for (const Vertex *const pVertex: (*pEnergyVertexList))
        {
            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = pVertex->GetPosition();
            parameters.m_vertexLabel = pVertex->GetVertexLabel();
            parameters.m_vertexType = pVertex->GetVertexType();
        
            const Vertex *pVertexClone(NULL);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertexClone));
            
            allVerticesList.insert(pVertexClone);
        }
    }
    
    if (!allVerticesList.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_allOtherVertexListName, allVerticesList));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FilterTop5Vertices(VertexScoringTool::VertexScoreList &top5VertexScoreList, VertexScoringTool::VertexScoreList &filteredTop5VertexScoreList)
{
    ClusterList clusterListU, clusterListV, clusterListW;
    
    try
    {
        this->GetClusters(clusterListU, clusterListV, clusterListW);
    }
    catch (StatusCodeException &statusCodeException)
    {
        filteredTop5VertexScoreList = top5VertexScoreList;
        return;
    }
    
    const Cluster *const pLongestClusterU(this->GetLongestCluster(clusterListU));
    const Cluster *const pLongestClusterV(this->GetLongestCluster(clusterListV));
    const Cluster *const pLongestClusterW(this->GetLongestCluster(clusterListW));
    
    try
    {
        const MCParticle *const mcParticleU(MCParticleHelper::GetMainMCParticle(pLongestClusterU));
        const MCParticle *const mcParticleV(MCParticleHelper::GetMainMCParticle(pLongestClusterV));
        const MCParticle *const mcParticleW(MCParticleHelper::GetMainMCParticle(pLongestClusterW));
    
        CartesianVector mcBeginPointU(mcParticleU->GetVertex()), mcEndPointU(mcParticleU->GetEndpoint()), mcBeginPointV(mcParticleV->GetVertex()), mcEndPointV(mcParticleV->GetEndpoint()), 
        mcBeginPointW(mcParticleW->GetVertex()), mcEndPointW(mcParticleW->GetEndpoint());
    
        if (m_directionFilter)
        {
            std::cout << "Employing direction filter." << std::endl;
            
            for (VertexScoringTool::VertexScore &vertexScore : top5VertexScoreList)
            {
                const Vertex *const pVertex(vertexScore.GetVertex());
                CartesianVector vertexPosition(pVertex->GetPosition());
                
                const CartesianVector vertexProjectionU(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_U));
                const CartesianVector vertexProjectionV(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_V));
                const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_W));
                
                int nBadViews(0);
                
                if ((vertexProjectionU - mcBeginPointU).GetMagnitude() > (vertexProjectionU - mcEndPointU).GetMagnitude())
                    nBadViews++;
                    
                if ((vertexProjectionV - mcBeginPointV).GetMagnitude() > (vertexProjectionV - mcEndPointV).GetMagnitude())
                    nBadViews++;
                    
                if ((vertexProjectionW - mcBeginPointW).GetMagnitude() > (vertexProjectionW - mcEndPointW).GetMagnitude())
                    nBadViews++;
                    
                if (nBadViews >= 2)
                {
                    std::cout << ">>> FILTERED" << std::endl;
                    continue;
                }
                    
                filteredTop5VertexScoreList.push_back(vertexScore);
            }
        }
        
        if (m_beamWeightFilter)
        {
            std::cout << "Employing beam weight filter." << std::endl;
            
            float smallestZinU(10000.f), largestZinU(0.f);
            for (const Cluster *const pCluster : clusterListU)
            {
                float smallestXCluster(0.f), largestXCluster(0.f), smallestZCluster(0.f), largestZCluster(0.f);
                LArClusterHelper::GetClusterSpanX(pCluster, smallestXCluster, largestXCluster);
                LArClusterHelper::GetClusterSpanZ(pCluster, smallestXCluster, largestXCluster, smallestZCluster, largestZCluster);
                smallestZinU = std::min(smallestZinU, smallestZCluster);
                largestZinU = std::max(largestZinU, largestZCluster);
            }
            
            float smallestZinV(10000.f), largestZinV(0.f);
            for (const Cluster *const pCluster : clusterListV)
            {
                float smallestXCluster(0.f), largestXCluster(0.f), smallestZCluster(0.f), largestZCluster(0.f);
                LArClusterHelper::GetClusterSpanX(pCluster, smallestXCluster, largestXCluster);
                LArClusterHelper::GetClusterSpanZ(pCluster, smallestXCluster, largestXCluster, smallestZCluster, largestZCluster);
                smallestZinV = std::min(smallestZinV, smallestZCluster);
                largestZinV = std::max(largestZinV, largestZCluster);
            }
            
            float smallestZinW(10000.f), largestZinW(0.f);
            for (const Cluster *const pCluster : clusterListW)
            {
                float smallestXCluster(0.f), largestXCluster(0.f), smallestZCluster(0.f), largestZCluster(0.f);
                LArClusterHelper::GetClusterSpanX(pCluster, smallestXCluster, largestXCluster);
                LArClusterHelper::GetClusterSpanZ(pCluster, smallestXCluster, largestXCluster, smallestZCluster, largestZCluster);
                smallestZinW = std::min(smallestZinW, smallestZCluster);
                largestZinW = std::max(largestZinW, largestZCluster);
            }
            

            for (VertexScoringTool::VertexScore &vertexScore : top5VertexScoreList)
            {
                const Vertex *const pVertex(vertexScore.GetVertex());
                CartesianVector vertexPosition(pVertex->GetPosition());
                
                const CartesianVector vertexProjectionU(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_U));
                const CartesianVector vertexProjectionV(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_V));
                const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition, TPC_VIEW_W));
                
                bool mcParticleUForwards(true), mcParticleVForwards(true), mcParticleWForwards(true);
                
                if (mcBeginPointU.GetZ() > mcEndPointU.GetZ())
                    mcParticleUForwards = false;
                if (mcBeginPointV.GetZ() > mcEndPointV.GetZ())
                    mcParticleUForwards = false;
                if (mcBeginPointW.GetZ() > mcEndPointW.GetZ())
                    mcParticleUForwards = false;
                
                int nBadViews(0);
                float USpan(std::abs(largestZinU - smallestZinU)), VSpan(std::abs(largestZinV - smallestZinV)), WSpan(std::abs(largestZinW - smallestZinW));
                
                int thresholdViews(0);
                
                if (USpan > 30.0)
                    thresholdViews++;
                if (VSpan > 30.0)
                    thresholdViews++;
                if (WSpan > 30.0)
                    thresholdViews++;
                
                if (((mcParticleUForwards == true && ((vertexProjectionU.GetZ() > smallestZinU + 0.75 * USpan)))
                || (mcParticleUForwards == false && ((vertexProjectionU.GetZ() < smallestZinU + 0.25 * USpan)))))
                    nBadViews++;
                    
                if ((mcParticleVForwards == true && ((vertexProjectionV.GetZ() > smallestZinV + 0.75 * VSpan)))
                || (mcParticleVForwards == false && ((vertexProjectionV.GetZ() < smallestZinV + 0.25 * VSpan))))
                    nBadViews++;
                    
                if ((mcParticleWForwards == true && ((vertexProjectionW.GetZ() > smallestZinW + 0.75 * WSpan)))
                || (mcParticleWForwards == false && ((vertexProjectionW.GetZ() < smallestZinW + 0.25 * WSpan))))
                    nBadViews++;
                
                if (nBadViews >= thresholdViews)
                {
                    std::cout << ">>> FILTERED" << std::endl;
                    continue;
                }
                    
                filteredTop5VertexScoreList.push_back(vertexScore);
            }
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "Caught exception" << std::endl;
        filteredTop5VertexScoreList = top5VertexScoreList;
        return;
    }
    
    if (filteredTop5VertexScoreList.size() == 0)
        filteredTop5VertexScoreList = top5VertexScoreList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::GetClusters(ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW)
{
    const ClusterList *pInputClusterListU(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "ClustersU", pInputClusterListU));

    if (!pInputClusterListU || pInputClusterListU->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CandidateVertexCreationAlgorithm: unable to find cluster list " << "ClustersU" << std::endl;

        return;
    }

    for (const Cluster *const pCluster : *pInputClusterListU)
            clusterListU.insert(pCluster);

    const ClusterList *pInputClusterListV(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "ClustersV", pInputClusterListV));

    if (!pInputClusterListV || pInputClusterListV->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CandidateVertexCreationAlgorithm: unable to find cluster list " << "ClustersV" << std::endl;

        return;
    }

    for (const Cluster *const pCluster : *pInputClusterListU)
        clusterListV.insert(pCluster);

    const ClusterList *pInputClusterListW(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "ClustersW", pInputClusterListW));

    if (!pInputClusterListW || pInputClusterListW->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CandidateVertexCreationAlgorithm: unable to find cluster list " << "ClustersW" << std::endl;

        return;
    }

    for (const Cluster *const pCluster : *pInputClusterListU)
        clusterListW.insert(pCluster);
}

//---------------------------------------------------------------------------------------------------------------------------------


const Cluster* VertexSelectionAlgorithm::GetLongestCluster(ClusterList &clusterList)
{
    float largestLength(0.f);
    
    for (const Cluster *const pCluster : clusterList)
    {
        if (LArClusterHelper::GetLength(pCluster) > largestLength)
            largestLength = LArClusterHelper::GetLength(pCluster);
    }
    
    for (const Cluster *const pCluster : clusterList)
    {
        if (LArClusterHelper::GetLength(pCluster) == largestLength)
            return pCluster;
    }
    
    return *(clusterList.begin());
}

//----------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmTool *pAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
        "VertexClustering", pAlgorithmTool));

    m_pVertexClusteringTool = dynamic_cast<VertexClusteringTool *>(pAlgorithmTool);

    AlgorithmTool *pAnotherAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
        "VertexScoring", pAnotherAlgorithmTool));

    m_pVertexScoringTool = dynamic_cast<VertexScoringTool *>(pAnotherAlgorithmTool);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TopologyVertexListName", m_topologyVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "EnergyVertexListName", m_energyVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Top5VertexListName", m_top5VertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "AllOtherVertexListName", m_allOtherVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamMode", m_beamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NDecayLengthsInZSpan", m_nDecayLengthsInZSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Kappa", m_kappa));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectSingleVertex", m_selectSingleVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreSelections", m_maxTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateDisplacement", m_minCandidateDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateScoreFraction", m_minCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDetectorGaps", m_useDetectorGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GapTolerance", m_gapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IsEmptyViewAcceptable", m_isEmptyViewAcceptable));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableClustering", m_enableClustering));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DirectionFilter", m_directionFilter));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamWeightFilter", m_beamWeightFilter));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
