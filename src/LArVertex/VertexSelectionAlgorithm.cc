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
    m_enableEnergyCandidates(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::Run()
{
    const VertexList *pTopologyVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_topologyVertexListName, pTopologyVertexList));
    
    const VertexList *pEnergyVertexList(NULL);
    if (m_enableEnergyCandidates)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_energyVertexListName, pEnergyVertexList));

    const VertexList inputVertexListCopy(*pTopologyVertexList);

    if (!pTopologyVertexList || pTopologyVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexSelectionAlgorithm: unable to find current vertex list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    if (!pEnergyVertexList || pEnergyVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexSelectionAlgorithm: there is no energy vertex list. Energy vertices are probably turned off." << std::endl;
    }


    std::vector<const VertexList*> vertexListVector = m_pVertexClusteringTool->ClusterVertices(pTopologyVertexList);
    
    std::vector<VertexScoringTool::VertexScoreList> scoredClusterCollection;
    m_pVertexScoringTool->ScoreVertices(this, pTopologyVertexList, vertexListVector, scoredClusterCollection);

    VertexScoringTool::VertexScoreList energyVertexScoreList;

    if (m_enableEnergyCandidates && !pEnergyVertexList->empty())
    {
        m_pVertexScoringTool->ScoreEnergyVertices(this, pEnergyVertexList, energyVertexScoreList);
        std::sort(energyVertexScoreList.begin(), energyVertexScoreList.end());
    }

    VertexList selectedVertexList;
    this->SelectTopScoreVertices(scoredClusterCollection, selectedVertexList);
    
    if (!selectedVertexList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, selectedVertexList));
    
    //---------------------------------------------------------------------------------------------------------------------------------------
    
    VertexScoringTool::VertexScoreList intermediateVertexScoreList;

    unsigned int clusterCounter(0);
    
    if (m_enableClustering)
    {
        if (m_enableEnergyCandidates && !energyVertexScoreList.empty())
        {
            intermediateVertexScoreList.push_back(*energyVertexScoreList.begin());
            clusterCounter++;
        }

        for (VertexScoringTool::VertexScoreList &tempVertexScoreList : scoredClusterCollection)
        {
            if (clusterCounter == 5 || clusterCounter == scoredClusterCollection.size())
                break;
                
            if (tempVertexScoreList.size() == 0)
                continue;
            
            std::sort(tempVertexScoreList.begin(), tempVertexScoreList.end());
            intermediateVertexScoreList.push_back(*tempVertexScoreList.begin());
            clusterCounter++;
        }
    }
    else
    {
        VertexScoringTool::VertexScoreList anotherVertexScoreList;
        
        for (VertexScoringTool::VertexScoreList &tempVertexScoreList : scoredClusterCollection)
        {
            if (tempVertexScoreList.size() == 0)
                continue;
            
            anotherVertexScoreList.insert(anotherVertexScoreList.begin(), tempVertexScoreList.begin(), tempVertexScoreList.end());
        }

        if (m_enableEnergyCandidates && !energyVertexScoreList.empty())
            anotherVertexScoreList.insert(anotherVertexScoreList.begin(), energyVertexScoreList.begin(), energyVertexScoreList.end());
        
        std::sort(anotherVertexScoreList.begin(), anotherVertexScoreList.end());
        
        for (VertexScoringTool::VertexScore &anotherTempVertexScore : anotherVertexScoreList)
        {
            if (clusterCounter == 5 || clusterCounter == anotherVertexScoreList.size())
                break;
            
            intermediateVertexScoreList.push_back(anotherTempVertexScore);
            clusterCounter++;
        }
    }
    
    std::sort(intermediateVertexScoreList.begin(), intermediateVertexScoreList.end());
    
    const VertexList *pTop5TemporaryList(NULL);
    std::string top5TemporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTop5TemporaryList, top5TemporaryListName));
    
    for (const VertexScoringTool::VertexScore &vertexScore: intermediateVertexScoreList)
    {
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = vertexScore.GetVertex()->GetPosition();
        parameters.m_vertexLabel = vertexScore.GetVertex()->GetVertexLabel();
        parameters.m_vertexType = vertexScore.GetVertex()->GetVertexType();

        const Vertex *pVertex(NULL);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
    }
    
    VertexList top5Vertices;
    for (const Vertex *const pVertex : (*pTop5TemporaryList))
        top5Vertices.insert(pVertex);
    
    if (!top5Vertices.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_top5VertexListName, top5Vertices));
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    const VertexList *pAllVerticesTemporaryList(NULL);
    std::string allVerticesTemporaryList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pAllVerticesTemporaryList, allVerticesTemporaryList));
    
    for (const Vertex *const pVertex : inputVertexListCopy)
    {
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = pVertex->GetPosition();
        parameters.m_vertexLabel = pVertex->GetVertexLabel();
        parameters.m_vertexType = pVertex->GetVertexType();

        const Vertex *pVertexClone(NULL);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertexClone));
    }
    
    VertexList allVertices;
    for (const Vertex *const pVertex : (*pAllVerticesTemporaryList))
        allVertices.insert(pVertex);
    
    if (!allVertices.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_allOtherVertexListName, allVertices));
    
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
    m_pVertexScoringTool->NormaliseVertexScores(vertexScoreList);
    
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
        "EnableEnergyCandidates", m_enableEnergyCandidates));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
