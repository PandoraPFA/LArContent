/**
 *  @file   larpandoracontent/LArVertex/MvaSelectionBaseAlgorithm.cc
 *
 *  @brief  Implementation of the mva vertex selection algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArVertex/EnergyKickFeatureTool.h"
#include "larpandoracontent/LArVertex/GlobalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/RPhiFeatureTool.h"
#include "larpandoracontent/LArVertex/ShowerAsymmetryFeatureTool.h"

#include "larpandoracontent/LArVertex/MvaVertexSelectionAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include <random>

using namespace pandora;

namespace lar_content
{

template <typename T>
MvaVertexSelectionAlgorithm<T>::MvaVertexSelectionAlgorithm() :
    TrainedVertexSelectionAlgorithm(),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void MvaVertexSelectionAlgorithm<T>::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    ClusterList clustersU, clustersV, clustersW;
    this->GetClusterLists(m_inputClusterListNames, clustersU, clustersV, clustersW);

    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    this->CalculateClusterSlidingFits(clustersU, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListU);
    this->CalculateClusterSlidingFits(clustersV, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListV);
    this->CalculateClusterSlidingFits(clustersW, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListW);

    ShowerClusterList showerClusterListU, showerClusterListV, showerClusterListW;
    this->CalculateShowerClusterList(clustersU, showerClusterListU);
    this->CalculateShowerClusterList(clustersV, showerClusterListV);
    this->CalculateShowerClusterList(clustersW, showerClusterListW);

    // Create maps from hit types to objects for passing to feature tools.
    const ClusterListMap clusterListMap{{TPC_VIEW_U, clustersU}, {TPC_VIEW_V, clustersV}, {TPC_VIEW_W, clustersW}};

    const SlidingFitDataListMap slidingFitDataListMap{
        {TPC_VIEW_U, slidingFitDataListU}, {TPC_VIEW_V, slidingFitDataListV}, {TPC_VIEW_W, slidingFitDataListW}};

    const ShowerClusterListMap showerClusterListMap{{TPC_VIEW_U, showerClusterListU}, {TPC_VIEW_V, showerClusterListV}, {TPC_VIEW_W, showerClusterListW}};

    const KDTreeMap kdTreeMap{{TPC_VIEW_U, kdTreeU}, {TPC_VIEW_V, kdTreeV}, {TPC_VIEW_W, kdTreeW}};

    // Calculate the event feature list and the vertex feature map.
    EventFeatureInfo eventFeatureInfo(this->CalculateEventFeatures(clustersU, clustersV, clustersW, vertexVector));

    LArMvaHelper::MvaFeatureVector eventFeatureList;
    this->AddEventFeaturesToVector(eventFeatureInfo, eventFeatureList);

    VertexFeatureInfoMap vertexFeatureInfoMap;
    for (const Vertex *const pVertex : vertexVector)
    {
        this->PopulateVertexFeatureInfoMap(
            beamConstants, clusterListMap, slidingFitDataListMap, showerClusterListMap, kdTreeMap, pVertex, vertexFeatureInfoMap);
    }

    // Use a simple score to get the list of vertices representing good regions.
    VertexScoreList initialScoreList;
    for (const Vertex *const pVertex : vertexVector)
        PopulateInitialScoreList(vertexFeatureInfoMap, pVertex, initialScoreList);

    VertexVector bestRegionVertices;
    this->GetBestRegionVertices(initialScoreList, bestRegionVertices);

    if (m_trainingSetMode)
        this->ProduceTrainingSets(vertexVector, bestRegionVertices, vertexFeatureInfoMap, eventFeatureList, kdTreeMap);

    if ((!m_trainingSetMode || m_allowClassifyDuringTraining) && !bestRegionVertices.empty())
    {
        // Use mva to choose the region.
        const Vertex *const pBestRegionVertex(
            this->CompareVertices(bestRegionVertices, vertexFeatureInfoMap, eventFeatureList, kdTreeMap, m_mvaRegion, m_useRPhiFeatureForRegion));

        // Get all the vertices in the best region.
        VertexVector regionalVertices{pBestRegionVertex};
        for (const Vertex *const pVertex : vertexVector)
        {
            if (pVertex == pBestRegionVertex)
                continue;

            if ((pBestRegionVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude() < m_regionRadius)
                regionalVertices.push_back(pVertex);
        }

        this->CalculateRPhiScores(regionalVertices, vertexFeatureInfoMap, kdTreeMap);

        if (!regionalVertices.empty())
        {
            // Use mva to choose the vertex and then fine-tune using the RPhi score.
            const Vertex *const pBestVertex(
                this->CompareVertices(regionalVertices, vertexFeatureInfoMap, eventFeatureList, kdTreeMap, m_mvaVertex, true));
            this->PopulateFinalVertexScoreList(vertexFeatureInfoMap, pBestVertex, vertexVector, vertexScoreList);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const pandora::Vertex *MvaVertexSelectionAlgorithm<T>::CompareVertices(const VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
    const LArMvaHelper::MvaFeatureVector &eventFeatureList, const KDTreeMap &kdTreeMap, const T &t, const bool useRPhi) const
{
    const Vertex *pBestVertex(vertexVector.front());
    LArMvaHelper::MvaFeatureVector chosenFeatureList;

    VertexFeatureInfo chosenVertexFeatureInfo(vertexFeatureInfoMap.at(pBestVertex));
    this->AddVertexFeaturesToVector(chosenVertexFeatureInfo, chosenFeatureList, useRPhi);

    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestVertex)
            continue;

        LArMvaHelper::MvaFeatureVector featureList;
        VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
        this->AddVertexFeaturesToVector(vertexFeatureInfo, featureList, useRPhi);

        if (!m_legacyVariables)
        {
            LArMvaHelper::MvaFeatureVector sharedFeatureList;
            float separation(0.f), axisHits(0.f);
            this->GetSharedFeatures(pVertex, pBestVertex, kdTreeMap, separation, axisHits);
            VertexSharedFeatureInfo sharedFeatureInfo(separation, axisHits);
            this->AddSharedFeaturesToVector(sharedFeatureInfo, sharedFeatureList);

            if (LArMvaHelper::Classify(t, LArMvaHelper::ConcatenateFeatureLists(eventFeatureList, featureList, chosenFeatureList, sharedFeatureList)))
            {
                pBestVertex = pVertex;
                chosenFeatureList = featureList;
            }
        }
        else
        {
	    if (LArMvaHelper::Classify(t, LArMvaHelper::ConcatenateFeatureLists(eventFeatureList, featureList, chosenFeatureList)))
            {
                pBestVertex = pVertex;
                chosenFeatureList = featureList;
            }
        }
    }

    return pBestVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode MvaVertexSelectionAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaFileName", m_mvaFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RegionMvaName", m_regionMvaName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexMvaName", m_vertexMvaName));

    // ATTN : Need access to base class member variables at this point, so call read settings prior to end of this function
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, TrainedVertexSelectionAlgorithm::ReadSettings(xmlHandle));

    if ((!m_trainingSetMode || m_allowClassifyDuringTraining))
    {
        if (m_mvaFileName.empty() || m_regionMvaName.empty() || m_vertexMvaName.empty())
        {
            std::cout << "MvaVertexSelectionAlgorithm: MvaFileName, RegionMvaName and VertexMvaName must be set if training set mode is"
                      << "off or we allow classification during training" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        const std::string fullMvaFileName(LArFileHelper::FindFileInPath(m_mvaFileName, m_filePathEnvironmentVariable));
        m_mvaRegion.Initialize(fullMvaFileName, m_regionMvaName);
        m_mvaVertex.Initialize(fullMvaFileName, m_vertexMvaName);
    }

    return STATUS_CODE_SUCCESS;
}

template class MvaVertexSelectionAlgorithm<AdaBoostDecisionTree>;
template class MvaVertexSelectionAlgorithm<SupportVectorMachine>;

} // namespace lar_content
