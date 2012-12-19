/**
 *  @file   VertexSeedFindingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex seed finding algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterHelper.h"
#include "LArPointingCluster.h"
#include "LArVertexHelper.h"
#include "VertexSeedFindingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode VertexSeedFindingAlgorithm::Run()
{
    // Store details of cluster inner and outer vertices
    const ClusterList *pInputClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pInputClusterList));

    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pInputClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    // Create list of clusters associated with the vertex
    ClusterList vertexSeedClusterList;

    const CartesianVector eventVertex(LArVertexHelper::GetCurrentVertex());

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        LArPointingCluster pointingCluster(*iter);

        if (LArClusterHelper::IsNode(eventVertex, pointingCluster.GetInnerVertex().GetPosition()) ||
            LArClusterHelper::IsNode(eventVertex, pointingCluster.GetOuterVertex().GetPosition()) ||
            LArClusterHelper::IsEmitted(eventVertex, pointingCluster.GetInnerVertex()) ||
            LArClusterHelper::IsEmitted(eventVertex, pointingCluster.GetOuterVertex()))
        {
            vertexSeedClusterList.insert(*iter);
        }
    }

    this->MakeVertexSeedMerges(eventVertex, vertexSeedClusterList);

    // Cluster list management
    ClusterList nonSeedClusterList(*pInputClusterList);

    for (ClusterList::const_iterator iter = vertexSeedClusterList.begin(), iterEnd = vertexSeedClusterList.end(); iter != iterEnd; ++iter)
    {
        nonSeedClusterList.erase(*iter);
    }

    if (!vertexSeedClusterList.empty())
    {
        if (!nonSeedClusterList.empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListName, nonSeedClusterList));
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName, vertexSeedClusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentClusterList(*this, m_seedClusterListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSeedFindingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        const unsigned int innerPseudoLayer(pCluster->GetInnerPseudoLayer());
        const unsigned int outerPseudoLayer(pCluster->GetOuterPseudoLayer());

        if (outerPseudoLayer - innerPseudoLayer < m_minClusterLayers)
            continue;

        const CartesianVector innerCentroid(pCluster->GetCentroid(innerPseudoLayer));
        const CartesianVector outerCentroid(pCluster->GetCentroid(outerPseudoLayer));

        if ((outerCentroid - innerCentroid).GetMagnitudeSquared() < m_minClusterLengthSquared)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSeedFindingAlgorithm::MakeVertexSeedMerges(const CartesianVector &eventVertex, ClusterList &vertexSeedClusterList) const
{
    for (ClusterList::iterator iterI = vertexSeedClusterList.begin(); iterI != vertexSeedClusterList.end(); ++iterI)
    {
        Cluster *pClusterI = *iterI;

        for (ClusterList::iterator iterJ = vertexSeedClusterList.begin(); iterJ != vertexSeedClusterList.end(); )
        {
            Cluster *pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
            {
                ++iterJ;
                continue;
            }

            // TODO, improve efficiency - currently need to recalculate every time
            LArPointingCluster pointingClusterI(*iterI);
            const bool isNodeInnerI(LArClusterHelper::IsNode(eventVertex, pointingClusterI.GetInnerVertex().GetPosition()));
            const bool isNodeOuterI(LArClusterHelper::IsNode(eventVertex, pointingClusterI.GetOuterVertex().GetPosition()));

            LArPointingCluster pointingClusterJ(*iterJ);
            const bool isNodeInnerJ(LArClusterHelper::IsNode(eventVertex, pointingClusterJ.GetInnerVertex().GetPosition()));
            const bool isNodeOuterJ(LArClusterHelper::IsNode(eventVertex, pointingClusterJ.GetOuterVertex().GetPosition()));
            const bool isEmittedInnerJ(LArClusterHelper::IsEmitted(eventVertex, pointingClusterJ.GetInnerVertex()));
            const bool isEmittedOuterJ(LArClusterHelper::IsEmitted(eventVertex, pointingClusterJ.GetOuterVertex()));

            // Whether to remove cluster J from list of vertex seeds
            bool removeDaughter(false);

            if (isNodeInnerI && !isNodeInnerJ && isEmittedInnerJ &&
                (2 * pClusterJ->GetOrderedCaloHitList().size() < pClusterI->GetOrderedCaloHitList().size()) &&
                ((pointingClusterJ.GetInnerVertex().GetPosition() - pointingClusterI.GetInnerVertex().GetPosition()).GetMagnitudeSquared() > 5.f * 5.f) &&
                (LArClusterHelper::GetClosestDistance(pointingClusterJ.GetInnerVertex().GetPosition(), pClusterI) < 2.f))
            {
                removeDaughter = true;
            }
            else if (isNodeOuterI && !isNodeOuterJ && isEmittedOuterJ &&
                (2 * pClusterJ->GetOrderedCaloHitList().size() < pClusterI->GetOrderedCaloHitList().size()) &&
                ((pointingClusterJ.GetOuterVertex().GetPosition() - pointingClusterI.GetOuterVertex().GetPosition()).GetMagnitudeSquared() > 5.f * 5.f) &&
                (LArClusterHelper::GetClosestDistance(pointingClusterJ.GetOuterVertex().GetPosition(), pClusterI) < 2.f))
            {
                removeDaughter = true;
            }

            // Whether to merge clusters I and J
            bool mergeDaughter(false);

            if (!isNodeOuterI && !isNodeInnerJ && isEmittedInnerJ &&
                ((pointingClusterI.GetOuterVertex().GetPosition() - pointingClusterJ.GetInnerVertex().GetPosition()).GetMagnitudeSquared() < 2.f * 2.f))
            {
                mergeDaughter = true;
            }
            else if (!isNodeInnerI && !isNodeOuterJ && isEmittedOuterJ &&
                ((pointingClusterI.GetInnerVertex().GetPosition() - pointingClusterJ.GetOuterVertex().GetPosition()).GetMagnitudeSquared() < 2.f * 2.f))
            {
                mergeDaughter = true;
            }

//if (removeDaughter || mergeDaughter)
//{
//    if (removeDaughter) std::cout << "Remove daughter " << std::endl;
//    if (mergeDaughter) std::cout << "Merge daughter " << std::endl;
//    if (removeDaughter && mergeDaughter) std::cout << "Both merge and remove daughter - PROBLEM! " << std::endl;
//    ClusterList parent, daughter; parent.insert(pClusterI); daughter.insert(pClusterJ);
//    PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//    PandoraMonitoringApi::VisualizeClusters(&vertexSeedClusterList, "vertexSeedClusterList", BLUE);
//    PandoraMonitoringApi::VisualizeClusters(&parent, "parent", RED);
//    PandoraMonitoringApi::VisualizeClusters(&daughter, "daughter", GREEN);
//    PandoraMonitoringApi::ViewEvent();
//}
            if (removeDaughter)
            {
                vertexSeedClusterList.erase(iterJ);
                iterJ = vertexSeedClusterList.begin();
            }
            else if (mergeDaughter)
            {
                vertexSeedClusterList.erase(iterJ);
                iterJ = vertexSeedClusterList.begin();
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClusterI, pClusterJ));
            }
            else
            {
                ++iterJ;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSeedFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    float minClusterLength = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    m_minClusterLayers = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterLayers", m_minClusterLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
