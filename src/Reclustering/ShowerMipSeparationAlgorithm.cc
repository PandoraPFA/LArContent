/**
 *  @file   ShowerMipSeparationAlgorithm.cc
 * 
 *  @brief  Implementation of the shower-mip separation algorithm class.
 * 
 *  $Log: $
 */

#include "Helpers/ParticleIdHelper.h"

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterHelper.h"
#include "ShowerMipSeparationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ShowerMipSeparationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (pCluster->GetNCaloHits() < 10)
            continue;

        const ClusterHelper::ClusterFitResult &inputClusterFit(pCluster->GetFitToAllHitsResult());

        if (!inputClusterFit.IsFitSuccessful())
            continue;

        //ClusterHelper::ClusterFitResult inputClusterFit;
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitLayerCentroids(pCluster, pCluster->GetInnerPseudoLayer(), pCluster->GetOuterPseudoLayer(), inputClusterFit));

        ClusterList inputClusterList;
        inputClusterList.insert(pCluster);
        std::string inputClusterListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, TrackList(), inputClusterList, inputClusterListName));

        const ClusterList *pReclusterList = NULL;
        std::string reclusterListName;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithmName, pReclusterList, reclusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_associationAlgorithmName));

        // Now examine resulting cluster list
        std::string chosenListName(inputClusterListName);

        if (pReclusterList->size() > 1)
        {
            //PandoraMonitoringApi::VisualizeClusters(pReclusterList, "ReclusterList", RED);
            //PandoraMonitoringApi::ViewEvent();

            ClusterVector subClusterVector(pReclusterList->begin(), pReclusterList->end());
            std::sort(subClusterVector.begin(), subClusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

            for (ClusterVector::const_iterator subIter = subClusterVector.begin(), subIterEnd = subClusterVector.end(); subIter != subIterEnd; ++subIter)
            {
                Cluster *pSubCluster = *subIter;

                if (pSubCluster->GetNCaloHits() < 10)
                    continue;

                const ClusterHelper::ClusterFitResult &subClusterFit(pSubCluster->GetFitToAllHitsResult());

                if (!subClusterFit.IsFitSuccessful())
                    continue;

                //std::cout << " cosTheta   " << subClusterFit.GetDirection().GetCosOpeningAngle(inputClusterFit.GetDirection()) << std::endl;
                //std::cout << " nOccLayers " << pSubCluster->GetOrderedCaloHitList().size() << std::endl;
                //std::cout << " layerSpan  " << pSubCluster->GetOuterPseudoLayer() - pSubCluster->GetInnerPseudoLayer() << std::endl;
                //std::cout << " rms        " << subClusterFit.GetRms() << std::endl;
                //std::cout << " chi2       " << subClusterFit.GetChi2() << std::endl;

                //ClusterList subClusterList;
                //subClusterList.insert(pSubCluster);
                //PandoraMonitoringApi::VisualizeClusters(&subClusterList, "SubCluster", BLUE);
                //PandoraMonitoringApi::ViewEvent();
            }
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, chosenListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerMipSeparationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "Clustering",
        m_clusteringAlgorithmName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "ClusterAssociation",
        m_associationAlgorithmName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
