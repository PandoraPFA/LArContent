/**
 *  @file   ExampleContent/src/ExampleAlgorithms/ExampleDirectionAlgorithm.cc
 * 
 *  @brief  Implementation of the access lists algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent.h"

#include "larpandoracontent/LArDirection/ExampleDirectionAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "Api/PandoraContentApi.h"

using namespace pandora;

namespace lar_content
{

StatusCode ExampleDirectionAlgorithm::Run()
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());
    std::cout << ">>>>>> Cluster fits." << std::endl;

    for (const Cluster *const pCluster : clusterVector)
    {    
        if (LArClusterHelper::GetClusterHitType(pCluster) == TPC_VIEW_W)
        {
            try
            {
                TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetClusterDirection(pCluster);
                fitResult.DrawFit();

                std::cout << "Probability: " << fitResult.GetProbability() << std::endl;
                std::cout << "Vertex position: (" << fitResult.GetBeginpoint().GetX() << ", " << fitResult.GetBeginpoint().GetY() << ", " << fitResult.GetBeginpoint().GetZ() << ")" << std::endl;
                std::cout << "Endpoint position: (" << fitResult.GetEndpoint().GetX() << ", " << fitResult.GetEndpoint().GetY() << ", " << fitResult.GetEndpoint().GetZ() << ")" << std::endl;
            }
            catch (...)
            {
                std::cout << "Skipping..." << std::endl;
            }
        }
    }

    pandora::PfoVector pfoVector(pPfoList->begin(), pPfoList->end());
    std::cout << ">>>>>> PFO fits." << std::endl;

    for (const pandora::ParticleFlowObject *const pPfo : pfoVector)
    {
        try
        {
            TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetPfoDirection(pPfo);
            
            std::cout << "Probability: " << fitResult.GetProbability() << std::endl;
            std::cout << "Vertex position: (" << fitResult.GetBeginpoint().GetX() << ", " << fitResult.GetBeginpoint().GetY() << ", " << fitResult.GetBeginpoint().GetZ() << ")" << std::endl;
            std::cout << "Endpoint position: (" << fitResult.GetEndpoint().GetX() << ", " << fitResult.GetEndpoint().GetY() << ", " << fitResult.GetEndpoint().GetZ() << ")" << std::endl;
        }
        catch (...)
        {
                std::cout << "Skipping..." << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ExampleDirectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PfoListName", m_pfoListName));

    //AlgorithmTool *pAlgorithmTool(nullptr);

    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "LArTrackDirectionTool", pAlgorithmTool));
 
    //if (!(this->m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool*>(pAlgorithmTool)))
    //   throw STATUS_CODE_FAILURE;
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "TrackDirection", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
      {
        TrackDirectionTool *const pTrackDirectionTool(dynamic_cast<TrackDirectionTool*>(pAlgorithmTool));
	if (!pTrackDirectionTool) return STATUS_CODE_INVALID_PARAMETER;
	m_trackDirectionToolVector.push_back(pTrackDirectionTool);
      }


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

