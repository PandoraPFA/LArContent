/**
 *  @file   larpandoracontent/LArReclustering/RandomFigureOfMeritTool.cc
 *
 *  @brief  Implementation file for the reclustering random figure of merit algorithm tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArReclustering/LArExample/RandomFigureOfMeritTool.h"

using namespace pandora;

namespace lar_content
{

RandomFigureOfMeritTool::RandomFigureOfMeritTool() :
    m_maxFomToRecluster {0.5f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

RandomFigureOfMeritTool::~RandomFigureOfMeritTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RandomFigureOfMeritTool::GetPfosToRecluster(const PfoList *pPfos, PfoList &pfosToRecluster)
{
    if (!pfosToRecluster.empty())
        return STATUS_CODE_FAILURE;

    for (const Pfo *const pPfo : *pPfos)
    {
        if (pPfo->GetNClusters() == 0)
            continue;
        if (GetRandomFom() < m_maxFomToRecluster)
            pfosToRecluster.push_back(pPfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RandomFigureOfMeritTool::CalcClusteringFom([[maybe_unused]] const ClusterList clusters, float &fom)
{
    fom = GetRandomFom();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RandomFigureOfMeritTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxFomToRecluster", m_maxFomToRecluster));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
