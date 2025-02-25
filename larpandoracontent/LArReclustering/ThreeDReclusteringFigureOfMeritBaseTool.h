/**
 *  @file   larpandoracontent/LArReclustering/ThreeDReclusteringFigureOfMeritBaseTool.h
 *
 *  @brief  Header file for the reclustering figure of merit algorithm tool base class.
 *
 *  $Log: $
 */

#ifndef LAR_THREE_D_RECLUSTERING_FIGURE_OF_MERIT_BASE_TOOL_H
#define LAR_THREE_D_RECLUSTERING_FIGURE_OF_MERIT_BASE_TOOL_H 1

#include "Pandora/AlgorithmTool.h"
#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
  *  @brief  ThreeDReclusteringFigureOfMeritBaseTool class
  */
class ThreeDReclusteringFigureOfMeritBaseTool : public pandora::AlgorithmTool
{
public:
    virtual pandora::StatusCode GetPfosToRecluster(const pandora::PfoList *pPfos, pandora::PfoList &pfosToRecluster) = 0;

    virtual pandora::StatusCode CalcClusteringFom(const pandora::ClusterList clusters, float &fom) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_RECLUSTERING_FIGURE_OF_MERIT_BASE_TOOL_H
