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
    /**
     *  @brief Identify pfos for which an attempt at 3D reclustering should be made
     *
     *  @param pPfos input list of all pfos
     *  @param pfosToRecluster output list of pfos that should be reclustered
     */
    virtual pandora::StatusCode GetPfosToRecluster(const pandora::PfoList *pPfos, pandora::PfoList &pfosToRecluster) = 0;

    /**
     *  @brief Calculate a measure of the goodness of a clustering
     *
     *  @param clusters input list of clusters
     *  @param fom output value of the clustering goodness, the figure of merit
     */
    virtual pandora::StatusCode CalcClusteringFom(const pandora::ClusterList &clusters, float &fom) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_RECLUSTERING_FIGURE_OF_MERIT_BASE_TOOL_H
