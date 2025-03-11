/**
 *  @file   larpandoracontent/LArReclustering/RandomFigureOfMeritTool.h
 *
 *  @brief  Header file for the reclustering random figure of merit algorithm tool class.
 *
 *  $Log: $
 */

#ifndef LAR_RANDOM_FIGURE_OF_MERIT_TOOL_H
#define LAR_RANDOM_FIGURE_OF_MERIT_TOOL_H 1

#include "larpandoracontent/LArReclustering/ThreeDReclusteringFigureOfMeritBaseTool.h"

namespace lar_content
{

/**
  *  @brief  RandomFigureOfMeritTool class
  */
class RandomFigureOfMeritTool : public ThreeDReclusteringFigureOfMeritBaseTool
{
public:

    /**
     *  @brief  Default constructor
     */
    RandomFigureOfMeritTool();

    /**
    *  @brief  Default destructor
    */
    ~RandomFigureOfMeritTool() = default;

    pandora::StatusCode GetPfosToRecluster(const pandora::PfoList *pPfos, pandora::PfoList &pfosToRecluster);

    pandora::StatusCode CalcClusteringFom(const pandora::ClusterList &clusters, float &fom);

private:
    float GetRandomFom();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxFomToRecluster; ///< threshold figure of merit for reclustering a pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float RandomFigureOfMeritTool::GetRandomFom() { return static_cast<float>(rand()) / RAND_MAX; }

} // namespace lar_content

#endif // #ifndef LAR_RANDOM_FIGURE_OF_MERIT_TOOL_H
