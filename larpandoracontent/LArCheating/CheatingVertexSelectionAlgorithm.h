/**
 *  @file   larpandoracontent/LArCheating/CheatingVertexSelectionAlgorithm.h
 *
 *  @brief  Header file for the cheating vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_VERTEX_SELECTION_ALGORITHM_H
#define LAR_CHEATING_VERTEX_SELECTION_ALGORITHM_H 1

#include "larpandoracontent/LArVertex/TrainedVertexSelectionAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingVertexSelectionAlgorithm class
 */
class CheatingVertexSelectionAlgorithm : public TrainedVertexSelectionAlgorithm
{
private:
    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_VERTEX_SELECTION_ALGORITHM_H
