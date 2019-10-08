/**
 *  @file   larpandoracontent/LArCheating/CheatingVertexSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the cheating vertex selection algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingVertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void CheatingVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &/*beamConstants*/,
    HitKDTree2D &/*kdTreeU*/, HitKDTree2D &/*kdTreeV*/, HitKDTree2D &/*kdTreeW*/, VertexScoreList &vertexScoreList) const
{
    const Vertex *pBestVertex(nullptr);
    float bestVertexDr(std::numeric_limits<float>::max());
    this->GetBestVertex(vertexVector, pBestVertex, bestVertexDr);
    vertexScoreList.emplace_back(pBestVertex, 1.0f);
}

} // namespace lar_content
