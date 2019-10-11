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
    if (pBestVertex)
        vertexScoreList.emplace_back(pBestVertex, 1.0f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // ATTN : Need access to base class member variables at this point, so call read settings prior to end of this function
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, TrainedVertexSelectionAlgorithm::ReadSettings(xmlHandle));

    if (m_mcParticleListName.empty())
    {
        std::cout << "CheatingVertexSelectionAlgorithm: MCParticleListName required for cheated vertex selection" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
