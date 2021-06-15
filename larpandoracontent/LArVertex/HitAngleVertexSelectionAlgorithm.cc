/**
 *  @file   larpandoracontent/LArVertex/HitAngleVertexSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the hit angle vertex selection algorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArVertex/RPhiFeatureTool.h"

#include "larpandoracontent/LArVertex/HitAngleVertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

HitAngleVertexSelectionAlgorithm::HitAngleVertexSelectionAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitAngleVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    const KDTreeMap kdTreeMap{{TPC_VIEW_U, kdTreeU}, {TPC_VIEW_V, kdTreeV}, {TPC_VIEW_W, kdTreeW}};

    float bestFastScore(0.f);
    for (const Vertex *const pVertex : vertexVector)
    {
        const float beamDeweightingScore(this->IsBeamModeOn() ? std::exp(this->GetBeamDeweightingScore(beamConstants, pVertex)) : 1.f);

        const float rPhiScore(LArMvaHelper::CalculateFeaturesOfType<RPhiFeatureTool>(m_featureToolVector, this, pVertex,
            SlidingFitDataListMap(), ClusterListMap(), kdTreeMap, ShowerClusterListMap(), beamDeweightingScore, bestFastScore)
                                  .at(0)
                                  .Get());

        vertexScoreList.emplace_back(pVertex, beamDeweightingScore * rPhiScore);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitAngleVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
