/**
 *  @file   larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the energy kick vertex selection algorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArVertex/EnergyKickFeatureTool.h"
#include "larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h"

#include "larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EnergyKickVertexSelectionAlgorithm::EnergyKickVertexSelectionAlgorithm() :
    m_minClusterCaloHits(12),
    m_slidingFitWindow(100),
    m_epsilon(0.06),
    m_asymmetryConstant(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D & /*kdTreeU*/, HitKDTree2D & /*kdTreeV*/, HitKDTree2D & /*kdTreeW*/, VertexScoreList &vertexScoreList) const
{
    ClusterList clustersU, clustersV, clustersW;
    this->GetClusterLists(m_inputClusterListNames, clustersU, clustersV, clustersW);

    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    this->CalculateClusterSlidingFits(clustersU, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListU);
    this->CalculateClusterSlidingFits(clustersV, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListV);
    this->CalculateClusterSlidingFits(clustersW, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListW);

    // Create maps from hit types to objects for passing to feature tools
    const SlidingFitDataListMap slidingFitDataListMap{
        {TPC_VIEW_U, slidingFitDataListU}, {TPC_VIEW_V, slidingFitDataListV}, {TPC_VIEW_W, slidingFitDataListW}};

    float bestFastScore(0.f); // not actually used - artefact of toolizing RPhi score and still using performance trick
    for (const Vertex *const pVertex : vertexVector)
    {
        const float beamDeweightingScore(this->IsBeamModeOn() ? this->GetBeamDeweightingScore(beamConstants, pVertex) : 0.f);

        const float energyKick(LArMvaHelper::CalculateFeaturesOfType<EnergyKickFeatureTool>(m_featureToolVector, this, pVertex,
            slidingFitDataListMap, ClusterListMap(), KDTreeMap(), ShowerClusterListMap(), beamDeweightingScore, bestFastScore)
                                   .at(0)
                                   .Get());

        const float energyAsymmetry(LArMvaHelper::CalculateFeaturesOfType<LocalAsymmetryFeatureTool>(m_featureToolVector, this, pVertex,
            slidingFitDataListMap, ClusterListMap(), KDTreeMap(), ShowerClusterListMap(), beamDeweightingScore, bestFastScore)
                                        .at(0)
                                        .Get());

        const float energyKickScore(-energyKick / m_epsilon);
        const float energyAsymmetryScore(energyAsymmetry / m_asymmetryConstant);

        vertexScoreList.push_back(VertexScore(pVertex, beamDeweightingScore + energyKickScore + energyAsymmetryScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyKickVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Epsilon", m_epsilon));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AsymmetryConstant", m_asymmetryConstant));

    if ((m_epsilon < std::numeric_limits<float>::epsilon()) || (m_asymmetryConstant < std::numeric_limits<float>::epsilon()))
    {
        std::cout << "EnergyKickVertexSelection: Invalid parameter(s), Epsilon " << m_epsilon << ", AsymmetryConstant "
                  << m_asymmetryConstant << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
