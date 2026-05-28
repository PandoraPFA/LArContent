/**
 *  @file   include/DLVertexingThreeDAlgorithm.h
 *
 *  @brief  Header file for the Deep Learning 3D Vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_VERTEXING_THREE_D_ALGORITHM_H
#define LAR_DL_VERTEXING_THREE_D_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningThreeDVertexingAlgorithm class
 */
class DlThreeDVertexingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief Default constructor
     */
    DlThreeDVertexingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode Infer();

    /**
     *  @brief Get the node data (positions and features) from the input CaloHit list.
     *
     *  @param caloHits The input calo hits.
     *  @param pos The positions of the nodes.
     *  @param node_features The features for each node.
     */
    pandora::StatusCode GetNodeData(
        const pandora::CaloHitList &caloHits, std::vector<pandora::CartesianVector> &pos, std::vector<std::array<float, 1>> &node_features);

    /**
     *  @brief Process the given node data into the expected tensor format for the model, and insert into the input vector.
     *
     *  @param inputs The input vector to be filled with the node data.
     *  @param pos The positions of the nodes.
     *  @param node_features The features for each node.
     */
    pandora::StatusCode BuildInput(
        LArDLHelper::TorchInputVector &inputs, std::vector<pandora::CartesianVector> &pos, std::vector<std::array<float, 1>> &node_features);

    LArDLHelper::TorchModel m_modelFile; ///< The model to use.

    float m_scalingFactor;               ///< The scaling factor for the input.
    std::vector<float> m_thresholds;     ///< Distance Class Thresholds.
    int m_nDistanceClasses;              ///< The number of distance classes (derived from thresholds).
    std::string m_caloHitListName;       ///< The name of the input CaloHit list.
    std::string m_outputVertexListName;  ///< The name of the output Vertex list to create with the predicted vertices.
};

} // namespace lar_dl_content

#endif // LAR_DL_THREE_D_VERTEXING_ALGORITHM_H
