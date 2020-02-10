/**
 *  @file   larpandoracontent/DLVertexFeatureTool.h
 *
 *  @brief  Header file for the DLVertex feature tool class.
 *
 *  $Log: $
 */
#ifndef LAR_DLVERTEX_FEATURE_TOOL_H
#define LAR_DLVERTEX_FEATURE_TOOL_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DLVertexFeatureTool class
 */
class DLVertexFeatureTool : public VertexSelectionBaseAlgorithm::VertexFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLVertexFeatureTool();

    /**
     *  @brief  Run the tool
     *
     *  @param  featureVector containing the DLVertex feature
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pVertex address of the vertex
     *
     */
    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm * const pAlgorithm, const pandora::Vertex * const pVertex,
        const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &/*slidingFitDataListMap*/,const VertexSelectionBaseAlgorithm::ClusterListMap &,
        const VertexSelectionBaseAlgorithm::KDTreeMap &, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_DLVertexConstant;          ///< The DLVertex constant
    std::string     m_inputVertexListName;       ///< The name of the input vertex list
    unsigned int    m_numClusterCaloHitsPar;     ///< The number of cluster calo hits parameter 
};

} // namespace lar_content

#endif // #ifndef LAR_DLVERTEX_FEATURE_TOOL_H
