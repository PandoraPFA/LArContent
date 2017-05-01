/**
 *  @file   larpandoracontent/LArVertex/HitAngleVertexSelectionAlgorithm.h
 *
 *  @brief  Header file for the hit angle vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_ANGLE_VERTEX_SELECTION_ALGORITHM_H
#define LAR_HIT_ANGLE_VERTEX_SELECTION_ALGORITHM_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  HitAngleVertexSelectionAlgorithm class
 */
class HitAngleVertexSelectionAlgorithm : public VertexSelectionBaseAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    HitAngleVertexSelectionAlgorithm();

private:
    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    VertexFeatureTool::FeatureToolVector m_featureToolVector; ///< The feature tool map
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *HitAngleVertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new HitAngleVertexSelectionAlgorithm();
}
} // namespace lar_content

#endif // #ifndef LAR_HIT_ANGLE_VERTEX_SELECTION_ALGORITHM_H
