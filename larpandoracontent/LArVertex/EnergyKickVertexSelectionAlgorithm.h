/**
 *  @file   larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.h
 *
 *  @brief  Header file for the energy kick vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H
#define LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EnergyKickVertexSelectionAlgorithm class
 */
class EnergyKickVertexSelectionAlgorithm : public VertexSelectionBaseAlgorithm
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
    EnergyKickVertexSelectionAlgorithm();

private:
    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    VertexFeatureTool::FeatureToolVector m_featureToolVector; ///< The feature tool map

    pandora::StringVector   m_inputClusterListNames;        ///< The list of cluster list names
    unsigned int            m_minClusterCaloHits;           ///< The min number of hits parameter in the energy score
    unsigned int            m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    float                   m_epsilon;                      ///< The epsilon parameter in the energy score
    float                   m_asymmetryConstant;            ///< The asymmetry constant parameter in the energy score
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EnergyKickVertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new EnergyKickVertexSelectionAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H