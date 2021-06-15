/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/ShowerPfoMopUpAlgorithm.h
 *
 *  @brief  Header file for the shower pfo mop up algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_PFO_MOP_UP_ALGORITHM_H
#define LAR_SHOWER_PFO_MOP_UP_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/VertexBasedPfoMopUpAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ShowerPfoMopUpAlgorithm::Algorithm class
 */
class ShowerPfoMopUpAlgorithm : public VertexBasedPfoMopUpAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerPfoMopUpAlgorithm();

private:
    bool IsVertexAssociated(const pandora::CartesianVector &vertex2D, const LArPointingCluster &pointingCluster) const;
    PfoAssociation GetPfoAssociation(const pandora::Pfo *const pVertexPfo, const pandora::Pfo *const pDaughterPfo,
        HitTypeToAssociationMap &hitTypeToAssociationMap) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxVertexLongitudinalDistance; ///< Vertex association check: max longitudinal distance cut
    float m_vertexAngularAllowance;        ///< Vertex association check: pointing angular allowance in degrees
};

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_PFO_MOP_UP_ALGORITHM_H
