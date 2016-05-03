/**
 *  @file   LArContent/include/LArThreeDReco/LArPfoMopUp/SplitShowerMergingAlgorithm.h
 * 
 *  @brief  Header file for the split shower merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SPLIT_SHOWER_MERGING_ALGORITHM_H
#define LAR_SPLIT_SHOWER_MERGING_ALGORITHM_H 1

#include "LArThreeDReco/LArPfoMopUp/VertexBasedPfoMergingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  SplitShowerMergingAlgorithm::Algorithm class
 */
class SplitShowerMergingAlgorithm : public VertexBasedPfoMergingAlgorithm
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
    SplitShowerMergingAlgorithm();

private:
    bool IsVertexAssociated(const pandora::CartesianVector &vertex2D, const LArPointingCluster &pointingCluster) const;
    PfoAssociation GetPfoAssociation(const pandora::Pfo *const pVertexPfo, const pandora::Pfo *const pDaughterPfo, HitTypeToAssociationMap &hitTypeToAssociationMap) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float                   m_maxVertexLongitudinalDistance;    ///< Vertex association check: max longitudinal distance cut
    float                   m_vertexAngularAllowance;           ///< Vertex association check: pointing angular allowance in degrees
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SplitShowerMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new SplitShowerMergingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_SPLIT_SHOWER_MERGING_ALGORITHM_H
