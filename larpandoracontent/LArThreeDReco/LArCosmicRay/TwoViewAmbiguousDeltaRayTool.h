/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewAmbiguousDeltaRayTool.h
 *
 *  @brief  Header file for the two view delta ray merge tool class
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_AMBIGUOUS_DELTA_RAY_TOOL_H
#define TWO_VIEW_AMBIGUOUS_DELTA_RAY_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoViewAmbiguousDeltaRayTool class
 */
class TwoViewAmbiguousDeltaRayTool : public DeltaRayMatrixTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    TwoViewAmbiguousDeltaRayTool();

private:
    bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Identify ambiguous matches (e.g. 3:2:1) and, if possible, create pfos out of the best 1:1:1 cluster match
     *
     *  @param  overlapMatrix the overlap matrix
     */
    void ExamineConnectedElements(MatrixType &overlapMatrix) const;

    /**
     *  @brief  Identify the best 1:1:1 match in a group of connected elements and from it create a pfo
     *
     *  @param  elementList the tensor element list
     *
     *  @return  whether any particles were created
     */
    bool PickOutGoodMatches(const MatrixType::ElementList &elementList) const;
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_AMBIGUOUS_DELTA_RAY_TOOL_H
