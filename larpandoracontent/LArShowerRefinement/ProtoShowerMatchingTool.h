/**
 *  @file   larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h
 *
 *  @brief  Header file for the ProtoShower matching tool class.
  *
 *  $Log: $
 */
#ifndef LAR_PROTO_SHOWER_MATCHING_TOOL_H
#define LAR_PROTO_SHOWER_MATCHING_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

namespace lar_content
{

class ProtoShowerMatchingTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ProtoShowerMatchingTool();

    pandora::StatusCode Run(const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV,
        const ProtoShowerVector &protoShowerVectorW, ProtoShowerMatchVector &protoShowerMatchVector);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Determine whether three 2D connection pathways form a consistent 3D connection pathway
     *
     *  @param  protoShowerU the U view ProtoShower
     *  @param  protoShowerV the V view ProtoShower
     *  @param  protoShowerW the W view ProtoShower
     *  @param  consistency the basis of the match
     *
     *  @return whether three 2D connection pathways form a consistent 3D connection pathway
     */
    bool ArePathwaysConsistent(
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, Consistency &consistency) const;

    /**
     *  @brief  Determine whether three 2D shower start positions correspond to the same 3D shower start position
     *
     *  @param  protoShowerU the U view ProtoShower
     *  @param  protoShowerV the V view ProtoShower
     *  @param  protoShowerW the W view ProtoShower
     *
     *  @return whether three 2D shower start positions correspond to the same 3D shower start position
     */
    bool AreShowerStartsConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW) const;

    /**
     *  @brief  Determine whether three 2D initial spine directions correspond to the same 3D initial spine direction
     *
     *  @param  protoShowerU the U view ProtoShower
     *  @param  protoShowerV the V view ProtoShower
     *  @param  protoShowerW the W view ProtoShower
     *
     *  @return whether three 2D initial spine directions correspond to the same 3D initial spine direction
     */
    bool AreDirectionsConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW) const;

    /**
     *  @brief  Determine whether three 2D initial spine directions correspond to the same 3D initial spine direction
     *
     *  @param  nuVertexU the U view projection of the neutrino vertex
     *  @param  nuVertexV the V view projection of the neutrino vertex
     *  @param  nuVertexW the W view projection of the neutrino vertex
     *  @param  directionU the U view initial spine direction
     *  @param  directionU the V view initial spine direction
     *  @param  directionU the W view initial spine direction
     *
     *  @return whether three 2D initial spine directions correspond to the same 3D initial spine direction
     */
    bool AreDirectionsConsistent(const pandora::CartesianVector &nuVertexU, const pandora::CartesianVector &nuVertexV, const pandora::CartesianVector &nuVertexW,
        const pandora::CartesianVector &directionU, const pandora::CartesianVector &directionV, const pandora::CartesianVector &directionW) const;

    unsigned int m_spineSlidingFitWindow; ///< The shower spine sliding fit window
    float m_maxXSeparation;               ///< The max. drift-coordinate separation between matched 2D shower start positions
    float m_maxSeparation;                ///< The max. average separation between true and projected 2D shower start positions for a match
    float m_xExtrapolation;               ///< Extrapolation distance in the x-direction
    float m_maxAngularDeviation;          ///< The max. opening angle between true and projected 2D initial directions for a match
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_PROTO_SHOWER_MATCHING_TOOL_H
