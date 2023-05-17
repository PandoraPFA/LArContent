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
    ProtoShowerMatchingTool();

    pandora::StatusCode Run(pandora::Algorithm *const pAlgorithm, const ElectronProtoShowerVector &protoShowerVectorU, 
        const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW, 
        ProtoShowerMatchVector &protoShowerMatchVector);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool ArePathwaysConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, 
        const ProtoShower &protoShowerW, Consistency &consistency);

    float m_maxXSeparation;
    float m_maxSeparation;
    float m_maxAngularDeviation;
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_PROTO_SHOWER_MATCHING_TOOL_H
