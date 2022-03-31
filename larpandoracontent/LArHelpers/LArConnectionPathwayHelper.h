/**
 *  @file   larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h
 *
 *  @brief  Header file for the connection pathway helper class.
 *
 *  $Log: $
 */
#ifndef LAR_CONNECTION_PATHWAY_HELPER_H
#define LAR_CONNECTION_PATHWAY_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "Objects/CartesianVector.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

namespace pandora
{
    class CartesianVector;
    class Pandora;
} // namespace pandora

namespace lar_content
{

/**
 *  @brief  LArConnectionPathwayHelper class
 */
class LArConnectionPathwayHelper
{
public:

   static bool AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxXSeparation, const float maxSeparation);

   static bool AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxXSeparation, const float maxSeparation, float &metric);

   static bool AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle);

   static bool AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle, float &metric);

   static bool FindShowerVertexFromPosition(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &showerStart3D);

   static bool FindShowerVertexFromDirection(pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector nuVertexPosition, 
       const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &showerStart3D);

    static bool FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector nuVertexPosition,
        const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxSeparation, 
        pandora::CartesianVector &showerStart3D);

};

} // namespace lar_content

#endif // #ifndef LAR_CONNECTION_PATHWAY_HELPER_H
