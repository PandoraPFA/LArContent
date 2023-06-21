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
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

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

    class SortByDistanceToPoint
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  referencePoint the point relative to which constituent hits are ordered
         */
        SortByDistanceToPoint(const pandora::CartesianVector referencePoint) : m_referencePoint(referencePoint)
        {
        }

        /**
         *  @brief  Sort constituent hits by their position relative to a referencePoint
         *
         *  @param  lhs first constituent hit
         *  @param  rhs second constituent hit
         *
         *  @return  whether lhs hit is closer to the referencePoint than the rhs hit
         */
        bool operator()(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs);
        bool operator()(const pandora::CaloHit *const lhs, const pandora::CaloHit *const rhs);

    private:
        const pandora::CartesianVector m_referencePoint; ///< The point relative to which constituent hits are ordered
    };


    static bool FindShowerStarts3D(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, 
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianVector &nuVertexPosition, const float maxSeparationFromHit, 
        const float maxProjectionSeparation, pandora::CartesianPointVector &showerStarts3D);

   static bool FindShowerStartFromPosition(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
       const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &showerStart3D);

    static bool FindShowerStartFromDirection(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &uShowerStart3D, 
        pandora::CartesianVector &vShowerStart3D, pandora::CartesianVector &wShowerStart3D);

   static bool FindShowerStartFromDirection(const pandora::Algorithm *const pAlgorithm, const pandora::CartesianVector nuVertexPosition, 
       const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &showerStart3D);

    static bool FindShowerStartFromXProjection(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxSeparation, 
        pandora::CartesianVector &showerStart3D);

    static bool FindClosestSpinePosition(const ProtoShower &protoShower, const pandora::CartesianVector &showerStart3D, 
        pandora::CartesianVector &foundShowerStart);

    static bool FindShowerStartFromXProjectionRelaxed(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
       const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, pandora::CartesianVector &showerStart3D);

    static void GetMinMiddleMax(const float value1, const float value2, const float value3, float &minValue, float &middleValue,
        float &maxValue);
};

} // namespace lar_content

#endif // #ifndef LAR_CONNECTION_PATHWAY_HELPER_H
