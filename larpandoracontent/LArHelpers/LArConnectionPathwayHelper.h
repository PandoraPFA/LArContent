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

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

namespace lar_content
{

/**
 *  @brief  LArConnectionPathwayHelper class
 */
class LArConnectionPathwayHelper
{
public:
    /**
     *  @brief  SortByDistanceToPoint class
     */
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

    /**
     *  @brief  Create 3D shower start position(s) from three input 2D positions
     *
     *  @param  pAlgorithm the calling algorithm
     *  @param  pShowerPfo the shower pfo
     *  @param  protoShowerMatch the ProtoShowerMatch object
     *  @param  nuVertexPosition the 3D neutrino vertex
     *  @param  maxSeparationFromHit the max. separation between a projected 3D shower start and the closest 2D shower hit
     *  @param  maxProjectionSeparation the max. separation between the projected 3D shower start and the shower start of that view
     *  @param  maxXSeparation the max. drift-coordinate separation between a 3D shower start and a matched 2D shower hit
     *  @param  showerStarts3D the output vector of 3D shower starts (ordered closest -> furthest from the neutrino vertex)
     *
     *  @return whether a consistent 3D shower start position could be created
     */
    static bool FindShowerStarts3D(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo,
        const ProtoShowerMatch &protoShowerMatch, const pandora::CartesianVector &nuVertexPosition, const float maxSeparationFromHit,
        const float maxProjectionSeparation, const float maxXSeparation, pandora::CartesianPointVector &showerStarts3D);

    /**
     *  @brief  Create 3D shower start position from three input 2D positions, assuming consistency of position
     *
     *  @param  pAlgorithm the calling algorithm
     *  @param  protoShowerU the U view ProtoShower
     *  @param  protoShowerV the V view ProtoShower
     *  @param  protoShowerW the W view ProtoShower
     *  @param  showerStart3D the output 3D shower start position
     */
    static void FindShowerStartFromPosition(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &showerStart3D);

    /**
     *  @brief  Create a 3D shower start position from each input 2D position, assuming consistency of initial direction
     *
     *  @param  pAlgorithm the calling algorithm
     *  @param  protoShowerU the U view ProtoShower
     *  @param  protoShowerV the V view ProtoShower
     *  @param  protoShowerW the W view ProtoShower
     *  @param  uShowerStart3D the output 3D shower start created from the U view shower start position
     *  @param  vShowerStart3D the output 3D shower start created from the V view shower start position
     *  @param  wShowerStart3D the output 3D shower start created from the W view shower start position
     *
     *  @return whether the 3D shower start positions could be created
     */
    static bool FindShowerStartFromDirection(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
        const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, pandora::CartesianVector &uShowerStart3D,
        pandora::CartesianVector &vShowerStart3D, pandora::CartesianVector &wShowerStart3D);

    /**
     *  @brief  Create a 3D shower start position from an input 2D position, assuming consistency of the drift coordinate
     *
     *  @param  pAlgorithm the calling algorithm
     *  @param  protoShower the ProtoShower from which the input 2D shower start is obtained
     *  @param  protoShower1 a matched ProtoShower from another view
     *  @param  protoShower2 the matched ProtoShower from the final view
     *  @param  maxSeparation the max. separation between the projected 3D shower start and the shower start of that view     
     *  @param  maxXSeparation the max. drift-coordinate separation between a 3D shower start and a matched 2D shower hit
     *  @param  showerStart3D the output 3D shower start
     *
     *  @return whether the 3D shower start positions could be created
     */
    static bool FindShowerStartFromXProjection(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, const ProtoShower &protoShower1,
        const ProtoShower &protoShower2, const float maxSeparation, const float maxXSeparation, pandora::CartesianVector &showerStart3D);

    /**
     *  @brief  Find the 2D spine hit that is closest to the neutrino vertex and shares a drift coordinate with an input 3D position
     *
     *  @param  protoShower the input ProtoShower in which to search
     *  @param  showerStart3D the input 3D shower start position
     *  @param  maxXSeparation the max. drift-coordinate separation between a 3D shower start and a matched 2D shower hit 
     *  @param  foundShowerStart the output matched 2D hit position
     *
     *  @return whether a 2D spine hit could be found
     */
    static bool FindClosestSpinePosition(const ProtoShower &protoShower, const pandora::CartesianVector &showerStart3D,
        const float maxXSeparation, pandora::CartesianVector &foundShowerStart);

    /**
     *  @brief  A relaxed approach to create a 3D shower start position from an input 2D position, assuming consistency of the drift coordinate
     *
     *  @param  pAlgorithm the calling algorithm
     *  @param  protoShower the ProtoShower from which the input 2D shower start is obtained
     *  @param  protoShower1 a matched ProtoShower from another view
     *  @param  protoShower2 the matched ProtoShower from the final view
     *  @param  maxSeparation the max. separation between the projected 3D shower start and the shower start of that view     
     *  @param  maxXSeparation the max. drift-coordinate separation between a 3D shower start and a matched 2D shower hit
     *  @param  ShowerStart3D the output 3D shower start
     *
     *  @return whether the 3D shower start position could be created
     */
    static bool FindShowerStartFromXProjectionRelaxed(const pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower,
        const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, const float maxXSeparation,
        pandora::CartesianVector &showerStart3D);

    /**
     *  @brief  Determine the lowest, median and highest value from an input of three numbers
     *
     *  @param  value1 the first value
     *  @param  value2 the second value
     *  @param  value3 the third value
     *  @param  minValue the minimum value
     *  @param  middleValue the median value
     *  @param  maxValue the maximum value
     */
    static void GetMinMiddleMax(const float value1, const float value2, const float value3, float &minValue, float &middleValue, float &maxValue);
};

} // namespace lar_content

#endif // #ifndef LAR_CONNECTION_PATHWAY_HELPER_H
