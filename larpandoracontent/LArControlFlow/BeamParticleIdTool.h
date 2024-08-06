/**
 *  @file   larpandoracontent/LArControlFlow/BeamParticleIdTool.h
 *
 *  @brief  Header file for the beam particle id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_BEAM_PARTICLE_ID_TOOL_H
#define LAR_BEAM_PARTICLE_ID_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/SliceIdBaseTool.h"

namespace lar_content
{

/**
 *  @brief  BeamParticleIdTool class
 */
class BeamParticleIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    BeamParticleIdTool();

    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &beamSliceHypotheses,
        const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    /**
     *  @brief  Plane class
     */
    class Plane
    {
    public:
        /**
         *  @brief  Constructor, using equation of plane: m_a*x + m_b*y + m_c*z + m_d = 0.
         *
         *  @param  normal a Cartesian vector that points in a direction that is normal to the plane
         *  @param  point a Cartesian vector that corresponds to any point on the plane
         */
        Plane(const pandora::CartesianVector &normal, const pandora::CartesianVector &point);

        /**
         *  @brief  Return the intersection between the plane and a line
         *
         *  @param  a0 point on the line
         *  @param  a vector pointing along the line
         */
        pandora::CartesianVector GetLineIntersection(const pandora::CartesianVector &a0, const pandora::CartesianVector &a) const;

    private:
        pandora::CartesianVector m_unitNormal; ///< Unit normal to plane
        pandora::CartesianVector m_point;      ///< A point on the plane
        float m_d;                             ///< Parameter defining a plane
    };

    pandora::StatusCode Initialize();

    /**
     *  @brief  Select a given fraction of a slice's calo hits that are closest to the beam spot
     *
     *  @param  inputCaloHitList all calo hits in slice
     *  @param  outputCaloHitList to receive the list of selected calo hits
     *  @param  closestHitToFaceDistance to receive the distance of closest hit to beam spot
     */
    void GetSelectedCaloHits(const pandora::CaloHitList &inputCaloHitList, pandora::CaloHitList &outputCaloHitList, float &closestHitToFaceDistance) const;

    /**
     *  @brief  Find the intercepts of a line with the protoDUNE detector
     *
     *  @param  a0 a point on the line in question
     *  @param  majorAxis the direction of the line in question
     *  @param  interceptOne to receive the first intersection between line and protoDUNE detector
     *  @param  interceptTwo to receive the second intersection between line and protoDUNE detector
     */
    void GetTPCIntercepts(const pandora::CartesianVector &a0, const pandora::CartesianVector &majorAxis,
        pandora::CartesianVector &interceptOne, pandora::CartesianVector &interceptTwo) const;

    /**
     *  @brief  Check if a given 3D spacepoint is inside the global TPC volume
     *
     *  @param  spacePoint
     */
    bool IsContained(const pandora::CartesianVector &spacePoint) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<Plane> PlaneVector;

    bool m_selectAllBeamParticles;            ///< First approach: select all beam particles, as opposed to selecting all cosmics
    bool m_selectOnlyFirstSliceBeamParticles; ///< First approach: select first slice beam particles, cosmics for all subsequent slices
    float m_tpcMinX;                          ///< Global TPC volume minimum x extent
    float m_tpcMaxX;                          ///< Global TPC volume maximum x extent
    float m_tpcMinY;                          ///< Global TPC volume minimum y extent
    float m_tpcMaxY;                          ///< Global TPC volume maximum y extent
    float m_tpcMinZ;                          ///< Global TPC volume minimum z extent
    float m_tpcMaxZ;                          ///< Global TPC volume maximum z extent
    pandora::CartesianVector m_beamTPCIntersection; ///< Intersection of beam and global TPC volume
    pandora::CartesianVector m_beamDirection;       ///< Beam direction
    PlaneVector m_tpcPlanes;                        ///< Vector of all planes making up global TPC volume

    float m_projectionIntersectionCut; ///< Projection intersection distance cut, used in beam event selection
    float m_closestDistanceCut;        ///< Closest distance (of hit to beam spot), used in beam event selection
    float m_angleToBeamCut;            ///< Angle between major axis and beam direction, used in beam event selection
    float m_selectedFraction;          ///< Fraction of hits to use in 3D cluster fits
    unsigned int m_nSelectedHits;      ///< Minimum number of hits to use in 3D cluster fits
};

} // namespace lar_content

#endif // #ifndef LAR_BEAM_PARTICLE_ID_TOOL_H
