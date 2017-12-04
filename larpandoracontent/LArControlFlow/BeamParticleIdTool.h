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

namespace lar_content
{

/**
 *  @brief  BeamParticleIdTool class
 */
class BeamParticleIdTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief Constructor
     */
    BeamParticleIdTool();

    /**
     *  @brief Initialize geometry related memeber variables 
     */
    pandora::StatusCode Initialize();

    /**
     *  @brief Select reconstructed beam Pfos from beam and cr slice hypotheses
     *
     *  @param beamSliceHypotheses beam slice hypotheses
     *  @param crSliceHypotheses cosmic ray slice hypotheses
     *  @param selectedPfos Pfos identified as beam particles 
     */
    void SelectOutputPfos(const SliceHypotheses &beamSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos);

private:
    class Plane
    {
    public:
        /**
         *  @brief Constructor for the plane class.  Equation of plane used here m_a*x + m_b*y + m_c*z + m_d = 0.
         *
         *  @param normal a Cartesian vector that points in a direction that is normal to the plane
         *  @param point a Cartesian vector that corresponds to any point on the plane
         */
        Plane(pandora::CartesianVector &normal, pandora::CartesianVector &point);

        /**
         *  @brief Return the intersection between the plane and a line 
         *
         *  @param a0 point on the line
         *  @param a vector pointing along the line
         */
        pandora::CartesianVector GetLineIntersection(pandora::CartesianVector &a0, pandora::CartesianVector &a) const;

    private:
        float                         m_a;                       ///< Parameter defining a plane
        float                         m_b;                       ///< Parameter defining a plane
        float                         m_c;                       ///< Parameter defining a plane
        float                         m_d;                       ///< Parameter defining a plane
        pandora::CartesianVector      m_unitNormal;              ///< Unit normal to plane
        pandora::CartesianVector      m_point;                   ///< A point on the plane
    };

    /**
     *  @brief Find the intercepts of a line with the protoDUNE detector
     *
     *  @param a0 a point on the line in question
     *  @param majorAxis the direction of the line in question
     *  @param interceptOne the first intersection between line and protoDUNE detector 
     *  @param interceptTwo the second intersection between line and protoDUNE detector 
     */
    pandora::StatusCode GetTPCIntercepts(pandora::CartesianVector &a0, pandora::CartesianVector &majorAxis, pandora::CartesianVector &interceptOne,  pandora::CartesianVector &interceptTwo) const;

    /**
     *  @brief Check if a given 3D spacepoint is inside the global TPC volume
     *
     *  @param spacePoint
     */
    bool IsContained(pandora::CartesianVector &spacePoint) const;

    /**
     *  @brief Read settings via xml
     *  
     *  @param xmlHanle
     */
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<const Plane*> PlaneVector;
    bool                      m_selectAllBeamParticles;               ///< First approach: select all beam particles, as opposed to selecting all cosmics
    bool                      m_selectOnlyFirstSliceBeamParticles;    ///< First approach: select first slice beam particles, cosmics for all subsequent slices
    bool                      m_visualizeID;                          ///< Produce an event display showing the ID on a per slice basis
    float                     m_projectionIntersectionCut;            ///< Distance cut used in beam event selection
    float                     m_tpcMinX;                              ///< Global TPC volume minimum x extent 
    float                     m_tpcMaxX;                              ///< Global TPC volume maximum x extent
    float                     m_tpcMinY;                              ///< Global TPC volume minimum y extent
    float                     m_tpcMaxY;                              ///< Global TPC volume maximum y extent
    float                     m_tpcMinZ;                              ///< Global TPC volume minimum z extent
    float                     m_tpcMaxZ;                              ///< Global TPC volume maximum z extent
    pandora::CartesianVector  m_beamTPCIntersection;                  ///< Intersection of beam and global TPC volume
    PlaneVector               m_tpcPlanes;                            ///< Vector of all planes making up global TPC volume
};

} // namespace lar_content

#endif // #ifndef LAR_BEAM_PARTICLE_ID_TOOL_H
