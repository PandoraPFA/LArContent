/**
 *  @file   larpandoracontent/LArCustomParticles/TrackParticleBuildingAlgorithm.h
 *
 *  @brief  Header file for the 3D track building algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_PARTICLE_BUILDING_ALGORITHM_H
#define LAR_TRACK_PARTICLE_BUILDING_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArTrackPfo.h"

#include "larpandoracontent/LArCustomParticles/CustomParticleCreationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TrackParticleBuildingAlgorithm class
 */
class TrackParticleBuildingAlgorithm : public CustomParticleCreationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackParticleBuildingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void CreatePfo(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ParticleFlowObject*& pOutputPfo) const;

    /**
     *  @brief  Apply 3D sliding fit to Pfo and return track trajectory
     *
     *  @param  pPfo the address of the input Pfo
     *  @param  trackTrajectory the output track trajectory
     */
    void GetSlidingFitTrajectory(const pandora::ParticleFlowObject *const pPfo, const pandora::Vertex *const pVertex,
        LArTrackStateVector &trackStateVector) const;

    /**
     *  @brief  Sort pfos by number of constituent hits
     *
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortByHitProjection(const LArTrackTrajectoryPoint &lhs, const LArTrackTrajectoryPoint &rhs);

    bool            m_cosmicMode;             ///<
    unsigned int    m_slidingFitHalfWindow;   ///<
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_PARTICLE_BUILDING_ALGORITHM_H
