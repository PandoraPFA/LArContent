/**
 *  @file   LArContent/include/LArHelpers/LArPfoHelper.h
 *
 *  @brief  Header file for the pfo helper class.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_HELPER_H
#define LAR_PFO_HELPER_H 1

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Vertex.h"

namespace lar_content
{

/**
 *  @brief  LArPfoHelper class
 */
class LArPfoHelper
{
public:
    /**
     *  @brief  Get a list of calo hits of a particular hit type from a list of pfos
     *
     *  @param  pfoList the input list of Pfos
     *  @param  hitType the cluster hit type
     *  @param  caloHitList the output list of calo hits
     */
    static void GetCaloHits(const pandora::PfoList &pfoList, const pandora::HitType &hitType, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get a list of calo hits of a particular hit type from a given pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  hitType the cluster hit type
     *  @param  caloHitList the output list of calo hits
     */
    static void GetCaloHits(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType &hitType, pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Get a list of clusters of a particular hit type from a list of pfos
     *
     *  @param  pfoList the input list of Pfos
     *  @param  hitType the cluster hit type
     *  @param  clusterList the output list of clusters
     */
    static void GetClusters(const pandora::PfoList &pfoList, const pandora::HitType &hitType, pandora::ClusterList &clusterList);

    /**
     *  @brief  Get a list of clusters of a particular hit type from a given pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  hitType the cluster hit type
     *  @param  clusterList the output list of clusters
     */
    static void GetClusters(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType &hitType, pandora::ClusterList &clusterList);

    /**
     *  @brief  Get a flat list of all pfos, recursively including all daughters and parents associated with those pfos in an input list
     *
     *  @param  inputPfoList the input pfo list
     *  @param  outputPfoList to receive the output pfo list
     */
    static void GetAllConnectedPfos(const pandora::PfoList &inputPfoList, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Get a flat list of all pfos, recursively including all daughters and parents associated with an input pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  outputPfoList to receive the output pfo list
     */
    static void GetAllConnectedPfos(pandora::ParticleFlowObject *const pPfo, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Get a flat list of all pfos, recursively, of all daughters associated with those pfos in an input list
     *
     *  @param  inputPfoList the input pfo list
     *  @param  outputPfoList to receive the output pfo list
     */
    static void GetAllDownstreamPfos(const pandora::PfoList &inputPfoList, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Get a flat list of all pfos, recursively, of all daughters and parents associated with an input pfo
     *
     *  @param  pPfo the input Pfo
     *  @param  outputPfoList to receive the output pfo list
     */
    static void GetAllDownstreamPfos(pandora::ParticleFlowObject *const pPfo, pandora::PfoList &outputPfoList);

    /**
     *  @brief  Calculate length of Pfo using 2D clusters
     *
     *  @param  pPfo the input Pfo
     *
     *  @return  length variable
     */
    static float GetTwoDLengthSquared(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Calculate length of Pfo using 2D clusters
     *
     *  @param  pPfo the input Pfo
     *
     *  @return  length variable
     */
    static float GetThreeDLengthSquared(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get closest distance between Pfo and cluster
     *
     *  @param  pPfo the address of the input Pfo
     *  @param  pCluster the address of the input cluster
     */
    static float GetClosestDistance(const pandora::ParticleFlowObject *const pPfo, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get distance between two Pfos using 2D clusters
     *
     *  @param  pPfo the address of the first Pfo
     *  @param  pPfo the address of the second Pfo
     */
    static float GetTwoDSeparation(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2);

    /**
     *  @brief  Get distance between two Pfos using 3D clusters
     *
     *  @param  pPfo the address of the first Pfo
     *  @param  pPfo the address of the second Pfo
     */
    static float GetThreeDSeparation(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2);

    /**
     *  @brief  Apply 3D sliding fit to Pfo and return vector of track states
     *
     *  @param  pPfo the address of the first Pfo
     *  @param  layerWindow the half-window used in the sliding fit
     *  @param  layerPitch the pitch used in the sliding fit
     *  @param  trackStateVector the output vector of track states
     */
    static void GetSlidingFitTrajectory(const pandora::ParticleFlowObject *const pPfo, const unsigned int layerWindow, const float layerPitch,
        std::vector<pandora::TrackState> &trackStateVector);

    /**
     *  @brief  Return track flag based on Pfo Particle ID
     *
     *  @param  pPfo the address of the Pfo
     */
    static bool IsTrack(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Return shower flag based on Pfo Particle ID
     *
     *  @param  pPfo the address of the Pfo
     */
    static bool IsShower(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get primary neutrino or antineutrino
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return pdg code of neutrino (or zero, otherwise)
     */
    static int GetPrimaryNeutrino(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Whether a pfo is a primary parent particle
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return boolean
     */
    static bool IsFinalState(const pandora::ParticleFlowObject *const pPfo);

     /**
     *  @brief  Whether a pfo is a final-state particle from a neutrino (or antineutrino) interaction
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return boolean
     */
    static bool IsNeutrinoFinalState(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Whether a pfo is a neutrino or (antineutrino)
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return boolean
     */
    static bool IsNeutrino(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Get the primary parent pfo
     *
     *  @param  pPfo the address of the Pfo
     *
     *  @return address of the primary parent pfo
     */
    static const pandora::ParticleFlowObject *GetParentPfo(const pandora::ParticleFlowObject *const pPfo);

     /**
     *  @brief  Get primary neutrino or antineutrino
     *
     *  @param   pPfo the address of the Pfo
     *
     *  @return address of primary neutrino pfo
     */
    static const pandora::ParticleFlowObject *GetParentNeutrino(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief  Sort pfos by number of constituent hits
     *
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortByNHits(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_HELPER_H
