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

namespace lar
{

/**
 *  @brief  LArPfoHelper class
 */
class LArPfoHelper
{
public:

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
     *  @brief  Sort pfos by number of constituent hits
     * 
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortByNHits(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

    /**
     *  @brief  Read the vertex helper settings
     *
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar

#endif // #ifndef LAR_PFO_HELPER_H
