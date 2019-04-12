/**
 *  @file   larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h
 *
 *  @brief  Header file for the pfo mop up algorithm base class.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_MOP_UP_BASE_ALGORITHM_H
#define LAR_PFO_MOP_UP_BASE_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/MopUpBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoMopUpBaseAlgorithm class
 */
class PfoMopUpBaseAlgorithm : public MopUpBaseAlgorithm
{
public:
    /**
     *  @brief  Merge and delete a pair of pfos, with a specific set of conventions for cluster merging, vertex use, etc.
     *
     *  @param  pPfoToEnlarge the address of the pfo to enlarge
     *  @param  pPfoToDelete the address of the pfo to delete (will become a dangling pointer)
     */
    virtual void MergeAndDeletePfos(const pandora::ParticleFlowObject *const pPfoToEnlarge, const pandora::ParticleFlowObject *const pPfoToDelete) const;

    /**
     *  @brief  Select the parent cluster (same hit type and most hits) using a provided cluster list and hit type
     *
     *  @param  clusterList the cluster list
     *  @param  hitType the hit type
     *
     *  @return the address of the parent cluster
     */
    static const pandora::Cluster *GetParentCluster(const pandora::ClusterList &clusterList, const pandora::HitType hitType);

protected:
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_MOP_UP_BASE_ALGORITHM_H
