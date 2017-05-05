/**
 *  @file   larpandoracontent/LArStitching/StitchingPfoMergingTool.h
 * 
 *  @brief  Header file for the stitching pfo merging tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_STITCHING_PFO_MERGING_TOOL_H
#define LAR_STITCHING_PFO_MERGING_TOOL_H 1

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"
#include "larpandoracontent/LArStitching/StitchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  StitchingPfoMergingTool class
 */
class StitchingPfoMergingTool : public StitchingTool
{
public:
    void Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo);

private:
    /**
     *  @brief  Merge and delete a pair of pfos, with a specific set of conventions for cluster merging, vertex use, etc.
     * 
     *  @param  pAlgorithm the address of the responsible algorithm
     *  @param  pPfoToEnlarge the address of the pfo to enlarge
     *  @param  pPfoToDelete the address of the pfo to delete (will become a dangling pointer)
     */
    void MergeAndDeletePfos(const StitchingAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfoToEnlarge,
        const pandora::ParticleFlowObject *const pPfoToDelete) const;

    /**
     *  @brief  Select the parent cluster (same hit type and most hits) using a provided cluster list and hit type
     * 
     *  @param  clusterList the cluster list
     *  @param  hitType the hit type
     * 
     *  @return the address of the parent cluster
     */
    const pandora::Cluster *GetParentCluster(const pandora::ClusterList &clusterList, const pandora::HitType hitType) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_PFO_MERGING_TOOL_H
