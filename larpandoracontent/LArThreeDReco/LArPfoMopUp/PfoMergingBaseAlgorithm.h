/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/PfoMergingBaseAlgorithm.h
 * 
 *  @brief  Header file for the pfo merging algorithm base class.
 * 
 *  $Log: $
 */
#ifndef LAR_PFO_MERGING_BASE_ALGORITHM_H
#define LAR_PFO_MERGING_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoMergingBaseAlgorithm class
 */
class PfoMergingBaseAlgorithm : public pandora::Algorithm
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
    virtual const pandora::Cluster *GetParentCluster(const pandora::ClusterList &clusterList, const pandora::HitType hitType) const;

    /**
     *  @brief  Find the name of the list hosting a specific object
     * 
     *  @param  pT the address of the object
     * 
     *  @return the name of the list
     */
    template <typename T>
    const std::string GetListName(const T *const pT) const;

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_daughterListNames;                ///< The list of potential daughter object list names
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_MERGING_BASE_ALGORITHM_H
