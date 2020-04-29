/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerFragments/MopUpRemnantsTool.h
 *
 *  @brief  Header file for the mop-up remnants tool class.
 *
 *  $Log: $
 */
#ifndef MOP_UP_REMNANTS_TOOL_H
#define MOP_UP_REMNANTS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ThreeViewRemnantsAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  MopUpRemnantsTool class
 */
class MopUpRemnantsTool : public RemnantTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    MopUpRemnantsTool();

    bool Run(ThreeViewRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Identify candidate particles
     *
     *  @param  overlapTensor the input overlap tensor
     *  @param  protoParticleVector the output vector of candidate particles
     *  @param  clusterMergeMap the output map of clusters to be merged
     */
    void FindBestShowers(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const;

    /**
     *  @brief  Copy used clusters into cluster list
     *
     *  @param  elementList the input list of elements
     *  @param  usedClusters to receive the list of used clusters
     */
    void GetUsedClusters(const TensorType::ElementList &elementList, pandora::ClusterSet &usedClusters) const;

    /**
     *  @brief  Select the best triplet of clusters
     *
     *  @param  elementList the input list of elements
     *  @param  usedClusters the list of cluster analysed so far
     *  @param  bestIter iterator to the best element in the list
     */
    void SelectBestElement(const TensorType::ElementList &elementList, const pandora::ClusterSet &usedClusters,
        TensorType::ElementList::const_iterator &bestIter) const;
};

} // namespace lar_content

#endif // #ifndef MOP_UP_REMNANTS_TOOL_H
