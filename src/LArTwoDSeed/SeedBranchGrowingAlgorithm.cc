/**
 *  @file   LArContent/src/LArTwoDSeed/SeedBranchGrowingAlgorithm.cc
 * 
 *  @brief  Implementation of the seed growing algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArTwoDSeed/SeedBranchGrowingAlgorithm.h"

using namespace pandora;

namespace lar
{

void SeedBranchGrowingAlgorithm::GetCandidateClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (pCluster->GetNCaloHits() < 10)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

SeedBranchGrowingAlgorithm::AssociationType SeedBranchGrowingAlgorithm::AreClustersAssociated(const Cluster *const pClusterSeed, const Cluster *const pCluster) const
{
    const bool currentVertexExists(LArVertexHelper::DoesCurrentVertexExist());

    // Direction of seed cluster
    const bool isSeedForward  = (currentVertexExists) ? LArVertexHelper::IsForwardInZ(pClusterSeed) : false;
    const bool isSeedBackward = (currentVertexExists) ? LArVertexHelper::IsBackwardInZ(pClusterSeed) : false;

    const bool checkSeedForward  = (isSeedForward  || !isSeedBackward);
    const bool checkSeedBackward = (isSeedBackward || !isSeedForward);

    // Direction of candidate cluster
    const bool isCandidateForward  = (currentVertexExists) ? LArVertexHelper::IsForwardInZ(pCluster) : false;
    const bool isCandidateBackward = (currentVertexExists) ? LArVertexHelper::IsBackwardInZ(pCluster) : false;

    const bool checkCandidateForward  = (isCandidateForward  || !isCandidateBackward);
    const bool checkCandidateBackward = (isCandidateBackward || !isCandidateForward);

    // Calculate distances of association
    const float sOuter(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetOuterPseudoLayer()), pCluster));
    const float cOuter(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()), pClusterSeed));
    const float sInner(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetInnerPseudoLayer()), pCluster));
    const float cInner(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()), pClusterSeed));

    // Association check 1(a), look for enclosed clusters
    if ((cOuter < 2.5 && cInner < 2.5) && (!checkSeedForward || sInner > 2.5) && (!checkSeedBackward || sOuter > 2.5))
        return STRONG;

    // Association check 1(b), look for overlapping clusters
    if ((checkSeedForward == checkCandidateForward) && (checkSeedBackward == checkCandidateBackward))
    {
        if ((cInner < 2.5 && sOuter < 2.5) && (!checkSeedForward || sInner > 2.5) && (!checkSeedBackward || cOuter > 2.5))
            return STRONG;

        if ((cOuter < 2.5 && sInner < 2.5) && (!checkSeedBackward || sOuter > 2.5) && (!checkSeedForward || cInner > 2.5))
            return STRONG;
    }

    // Association check 2, look for branching clusters
    if ((!checkSeedForward || sInner > 10.) && (!checkSeedBackward || sOuter > 10.) && ((checkCandidateForward && cInner < 2.5) || (checkCandidateBackward && cOuter < 2.5)))
        return STANDARD;

    // Association check 3, look any distance below threshold
    //if ((sOuter < 2.5) || (cOuter < 2.5) || (sInner < 2.5) || (cInner < 2.5))
    //    return SINGLE_ORDER;

    return NONE;
}

} // namespace lar
