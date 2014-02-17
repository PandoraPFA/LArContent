/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/UndershootTracksTool.cc
 * 
 *  @brief  Implementation of the undershoot tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/UndershootTracksTool.h"

using namespace pandora;

namespace lar
{

bool UndershootTracksTool::Run(const SlidingFitResultMap &slidingFitResultMap, TrackOverlapTensor &overlapTensor, ProtoParticleVector &protoParticleVector)
{
    std::cout << "UndershootTracksTool::Run() " << std::endl;

    //    if (protoParticleVector.empty())
    //        return false;
    //
    //    if (protoParticleVector.size() > 1)
    //        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    //
    //    ProtoParticle &protoParticle(protoParticleVector.at(0));
    //    this->BuildProtoParticle(ParticleComponent(pBestClusterU, pBestClusterV, pBestClusterW, bestTrackOverlapResult), protoParticle);

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

//void UndershootTracksTool::BuildProtoParticle(const ParticleComponent &firstComponent, ProtoParticle &protoParticle) const
//{
//    Cluster *pClusterU(firstComponent.GetClusterU()), *pClusterV(firstComponent.GetClusterV()), *pClusterW(firstComponent.GetClusterW());
//    protoParticle.m_clusterListU.insert(pClusterU);
//    protoParticle.m_clusterListV.insert(pClusterV);
//    protoParticle.m_clusterListW.insert(pClusterW);
//
//    ParticleComponentList particleComponentList;
//    const ClusterList &clusterListU(m_overlapTensor.GetClusterListU());
//    const ClusterList &clusterListV(m_overlapTensor.GetClusterListV());
//    const ClusterList &clusterListW(m_overlapTensor.GetClusterListW());
//
//    for (ClusterList::const_iterator iterU = clusterListU.begin(), iterUEnd = clusterListU.end(); iterU != iterUEnd; ++iterU)
//    {
//        for (ClusterList::const_iterator iterV = clusterListV.begin(), iterVEnd = clusterListV.end(); iterV != iterVEnd; ++iterV)
//        {
//            for (ClusterList::const_iterator iterW = clusterListW.begin(), iterWEnd = clusterListW.end(); iterW != iterWEnd; ++iterW)
//            {
//                try
//                {
//                    if ((pClusterU != *iterU) && (pClusterV != *iterV) && (pClusterW != *iterW))
//                        continue;
//
//                    const TrackOverlapResult &trackOverlapResult(m_overlapTensor.GetOverlapResult(*iterU, *iterV, *iterW));
//                    particleComponentList.push_back(ParticleComponent(*iterU, *iterV, *iterW, trackOverlapResult));
//                }
//                catch (StatusCodeException &)
//                {
//                }
//            }
//        }
//    }
//
//    for (ParticleComponentList::const_iterator iter = particleComponentList.begin(), iterEnd = particleComponentList.end(); iter != iterEnd; ++iter)
//    {
//        if (protoParticle.m_clusterListU.count(iter->GetClusterU()) && protoParticle.m_clusterListV.count(iter->GetClusterV()) &&
//            protoParticle.m_clusterListW.count(iter->GetClusterW()))
//        {
//            continue;
//        }
//
//        if (this->IsParticleMatch(firstComponent, *iter))
//        {
//            this->BuildProtoParticle(*iter, protoParticle);
//        }
//    }
//}
//
//------------------------------------------------------------------------------------------------------------------------------------------
//
//bool UndershootTracksTool::IsParticleMatch(const ParticleComponent &firstComponent, const ParticleComponent &secondComponent) const
//{
//    return (this->IsPossibleMatch(firstComponent.GetClusterU(), secondComponent.GetClusterU()) &&
//            this->IsPossibleMatch(firstComponent.GetClusterV(), secondComponent.GetClusterV()) &&
//            this->IsPossibleMatch(firstComponent.GetClusterW(), secondComponent.GetClusterW()));
//}
//
//------------------------------------------------------------------------------------------------------------------------------------------
//
//bool UndershootTracksTool::IsPossibleMatch(Cluster *const pFirstCluster, Cluster *const pSecondCluster) const
//{
//    if (pFirstCluster == pSecondCluster)
//        return true;
//
//    SlidingFitResultMap::const_iterator iter1 = m_slidingFitResultMap.find(pFirstCluster);
//    SlidingFitResultMap::const_iterator iter2 = m_slidingFitResultMap.find(pSecondCluster);
//
//    if ((m_slidingFitResultMap.end() == iter1) || (m_slidingFitResultMap.end() == iter2))
//        throw StatusCodeException(STATUS_CODE_FAILURE);
//
//    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult1(iter1->second);
//    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult2(iter2->second);
//
//    // Check there is no significant overlap between clusters
//    const CartesianVector minLayerPosition1(slidingFitResult1.GetGlobalMinLayerPosition());
//    const CartesianVector maxLayerPosition1(slidingFitResult1.GetGlobalMaxLayerPosition());
//    const CartesianVector minLayerPosition2(slidingFitResult2.GetGlobalMinLayerPosition());
//    const CartesianVector maxLayerPosition2(slidingFitResult2.GetGlobalMaxLayerPosition());
//
//    const CartesianVector linearDirection1((maxLayerPosition1 - minLayerPosition1).GetUnitVector());
//    const CartesianVector linearDirection2((maxLayerPosition2 - minLayerPosition2).GetUnitVector());
//
//    const float clusterOverlap1_Pos0((maxLayerPosition1 - minLayerPosition1).GetMagnitude());
//    const float clusterOverlap1_Pos1(linearDirection1.GetDotProduct(minLayerPosition2 - minLayerPosition1));
//    const float clusterOverlap1_Pos2(linearDirection1.GetDotProduct(maxLayerPosition2 - minLayerPosition1));
//
//    const float clusterOverlap2_Pos0((maxLayerPosition2 - minLayerPosition2).GetMagnitude());
//    const float clusterOverlap2_Pos1(linearDirection2.GetDotProduct(minLayerPosition1 - minLayerPosition2));
//    const float clusterOverlap2_Pos2(linearDirection2.GetDotProduct(maxLayerPosition1 - minLayerPosition2));
//
//    return ((std::min(clusterOverlap1_Pos1, clusterOverlap1_Pos2) > clusterOverlap1_Pos0 - 2.f) ||
//            (std::min(clusterOverlap2_Pos1, clusterOverlap2_Pos2) > clusterOverlap2_Pos0 - 2.f) ||
//            (std::max(clusterOverlap1_Pos1, clusterOverlap1_Pos2) < 2.f) ||
//            (std::max(clusterOverlap2_Pos1, clusterOverlap2_Pos2) < 2.f));
//}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode UndershootTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
