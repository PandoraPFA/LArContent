/**
 *  @file   LArContent/src/LArControlFlow/StitchingCosmicRayMergingTool.cc
 *
 *  @brief  Implementation of the stitching cosmic ray merging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "StitchingCosmicRayMergingTool.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>

#include "larpandoracontent/LArControlFlow/StitchingCosmicRayMergingTool.h"

using namespace pandora;

namespace lar_content
{

StitchingCosmicRayMergingTool::StitchingCosmicRayMergingTool() :
    m_useXcoordinate(false),
    m_alwaysApplyT0Calculation(true),
    m_stitchPfosZMatches(false),
    m_maxXYdisplacementPfoZmatching(1.f),
    m_maxZdistancePfoZmatching(10.f),
    m_halfWindowLayers(30),
    m_minLengthSquared(50.f),
    m_minCosRelativeAngle(0.966),
    m_relaxMinLongitudinalDisplacement(-5.f),
    m_maxLongitudinalDisplacementX(15.f),
    m_maxTransverseDisplacement(5.f),
    m_relaxCosRelativeAngle(0.906),
    m_relaxTransverseDisplacement(2.5f),
    m_minNCaloHits3D(0),
    m_maxX0FractionalDeviation(0.3f),
    m_boundaryToleranceWidth(10.f)
{
}

void StitchingCosmicRayMergingTool::Run(const MasterAlgorithm *const pAlgorithm, const PfoList *const pMultiPfoList,
    PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    if (this->GetPandora().GetGeometry()->GetLArTPCMap().size() < 2)
        return;

    if (pfoToLArTPCMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    PfoList primaryPfos;
    this->SelectPrimaryPfos(pMultiPfoList, pfoToLArTPCMap, primaryPfos);

    ThreeDPointingClusterMap pointingClusterMap;
    this->BuildPointingClusterMaps(primaryPfos, pfoToLArTPCMap, pointingClusterMap);

    LArTPCToPfoMap larTPCToPfoMap;
    this->BuildTPCMaps(primaryPfos, pfoToLArTPCMap, larTPCToPfoMap);

    // std::cout << larTPCToPfoMap.size() << "\n";
    // for (const auto & pair : larTPCToPfoMap)
    // {
    //   std::cout << "TPC " << pair.first->GetLArTPCVolumeId() 
    //             <<" has " <<pair.second.size()
    //             <<" pfos \n";
    //   for (const auto &pfo_in_TPC : pair.second)
    //   {
    //
    //     ThreeDPointingClusterMap::const_iterator iter1 = pointingClusterMap.find(pfo_in_TPC);
    //     if (iter1 == pointingClusterMap.end() ) continue;
    //     const LArPointingCluster &pointingCluster(iter1->second);
    //     const LArPointingCluster::Vertex inner = pointingCluster.GetInnerVertex();
    //     const LArPointingCluster::Vertex outer = pointingCluster.GetOuterVertex();
    //     auto pfo_in_TPC_vertex = pfo_in_TPC->GetVertexList().begin();
    //     if (pfo_in_TPC_vertex == pfo_in_TPC->GetVertexList().end()) continue;
    //     const pandora::Vertex *vertex = *pfo_in_TPC_vertex;
    //
    //     std::cout << "{ \"pfo_id\": " << pfo_in_TPC
    //       << ", \"tpc_id\": " << pair.first->GetLArTPCVolumeId()
    //       << ", \"vertex_x\": " << vertex->GetPosition().GetX()
    //       << ", \"vertex_y\": " << vertex->GetPosition().GetY()  
    //       << ", \"vertex_z\": " << vertex->GetPosition().GetZ()
    //       << ", \"inner_x\": " << inner.GetPosition().GetX()
    //       << ", \"inner_y\": " << inner.GetPosition().GetY()
    //       << ", \"inner_z\": " << inner.GetPosition().GetZ()
    //       << ", \"outer_x\": " << outer.GetPosition().GetX()
    //       << ", \"outer_y\": " << outer.GetPosition().GetY()
    //       << ", \"outer_z\": " << outer.GetPosition().GetZ()
    //       << " },\n";
    //   }
    // }

    // Functions to enable Z stitching
    PfoAssociationMatrix pfoAssociationMatrixZ;
    this->CreatePfoMatchesZ(larTPCToPfoMap, pointingClusterMap, pfoAssociationMatrixZ);
   
    PfoAssociationMatrix pfoBestMatchesZ;
    this->SelectBestMatchesZ(pfoAssociationMatrixZ, pfoBestMatchesZ);

    PfoAssociationMatrix pfoUniqueMatchesZ;
    this->SelectUniqueMatchesZ(pfoBestMatchesZ, pfoUniqueMatchesZ);

    std::vector<PfoVector> pfosGroupedAlongZ;
    this->GroupPfoAlongZ(pfoUniqueMatchesZ, pfosGroupedAlongZ);

    this->SortGroupedPfosZ(pfosGroupedAlongZ, pointingClusterMap);
    // end functions to enable Z stitching

    PfoAssociationMatrix pfoAssociationMatrix;
    this->CreatePfoMatches(larTPCToPfoMap, pointingClusterMap, pfoAssociationMatrix);

    PfoMergeMap pfoSelectedMatches;
    this->SelectPfoMatches(pfoAssociationMatrix, pfoSelectedMatches);

    PfoMergeMap pfoSelectedMerges;
    this->SelectPfoMerges(pfoSelectedMatches, pfoSelectedMerges);

    PfoMergeMap pfoOrderedMerges;
    this->OrderPfoMerges(pfoToLArTPCMap, pointingClusterMap, pfoSelectedMerges, pfoOrderedMerges);

    this->StitchPfos(pAlgorithm, pointingClusterMap, pfoOrderedMerges, pfoToLArTPCMap, stitchedPfosToX0Map, pfosGroupedAlongZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SelectPrimaryPfos(const PfoList *pInputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, PfoList &outputPfoList) const
{
    for (const ParticleFlowObject *const pPfo : *pInputPfoList)
    {
        ClusterList clusterList;
        LArPfoHelper::GetThreeDClusterList(pPfo, clusterList);
        if(clusterList.size()==0) continue;
        CaloHitList caloHitList;
        clusterList.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        if (!LArPfoHelper::IsFinalState(pPfo) || !LArPfoHelper::IsTrack(pPfo))
            continue;

        if (!pfoToLArTPCMap.count(pPfo))
            continue;

        if (caloHitList.size()<2) // avoid crash when performing the fit on clusters with 1 hit
        {
          continue;
        }

        outputPfoList.push_back(pPfo);
    }

    outputPfoList.sort(LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::BuildPointingClusterMaps(
    const PfoList &inputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, ThreeDPointingClusterMap &pointingClusterMap) const
{
    for (const ParticleFlowObject *const pPfo : inputPfoList)
    {
        try
        {
            PfoToLArTPCMap::const_iterator tpcIter(pfoToLArTPCMap.find(pPfo));

            if (pfoToLArTPCMap.end() == tpcIter)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            const float slidingFitPitch(tpcIter->second->GetWirePitchW());

            ClusterList clusterList;
            LArPfoHelper::GetThreeDClusterList(pPfo, clusterList);

            if (1 != clusterList.size())
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            const ThreeDSlidingFitResult slidingFitResult(clusterList.front(), m_halfWindowLayers, slidingFitPitch);
            (void)pointingClusterMap.insert(ThreeDPointingClusterMap::value_type(pPfo, LArPointingCluster(slidingFitResult)));

        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::BuildTPCMaps(const PfoList &inputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, LArTPCToPfoMap &larTPCToPfoMap) const
{
    for (const ParticleFlowObject *const pPfo : inputPfoList)
    {
        PfoToLArTPCMap::const_iterator iter(pfoToLArTPCMap.find(pPfo));

        if (pfoToLArTPCMap.end() != iter)
            larTPCToPfoMap[iter->second].push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CreatePfoMatchesZ(const LArTPCToPfoMap &larTPCToPfoMap, const ThreeDPointingClusterMap &pointingClusterMap, PfoAssociationMatrix &pfoAssociationMatrix) const
{

    LArTPCVector larTPCVector;

    if(larTPCToPfoMap.empty()) return;

    for (const auto &mapEntry : larTPCToPfoMap) larTPCVector.push_back(mapEntry.first); 
    
    // Order TPCs in ascending order firs on their x coordinates than y and z. 
    std::sort(larTPCVector.begin(), larTPCVector.end(), LArStitchingHelper::SortTPCs); 

    for (LArTPCVector::const_iterator tpcIter = larTPCVector.begin(), tpcIterEnd = larTPCVector.end(); tpcIter != tpcIterEnd-1; ++tpcIter)
    {
        const LArTPC *const pLArTPC1(*tpcIter);
        const LArTPC *const pLArTPC2(*(tpcIter+1));

        const PfoList &pfoList1(larTPCToPfoMap.at(pLArTPC1));
        const PfoList &pfoList2(larTPCToPfoMap.at(pLArTPC2));

        std::cout << "\n";
        std::cout << "_____________ considered pairs of tpc (" 
          << pLArTPC1->GetLArTPCVolumeId() << ", " <<pLArTPC1->GetLArTPCVolumeId()<<")____________________ \n";

        // find pfo z matches within the same tpc 
        this->CreatePfoMatchesZ(0.f, pfoList1, pfoList1, pointingClusterMap, pfoAssociationMatrix);

        std::cout << "\n";
        std::cout << "_____________ considered pairs of tpc (" 
          << pLArTPC1->GetLArTPCVolumeId() << ", " <<pLArTPC2->GetLArTPCVolumeId()<<")____________________ \n";

        if ( fabs(pLArTPC1->GetCenterY() - pLArTPC2->GetCenterY()) > std::numeric_limits<float>::epsilon() ) continue;
        if ( fabs(pLArTPC1->GetCenterX() - pLArTPC2->GetCenterX()) > std::numeric_limits<float>::epsilon() ) continue;

        // find pfo z matches with the clostest TPC along z
        const float TPCsGap = fabs(pLArTPC1->GetCenterZ() - pLArTPC2->GetCenterZ()) - (0.5f*pLArTPC1->GetWidthZ() + 0.5f*pLArTPC2->GetWidthZ());
        this->CreatePfoMatchesZ(TPCsGap, pfoList1, pfoList2, pointingClusterMap, pfoAssociationMatrix);
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CreatePfoMatchesZ(const float TPCsGap, const PfoList& pfoList1, const PfoList& pfoList2, const ThreeDPointingClusterMap &pointingClusterMap, PfoAssociationMatrix &pfoAssociationMatrix) const
{

      for (const ParticleFlowObject *const pPfo1 : pfoList1){

          ThreeDPointingClusterMap::const_iterator iter1 = pointingClusterMap.find(pPfo1);
          if(iter1 == pointingClusterMap.end()) continue;
          const LArPointingCluster &pointingCluster1(iter1->second);

          if(pointingCluster1.GetLengthSquared() < m_minLengthSquared) continue;

          CaloHitList caloHitList3D1;
          LArPfoHelper::GetCaloHits(pPfo1, TPC_3D, caloHitList3D1);
          
          if (caloHitList3D1.size() < m_minNCaloHits3D) continue;

          const LArPointingCluster::Vertex outer = pointingCluster1.GetOuterVertex();
          const LArPointingCluster::Vertex inner1 = pointingCluster1.GetInnerVertex();

          // if pointingCluster outer vertex is far from the tpc z boundary
          // we don't need to match it
          //
          // TODO come up with non-hardcoded cut
          // if ( fabs(outer.GetPosition().GetZ() - (pLArTPC1->GetCenterZ() + 0.5*pLArTPC1->GetWidthZ())) > 20) continue;

          for (const ParticleFlowObject *const pPfo2 : pfoList2)
          {

              ThreeDPointingClusterMap::const_iterator iter2 = pointingClusterMap.find(pPfo2);
              if(iter2 == pointingClusterMap.end()) continue;
              const LArPointingCluster &pointingCluster2(iter2->second);

              if(pointingCluster2.GetLengthSquared() < m_minLengthSquared) continue;

              CaloHitList caloHitList3D2;
              LArPfoHelper::GetCaloHits(pPfo1, TPC_3D, caloHitList3D2);

              if (caloHitList3D2.size() < m_minNCaloHits3D) continue;

              const LArPointingCluster::Vertex inner = pointingCluster2.GetInnerVertex();
              const LArPointingCluster::Vertex outer2 = pointingCluster2.GetOuterVertex();

              std::cout << "\n";
              std::cout<< "new pair of pfo considered:\n";
              std::cout << pPfo1
                        // <<", inner vertex (x = " << inner1.GetPosition().GetX() 
                        //               << ", y = " << inner1.GetPosition().GetY()
                        //               << ", z = " << inner1.GetPosition().GetZ() 
                        // <<"), outer vertex (x = " << outer.GetPosition().GetX() 
                        //               << ", y = " << outer.GetPosition().GetY()
                        //               << ", z = " << outer.GetPosition().GetZ() 
                                      <<")\n";    
              std::cout << pPfo2
                        // <<", inner vertex (x = " << inner.GetPosition().GetX() 
                        //               << ", y = " << inner.GetPosition().GetY()
                        //               << ", z = " << inner.GetPosition().GetZ() 
                        // <<"), outer vertex (x = " << outer2.GetPosition().GetX() 
                        //               << ", y = " << outer2.GetPosition().GetY()
                        //               << ", z = " << outer2.GetPosition().GetZ() 
                                      <<")\n";

             LArPointingCluster::Vertex pointingVertex1, pointingVertex2;

             LArStitchingHelper::GetClosestVertices(pointingCluster1, pointingCluster2, pointingVertex1, pointingVertex2);

              if (this->DoVerticesMatchZ(TPCsGap, pointingVertex1, pointingVertex2))
              {
                std::cout << "vertices match \n";
                // Store this association
                const PfoAssociation::VertexType vertexType1(pointingVertex1.IsInnerVertex() ? PfoAssociation::INNER 
                                                                                             : PfoAssociation::OUTER);
                const PfoAssociation::VertexType vertexType2(pointingVertex2.IsInnerVertex() ? PfoAssociation::INNER 
                                                                                             : PfoAssociation::OUTER);

                const float particleLength1(pointingCluster1.GetLengthSquared());
                const float particleLength2(pointingCluster2.GetLengthSquared());

                pfoAssociationMatrix[pPfo1].insert(PfoAssociationMap::value_type(pPfo2, PfoAssociation(vertexType1, vertexType2, particleLength2)));
                pfoAssociationMatrix[pPfo2].insert(PfoAssociationMap::value_type(pPfo1, PfoAssociation(vertexType2, vertexType1, particleLength1)));
              }
              } // pfo2 loop 
          } // pfo 1 loop
}
//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CreatePfoMatches(const LArTPCToPfoMap &larTPCToPfoMap,
    const ThreeDPointingClusterMap &pointingClusterMap, PfoAssociationMatrix &pfoAssociationMatrix) const
{
    LArTPCVector larTPCVector;
    for (const auto &mapEntry : larTPCToPfoMap)
        larTPCVector.push_back(mapEntry.first);
    std::sort(larTPCVector.begin(), larTPCVector.end(), LArStitchingHelper::SortTPCs);

    for (LArTPCVector::const_iterator tpcIter1 = larTPCVector.begin(), tpcIterEnd = larTPCVector.end(); tpcIter1 != tpcIterEnd; ++tpcIter1)
    {
        const LArTPC *const pLArTPC1(*tpcIter1);
        const PfoList &pfoList1(larTPCToPfoMap.at(pLArTPC1));

        for (LArTPCVector::const_iterator tpcIter2 = tpcIter1; tpcIter2 != tpcIterEnd; ++tpcIter2)
        {
            const LArTPC *const pLArTPC2(*tpcIter2);
            const PfoList &pfoList2(larTPCToPfoMap.at(pLArTPC2));

            if (!LArStitchingHelper::CanTPCsBeStitched(*pLArTPC1, *pLArTPC2))
                continue;

            for (const ParticleFlowObject *const pPfo1 : pfoList1)
            {
                for (const ParticleFlowObject *const pPfo2 : pfoList2)
                    this->CreatePfoMatches(*pLArTPC1, *pLArTPC2, pPfo1, pPfo2, pointingClusterMap, pfoAssociationMatrix);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool StitchingCosmicRayMergingTool::DoVerticesMatchZ(
                                                     const float TPCsGap, 
                                                     const LArPointingCluster::Vertex &pointingVertex1,
                                                     const LArPointingCluster::Vertex &pointingVertex2) const
{
    // Input pointing vertices from adjacent TPCs match across z if:
    // - they are not more than m_maxZdistancePfoZmatching apart along Z
    // - they are colinear
    // - their projections onto a plane at z (in the middle of the TPC) are consistent

    const CartesianVector vertex1 = pointingVertex1.GetPosition(); 
    const CartesianVector vertex2 = pointingVertex2.GetPosition(); 
    const CartesianVector direction1 = pointingVertex1.GetDirection(); 
    const CartesianVector direction2 = pointingVertex2.GetDirection(); 

    // Pointing vertices must not be further than m_maxZdistancePfoZmatching (not including the gap between TPCs' boundaries)
    float maxPfoZdistance = m_maxZdistancePfoZmatching + TPCsGap;
    if (fabs(vertex2.GetZ() - vertex1.GetZ()) > maxPfoZdistance) return false;

    // Pointing clusters must have a parallel direction
    const float cosRelativeAngle(-direction1.GetDotProduct(direction2));

    if (cosRelativeAngle < m_relaxCosRelativeAngle) return false;

    const CartesianVector diff = {vertex2.GetX() - vertex2.GetX(), vertex2.GetY() - vertex1.GetY(), vertex2.GetZ() - vertex1.GetZ()};
    const CartesianVector middle_point = vertex1 + diff.GetUnitVector() * 0.5f * diff.GetMagnitude();
    const float z_plane = middle_point.GetZ();
 
    const float t1 = (z_plane - vertex1.GetZ())/direction1.GetZ();    
    const float x1 = vertex1.GetX() + t1*direction1.GetX(); 
    const float y1 = vertex1.GetY() + t1*direction1.GetY(); 

    const float t2 = (z_plane - vertex2.GetZ())/direction2.GetZ();    
    const float x2 = vertex2.GetX() + t2*direction2.GetX(); 
    const float y2 = vertex2.GetY() + t2*direction2.GetY(); 
    
    const float relative_XYdisplacement = std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

    // Pointing clusters projected onto XY planes must be not further than m_maxXYdisplacementPfoZmatching 
    if(relative_XYdisplacement > m_maxXYdisplacementPfoZmatching) return false;

    return true;
} 

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CreatePfoMatches(const LArTPC &larTPC1, const LArTPC &larTPC2, const ParticleFlowObject *const pPfo1,
    const ParticleFlowObject *const pPfo2, const ThreeDPointingClusterMap &pointingClusterMap, PfoAssociationMatrix &pfoAssociationMatrix) const
{
    // Get centre and width of boundary between tpcs
    const float boundaryCenterX(LArStitchingHelper::GetTPCBoundaryCenterX(larTPC1, larTPC2));
    const float boundaryWidthX(LArStitchingHelper::GetTPCBoundaryWidthX(larTPC1, larTPC2));
    const float maxLongitudinalDisplacementX(m_maxLongitudinalDisplacementX + boundaryWidthX);

    // Get the pointing cluster corresponding to each of these Pfos
    ThreeDPointingClusterMap::const_iterator iter1 = pointingClusterMap.find(pPfo1);
    ThreeDPointingClusterMap::const_iterator iter2 = pointingClusterMap.find(pPfo2);

    if (pointingClusterMap.end() == iter1 || pointingClusterMap.end() == iter2)
        return;

    const LArPointingCluster &pointingCluster1(iter1->second);
    const LArPointingCluster &pointingCluster2(iter2->second);

    // Check length of pointing clusters
    if (pointingCluster1.GetLengthSquared() < m_minLengthSquared || pointingCluster2.GetLengthSquared() < m_minLengthSquared)
        return;

    // Get number of 3D hits in each of the pfos
    CaloHitList caloHitList3D1;
    LArPfoHelper::GetCaloHits(pPfo1, TPC_3D, caloHitList3D1);

    CaloHitList caloHitList3D2;
    LArPfoHelper::GetCaloHits(pPfo2, TPC_3D, caloHitList3D2);

    // Check number of 3D hits in each of the pfos
    if (caloHitList3D1.size() < m_minNCaloHits3D || caloHitList3D2.size() < m_minNCaloHits3D)
        return;

    // Get closest pair of vertices
    LArPointingCluster::Vertex pointingVertex1, pointingVertex2;

    try
    {
        LArStitchingHelper::GetClosestVertices(larTPC1, larTPC2, pointingCluster1, pointingCluster2, pointingVertex1, pointingVertex2);
    }
    catch (const pandora::StatusCodeException &)
    {
        return;
    }

    // Pointing clusters must have a parallel direction
    const float cosRelativeAngle(-pointingVertex1.GetDirection().GetDotProduct(pointingVertex2.GetDirection()));

    if (cosRelativeAngle < m_relaxCosRelativeAngle)
        return;
    // Pointing clusters must have a non-zero X direction (so that they point across drift volume boundary)
    const float pX1(std::fabs(pointingVertex1.GetDirection().GetX()));
    const float pX2(std::fabs(pointingVertex2.GetDirection().GetX()));

    if (pX1 < std::numeric_limits<float>::epsilon() || pX2 < std::numeric_limits<float>::epsilon())
        return;

    // Pointing clusters must intersect at a drift volume boundary
    const float intersectX(0.5 * (pointingVertex1.GetPosition().GetX() + pointingVertex2.GetPosition().GetX()));

    if (std::fabs(intersectX - boundaryCenterX) > maxLongitudinalDisplacementX)
        return;

    // Impact parameters
    float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);

    try
    {
        if (m_useXcoordinate)
        {
            LArPointingClusterHelper::GetImpactParameters(pointingVertex1, pointingVertex2, rL1, rT1);
            LArPointingClusterHelper::GetImpactParameters(pointingVertex2, pointingVertex1, rL2, rT2);
        }
        else
        {
            LArPointingClusterHelper::GetImpactParametersInYZ(pointingVertex1, pointingVertex2, rL1, rT1);
            LArPointingClusterHelper::GetImpactParametersInYZ(pointingVertex2, pointingVertex1, rL2, rT2);
        }
    }
    catch (const pandora::StatusCodeException &)
    {
        return;
    }

    const float minL((!LArGeometryHelper::IsInGap(this->GetPandora(), pointingVertex1.GetPosition(), TPC_3D) ||
                         !LArGeometryHelper::IsInGap(this->GetPandora(), pointingVertex2.GetPosition(), TPC_3D))
            ? -1.f
            : m_relaxMinLongitudinalDisplacement);
    const float dXdL1(m_useXcoordinate                                  ? pX1
            : (1.f - pX1 * pX1 > std::numeric_limits<float>::epsilon()) ? pX1 / std::sqrt(1.f - pX1 * pX1)
                                                                        : minL);
    const float dXdL2(m_useXcoordinate                                  ? pX2
            : (1.f - pX2 * pX2 > std::numeric_limits<float>::epsilon()) ? pX2 / std::sqrt(1.f - pX2 * pX2)
                                                                        : minL);
    const float maxL1(maxLongitudinalDisplacementX / dXdL1);
    const float maxL2(maxLongitudinalDisplacementX / dXdL2);

    if (rL1 < minL || rL1 > maxL1 || rL2 < minL || rL2 > maxL2)
        return;

    // Selection cuts on transverse impact parameters
    const bool minPass(std::min(rT1, rT2) < m_relaxTransverseDisplacement && cosRelativeAngle > m_relaxCosRelativeAngle);
    const bool maxPass(std::max(rT1, rT2) < m_maxTransverseDisplacement && cosRelativeAngle > m_minCosRelativeAngle);

    if (!minPass && !maxPass)
        return;

    // Store this association
    const PfoAssociation::VertexType vertexType1(pointingVertex1.IsInnerVertex() ? PfoAssociation::INNER : PfoAssociation::OUTER);
    const PfoAssociation::VertexType vertexType2(pointingVertex2.IsInnerVertex() ? PfoAssociation::INNER : PfoAssociation::OUTER);

    const float particleLength1(pointingCluster1.GetLengthSquared());
    const float particleLength2(pointingCluster2.GetLengthSquared());

    pfoAssociationMatrix[pPfo1].insert(PfoAssociationMap::value_type(pPfo2, PfoAssociation(vertexType1, vertexType2, particleLength2)));
    pfoAssociationMatrix[pPfo2].insert(PfoAssociationMap::value_type(pPfo1, PfoAssociation(vertexType2, vertexType1, particleLength1)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool StitchingCosmicRayMergingTool::GetBestMatch(const PfoAssociationMatrix &bestAssociationMatrix, const ParticleFlowObject*& current_pfo, PfoAssociation::VertexType vertexType){
  
  if(!current_pfo) 
    return false;

  const auto mapIt = bestAssociationMatrix.find(current_pfo);

  if (mapIt == bestAssociationMatrix.end())
    return false;

  const PfoAssociationMap &associationMap = mapIt->second;

  // TODO : check input PfoAssociationMap as only two keys: best inner, best outer
  for(const auto &entry: associationMap)
  {
    const ParticleFlowObject* pfo = entry.first;
    const PfoAssociation& association = entry.second;
    
    if(association.GetDaughter() == vertexType)
    {
      current_pfo = pfo;
      return true;
    }
  }
  current_pfo = nullptr;
  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::GroupPfoAlongZ(const PfoAssociationMatrix &bestAssociationMatrix, std::vector<PfoVector>& pfosGroupedAlongZ){

    PfoVector alreadyProcessed;
    for (const auto &mapEntry : bestAssociationMatrix)
    {
        // this will contain the list of pfo to be merged to pPfo
        // based on the best inner/outer matches
        PfoVector pfosToBeMerged;

        const ParticleFlowObject* this_pfo = mapEntry.first;

        auto it = std::find(alreadyProcessed.begin(), alreadyProcessed.end(), this_pfo);
        
        if(it != alreadyProcessed.end()){
          continue; // skip pfo already processed
        }else{
          pfosToBeMerged.push_back(this_pfo);
          alreadyProcessed.push_back(this_pfo);
        }

        // ptr to the best inner match for the current pfo
        const ParticleFlowObject* inner = this_pfo;
        // ptr to the best outer match for the current pfo
        const ParticleFlowObject* outer = this_pfo;

        while(GetBestMatch(bestAssociationMatrix, outer, PfoAssociation::VertexType::OUTER))
        {
          pfosToBeMerged.push_back(outer);
          alreadyProcessed.push_back(outer);
        }

        while(GetBestMatch(bestAssociationMatrix, inner, PfoAssociation::VertexType::INNER))
        {
          pfosToBeMerged.push_back(inner);
          alreadyProcessed.push_back(inner);
        }

        pfosGroupedAlongZ.push_back(pfosToBeMerged);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SortGroupedPfosZ(std::vector<PfoVector>& pfosGroupedAlongZ, const ThreeDPointingClusterMap &pointingClusterMap){

    // loop over groups of pfo
    for (PfoVector& pfoVector: pfosGroupedAlongZ)
    {
      std::sort(pfoVector.begin(), pfoVector.end(), 
          [&pointingClusterMap](const ParticleFlowObject* pfo1, const ParticleFlowObject* pfo2)
          {
            const LArPointingCluster& p1 = pointingClusterMap.at(pfo1);
            const LArPointingCluster& p2 = pointingClusterMap.at(pfo2);

            return (p1.GetInnerVertex().GetPosition().GetZ() < p2.GetInnerVertex().GetPosition().GetZ());
          });
     }
  }

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SelectBestMatchesZ(const PfoAssociationMatrix &pfoAssociationMatrix, PfoAssociationMatrix &bestAssociationMatrix) const
{

    PfoVector pfos;
    for (const auto &mapEntry : pfoAssociationMatrix)
        pfos.push_back(mapEntry.first);
    std::sort(pfos.begin(), pfos.end(), LArPfoHelper::SortByNHits);

    // loop over pfos
    for (const ParticleFlowObject *const pPfo : pfos)
    {
        const PfoAssociationMap &pfoAssociationMap(pfoAssociationMatrix.at(pPfo));

        const ParticleFlowObject *pBestPfoInner(nullptr);
        const ParticleFlowObject *pBestPfoOuter(nullptr);
        PfoAssociation bestAssociationInner(PfoAssociation::UNDEFINED, PfoAssociation::UNDEFINED, 0.f);
        PfoAssociation bestAssociationOuter(PfoAssociation::UNDEFINED, PfoAssociation::UNDEFINED, 0.f);

        PfoVector matchedPfos;
        for (const auto &mapEntry : pfoAssociationMap)
            matchedPfos.push_back(mapEntry.first);
        std::sort(matchedPfos.begin(), matchedPfos.end(), LArPfoHelper::SortByNHits);

        // loop over matched pfo
        for (const ParticleFlowObject *const matchedPfo : matchedPfos)
        {
            const PfoAssociation &thisPfoAssociation(pfoAssociationMap.at(matchedPfo));

            if (thisPfoAssociation.GetParent() == PfoAssociation::INNER &&
                thisPfoAssociation.GetFigureOfMerit() > bestAssociationInner.GetFigureOfMerit())
            {
              bestAssociationInner = thisPfoAssociation;
              pBestPfoInner = matchedPfo;

            }else if (thisPfoAssociation.GetParent() == PfoAssociation::OUTER &&
                      thisPfoAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit()){
              bestAssociationOuter = thisPfoAssociation;
              pBestPfoOuter = matchedPfo;
            }else{
              continue;
            }

         }

         (void)bestAssociationMatrix[pPfo].insert(PfoAssociationMap::value_type(pBestPfoInner, bestAssociationInner));

         (void)bestAssociationMatrix[pPfo].insert(PfoAssociationMap::value_type(pBestPfoOuter, bestAssociationOuter));

         std::cout << "PFO : " << pPfo << ", best inner : " << pBestPfoInner << ", best outer : " << pBestPfoOuter  << "\n";
         std::cout << "is inner null " << (pBestPfoInner==nullptr) << ", is outer null " << (pBestPfoOuter==nullptr) << "\n";

      }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<const pandora::ParticleFlowObject *, const pandora::ParticleFlowObject *> StitchingCosmicRayMergingTool::GetInnerOuterDaughters(const pandora::ParticleFlowObject *pfo, const PfoAssociationMatrix& associationMatrix) const 
{
    const pandora::ParticleFlowObject *inner = nullptr;
    const pandora::ParticleFlowObject *outer = nullptr;

    const auto it = associationMatrix.find(pfo);
    if(it == associationMatrix.end())
      return {nullptr, nullptr};

    const auto &associationMap = it->second;
    if(associationMap.size() != 2)
      return {nullptr, nullptr};

    for (const auto &[daughter, assoc] : associationMap)
    {
      if(assoc.GetDaughter() == PfoAssociation::INNER)
        inner = daughter;
      else
        outer = daughter;
    }

    return {inner, outer};
}

void StitchingCosmicRayMergingTool::SelectUniqueMatchesZ(const PfoAssociationMatrix &bestAssociationMatrix, PfoAssociationMatrix &uniqueAssociationMatrix) const
{
    // P1 = parent1 , P11 = parent1 best inner match, P12 = parent1 best outer match
    // P2 = parent2,  P21 = parent2 best inner match, P22 = parent2 best outer match
    // unique match between P1 and P2 if : (P12 == P2 && P21 == P1) || (P11 == P2 && P22 == P1) 
    PfoVector pParentPfos;

    for (const auto &mapEntry : bestAssociationMatrix)
        pParentPfos.push_back(mapEntry.first);
    std::sort(pParentPfos.begin(), pParentPfos.end(), LArPfoHelper::SortByNHits);

    // loop over parents
    for (const ParticleFlowObject *const P1 : pParentPfos)
    {
      const PfoAssociationMap &parentAssociationMap(bestAssociationMatrix.at(P1));

      if(parentAssociationMap.size()!=2)
      {
        std::cout << "pfo " << P1 << " association map should have two entries (inner and outer). Skipping this pfo \n";
        continue;
      }else
      {
        // P11 = parent1 inner match
        // P12 = parent1 outer match
        auto [P11, P12] = GetInnerOuterDaughters(P1, bestAssociationMatrix);

        if(P11)
        {
            auto P2 = P11; 
            // P21 = parent2 inner match
            // P22 = parent2 outer match
            auto [P21, P22] = GetInnerOuterDaughters(P2, bestAssociationMatrix); 
            if ((P22 == P1) || ((P21 == P1) && (P12 == P2)))
            {
                const PfoAssociation &P11Association(parentAssociationMap.at(P11));
                (void)uniqueAssociationMatrix[P1].insert(PfoAssociationMap::value_type(P11, P11Association));
            }
        }
        
        if(P12)
        {
            auto P2 = P12; 
            // P21 = parent2 inner match
            // P22 = parent2 outer match
            auto [P21, P22] = GetInnerOuterDaughters(P2, bestAssociationMatrix); 
            if ((P21 == P1) || ((P21 == P1) && (P12 == P2)))
            {
                const PfoAssociation &P12Association(parentAssociationMap.at(P12));
                (void)uniqueAssociationMatrix[P1].insert(PfoAssociationMap::value_type(P12, P12Association));
            }
        }
      }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
void StitchingCosmicRayMergingTool::SelectPfoMatches(const PfoAssociationMatrix &pfoAssociationMatrix, PfoMergeMap &pfoMatches) const
{
    // First step: loop over association matrix and find best associations A -> X and B -> Y
    // =====================================================================================
    PfoAssociationMatrix bestAssociationMatrix;

    PfoVector pfoVector1;
    for (const auto &mapEntry : pfoAssociationMatrix)
        pfoVector1.push_back(mapEntry.first);
    std::sort(pfoVector1.begin(), pfoVector1.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo1 : pfoVector1)
    {
        const PfoAssociationMap &pfoAssociationMap(pfoAssociationMatrix.at(pPfo1));

        const ParticleFlowObject *pBestPfoInner(nullptr);
        PfoAssociation bestAssociationInner(PfoAssociation::UNDEFINED, PfoAssociation::UNDEFINED, 0.f);

        const ParticleFlowObject *pBestPfoOuter(nullptr);
        PfoAssociation bestAssociationOuter(PfoAssociation::UNDEFINED, PfoAssociation::UNDEFINED, 0.f);

        PfoVector pfoVector2;
        for (const auto &mapEntry : pfoAssociationMap)
            pfoVector2.push_back(mapEntry.first);
        std::sort(pfoVector2.begin(), pfoVector2.end(), LArPfoHelper::SortByNHits);

        for (const ParticleFlowObject *const pPfo2 : pfoVector2)
        {
            const PfoAssociation &pfoAssociation(pfoAssociationMap.at(pPfo2));

            // Inner associations
            if (pfoAssociation.GetParent() == PfoAssociation::INNER)
            {
                if (pfoAssociation.GetFigureOfMerit() > bestAssociationInner.GetFigureOfMerit())
                {
                    bestAssociationInner = pfoAssociation;
                    pBestPfoInner = pPfo2;
                }
            }

            // Outer associations
            if (pfoAssociation.GetParent() == PfoAssociation::OUTER)
            {
                if (pfoAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit())
                {
                    bestAssociationOuter = pfoAssociation;
                    pBestPfoOuter = pPfo2;
                }
            }
        }

        if (pBestPfoInner)
            (void)bestAssociationMatrix[pPfo1].insert(PfoAssociationMap::value_type(pBestPfoInner, bestAssociationInner));

        if (pBestPfoOuter)
            (void)bestAssociationMatrix[pPfo1].insert(PfoAssociationMap::value_type(pBestPfoOuter, bestAssociationOuter));
    }

    // Second step: make the merge if A -> X and B -> Y is in fact A -> B and B -> A
    // =============================================================================
    PfoVector pfoVector3;
    for (const auto &mapEntry : bestAssociationMatrix)
        pfoVector3.push_back(mapEntry.first);
    std::sort(pfoVector3.begin(), pfoVector3.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pParentPfo : pfoVector3)
    {
        const PfoAssociationMap &parentAssociationMap(bestAssociationMatrix.at(pParentPfo));

        PfoVector pfoVector4;
        for (const auto &mapEntry : parentAssociationMap)
            pfoVector4.push_back(mapEntry.first);
        std::sort(pfoVector4.begin(), pfoVector4.end(), LArPfoHelper::SortByNHits);

        for (const ParticleFlowObject *const pDaughterPfo : pfoVector4)
        {
            const PfoAssociation &parentToDaughterAssociation(parentAssociationMap.at(pDaughterPfo));
            PfoAssociationMatrix::const_iterator iter5 = bestAssociationMatrix.find(pDaughterPfo);

            if (bestAssociationMatrix.end() == iter5)
                continue;

            const PfoAssociationMap &daughterAssociationMap(iter5->second);

            PfoAssociationMap::const_iterator iter6 = daughterAssociationMap.find(pParentPfo);
            if (daughterAssociationMap.end() == iter6)
                continue;

            const PfoAssociation &daughterToParentAssociation(iter6->second);

            if (parentToDaughterAssociation.GetParent() == daughterToParentAssociation.GetDaughter() &&
                parentToDaughterAssociation.GetDaughter() == daughterToParentAssociation.GetParent())
            {
                pfoMatches[pParentPfo].push_back(pDaughterPfo);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SelectPfoMerges(const PfoMergeMap &pfoMatches, PfoMergeMap &pfoMerges) const
{
    PfoSet vetoSet;

    PfoVector inputPfoVector;
    for (const auto &mapEntry : pfoMatches)
        inputPfoVector.push_back(mapEntry.first);
    std::sort(inputPfoVector.begin(), inputPfoVector.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pInputPfo : inputPfoVector)
    {
        const PfoList &pfoList(pfoMatches.at(pInputPfo));

        for (const ParticleFlowObject *const pSeedPfo : pfoList)
        {
            if (vetoSet.count(pSeedPfo))
                continue;

            PfoList mergeList;
            this->CollectAssociatedPfos(pSeedPfo, pSeedPfo, pfoMatches, vetoSet, mergeList);

            vetoSet.insert(pSeedPfo);
            PfoList &selectedPfoList(pfoMerges[pSeedPfo]);
            selectedPfoList.push_back(pSeedPfo);

            for (const ParticleFlowObject *const pAssociatedPfo : mergeList)
            {
                // Check if particle has already been counted
                if (vetoSet.count(pAssociatedPfo) || (selectedPfoList.end() != std::find(selectedPfoList.begin(), selectedPfoList.end(), pAssociatedPfo)))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                vetoSet.insert(pAssociatedPfo);
                selectedPfoList.push_back(pAssociatedPfo);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CollectAssociatedPfos(const ParticleFlowObject *const pSeedPfo,
    const ParticleFlowObject *const pCurrentPfo, const PfoMergeMap &pfoMergeMap, const PfoSet &vetoSet, PfoList &associatedList) const
{
    if (vetoSet.count(pCurrentPfo))
        return;

    PfoMergeMap::const_iterator iter1 = pfoMergeMap.find(pCurrentPfo);

    if (pfoMergeMap.end() == iter1)
        return;

    for (PfoList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        const ParticleFlowObject *const pAssociatedPfo = *iter2;

        if (pAssociatedPfo == pSeedPfo)
            continue;

        if (associatedList.end() != std::find(associatedList.begin(), associatedList.end(), pAssociatedPfo))
            continue;

        associatedList.push_back(pAssociatedPfo);

        this->CollectAssociatedPfos(pSeedPfo, pAssociatedPfo, pfoMergeMap, vetoSet, associatedList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::OrderPfoMerges(const PfoToLArTPCMap &pfoToLArTPCMap, const ThreeDPointingClusterMap &pointingClusterMap,
    const PfoMergeMap &inputPfoMerges, PfoMergeMap &outputPfoMerges) const
{
    PfoVector inputPfoVector;
    for (const auto &mapEntry : inputPfoMerges)
        inputPfoVector.push_back(mapEntry.first);
    std::sort(inputPfoVector.begin(), inputPfoVector.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pInputPfo : inputPfoVector)
    {
        const PfoList &pfoList(inputPfoMerges.at(pInputPfo));

        float bestLength(0.f);
        const ParticleFlowObject *pVertexPfo = nullptr;

        for (PfoList::const_iterator iter1 = pfoList.begin(), iterEnd = pfoList.end(); iter1 != iterEnd; ++iter1)
        {
            const ParticleFlowObject *const pPfo1(*iter1);
            PfoToLArTPCMap::const_iterator tpcIter1 = pfoToLArTPCMap.find(pPfo1);
            ThreeDPointingClusterMap::const_iterator pointingIter1 = pointingClusterMap.find(pPfo1);

            if (pfoToLArTPCMap.end() == tpcIter1 || pointingClusterMap.end() == pointingIter1)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const LArTPC *const pLArTPC1(tpcIter1->second);
            const LArPointingCluster &pointingCluster1(pointingIter1->second);

            for (PfoList::const_iterator iter2 = iter1; iter2 != iterEnd; ++iter2)
            {
                const ParticleFlowObject *const pPfo2(*iter2);
                PfoToLArTPCMap::const_iterator tpcIter2 = pfoToLArTPCMap.find(pPfo2);
                ThreeDPointingClusterMap::const_iterator pointingIter2 = pointingClusterMap.find(pPfo2);

                if (pfoToLArTPCMap.end() == tpcIter2 || pointingClusterMap.end() == pointingIter2)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const LArTPC *const pLArTPC2(tpcIter2->second);
                const LArPointingCluster &pointingCluster2(pointingIter2->second);

                if (pLArTPC1 == pLArTPC2)
                    continue;

                const float thisLength(LArStitchingHelper::GetTPCDisplacement(*pLArTPC1, *pLArTPC2));

                if (thisLength < bestLength)
                    continue;

                bestLength = thisLength;

                try
                {
                    pVertexPfo = nullptr;

                    LArPointingCluster::Vertex nearVertex1, nearVertex2;
                    LArStitchingHelper::GetClosestVertices(*pLArTPC1, *pLArTPC2, pointingCluster1, pointingCluster2, nearVertex1, nearVertex2);

                    const LArPointingCluster::Vertex &farVertex1(
                        nearVertex1.IsInnerVertex() ? pointingCluster1.GetOuterVertex() : pointingCluster1.GetInnerVertex());
                    const LArPointingCluster::Vertex &farVertex2(
                        nearVertex2.IsInnerVertex() ? pointingCluster2.GetOuterVertex() : pointingCluster2.GetInnerVertex());
                    const float deltaY(farVertex1.GetPosition().GetY() - farVertex2.GetPosition().GetY());

                    if (std::fabs(deltaY) < std::numeric_limits<float>::epsilon())
                        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                    pVertexPfo = ((deltaY > 0.f) ? pPfo1 : pPfo2);
                }
                catch (const pandora::StatusCodeException &)
                {
                }
            }
        }

        if (pVertexPfo)
            outputPfoMerges[pVertexPfo].insert(outputPfoMerges[pVertexPfo].begin(), pfoList.begin(), pfoList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::StitchPfos(const MasterAlgorithm *const pAlgorithm, const ThreeDPointingClusterMap &pointingClusterMap,
    const PfoMergeMap &pfoMerges, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map, const std::vector<PfoVector>& pfosGroupedAlongZ) const
{
    std::vector<bool> isZgroupStitched(pfosGroupedAlongZ.size(), false);

    PfoVector pfoVectorToEnlarge;
    for (const auto &mapEntry : pfoMerges)
        pfoVectorToEnlarge.push_back(mapEntry.first);
    std::sort(pfoVectorToEnlarge.begin(), pfoVectorToEnlarge.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfoToEnlarge : pfoVectorToEnlarge)
    {

        const PfoList &pfoList(pfoMerges.at(pPfoToEnlarge));
        const PfoVector pfoVector(pfoList.begin(), pfoList.end());
        PfoToPointingVertexMatrix pfoToPointingVertexMatrix;
        float x0(0.f);

        if (!m_useXcoordinate || m_alwaysApplyT0Calculation)
        {
            try
            {
                // If stitching contributions are inconsistent, abort
                if (!this->CalculateX0(pfoToLArTPCMap, pointingClusterMap, pfoVector, x0, pfoToPointingVertexMatrix))
                    continue;
            }
            catch (const pandora::StatusCodeException &)
            {
                continue;
            }
        }

        // ATTN: shift the pfos one at a time
        PfoSet shiftedPfos;
        for (PfoVector::const_iterator iterI = pfoVector.begin(); iterI != pfoVector.end(); ++iterI)
        {
            const ParticleFlowObject *const pPfoI(*iterI);
            const LArTPC *const pLArTPCI(pfoToLArTPCMap.at(pPfoI));

            // std::cout << "\n";
            // std::cout << "pPfoI : " << pPfoI << "\n";
            for (PfoVector::const_iterator iterJ = std::next(iterI); iterJ != pfoVector.end(); ++iterJ)
            {
                const ParticleFlowObject *const pPfoJ(*iterJ);
                const LArTPC *const pLArTPCJ(pfoToLArTPCMap.at(pPfoJ));

                if (!LArStitchingHelper::CanTPCsBeStitched(*pLArTPCI, *pLArTPCJ))
                    continue;

                // std::cout << "pPfoJ : " << pPfoJ << " can be stitched\n";

                if (std::find(shiftedPfos.begin(), shiftedPfos.end(), pPfoI) == shiftedPfos.end())
                { // shift pPfoI
                    const float signedX0 = this->GetSignedX0(pPfoI, pPfoJ, x0, pfoToLArTPCMap, pfoToPointingVertexMatrix);

                    if (!m_useXcoordinate || m_alwaysApplyT0Calculation)
                        this->ShiftPfo(pAlgorithm, pPfoI, signedX0, pfoToLArTPCMap);

                    // std::cout << pPfoI  << " shifted by X0 " << signedX0 << "\n";
                    shiftedPfos.insert(pPfoI);
                    
                    if (m_stitchPfosZMatches && pfosGroupedAlongZ.size()!=0)
                    { // z stitching
                        const PfoVector* zMatches = nullptr;

                        // find pfoI z matches
                        // std::cout << "checking if pfo "<< pPfoI << " has a group of Z matches in pfosGroupedAlongZ\n";
                        if (LArPfoHelper::FindPfoGroup(pPfoI, pfosGroupedAlongZ, zMatches))
                        {
                          size_t index = std::distance(&pfosGroupedAlongZ[0], zMatches);
                          isZgroupStitched[index] = true;

                          for (const ParticleFlowObject* const zMatch : *zMatches)
                          { // loop on z matches
                            if (zMatch == pPfoI) continue; // already shifted
                            this->ShiftPfo(pAlgorithm, zMatch, signedX0, pfoToLArTPCMap);
                            shiftedPfos.insert(zMatch);
                            // std::cout << "found zMatch " << zMatch << " to shift\n";
                          }
                        }
                    } // z stitching 
                } // shift pPfoI

                if (std::find(shiftedPfos.begin(), shiftedPfos.end(), pPfoJ) == shiftedPfos.end())
                { // shift pPfoJ
                    const float signedX0 = this->GetSignedX0(pPfoJ, pPfoI, x0, pfoToLArTPCMap, pfoToPointingVertexMatrix);

                    if (!m_useXcoordinate || m_alwaysApplyT0Calculation)
                        this->ShiftPfo(pAlgorithm, pPfoJ, signedX0, pfoToLArTPCMap);

                    std::cout << pPfoJ  << " shifted by X0 " << signedX0 << "\n";
                    shiftedPfos.insert(pPfoJ);
                    
                    if (m_stitchPfosZMatches && pfosGroupedAlongZ.size()!=0)
                    { // z stitching
                        const PfoVector* zMatches = nullptr;

                        // find pfoJ z matches
                        // std::cout << "checking if pfo "<< pPfoJ << " has a group of Z matches in pfosGroupedAlongZ\n";
                        if (LArPfoHelper::FindPfoGroup(pPfoJ, pfosGroupedAlongZ, zMatches))
                        {
                          size_t index = std::distance(&pfosGroupedAlongZ[0], zMatches);
                          isZgroupStitched[index] = true;

                          for (const ParticleFlowObject* const zMatch : *zMatches)
                          { // loop on z matches
                            if (zMatch == pPfoJ) continue; // already shifted
                            this->ShiftPfo(pAlgorithm, zMatch, signedX0, pfoToLArTPCMap);
                            shiftedPfos.insert(zMatch);
                            // std::cout << "found zMatch " << zMatch << "to shift by " << signedX0 << "\n";
                          }
                        }
                    } // z stitching pPfoJ
                } // shift pPfoJ
            }
        }

        // now merge all pfos
        for (const ParticleFlowObject *const pPfoToDelete : shiftedPfos)
        {
            if (pPfoToDelete == pPfoToEnlarge)
                continue;

            pAlgorithm->StitchPfos(pPfoToEnlarge, pPfoToDelete, pfoToLArTPCMap);
            std::cout << pPfoToDelete << " stitched to " << pPfoToEnlarge << "\n";
        }

        stitchedPfosToX0Map.insert(PfoToFloatMap::value_type(pPfoToEnlarge, x0));
    }

    // stitch across z remaining pfos
    std::cout << "\n";
    std::cout << "stitch remaining groups of pfos along z if any\n";
    if (m_stitchPfosZMatches)
    {
      for (size_t i = 0; i < pfosGroupedAlongZ.size(); ++i )
      {
        if (isZgroupStitched[i] == false && pfosGroupedAlongZ[i].size() > 1)
        {
          std::cout << "found a group to stitch \n";
          const ParticleFlowObject* const pPfoToEnlarge = pfosGroupedAlongZ[i][0];

          std::cout << "pfo to enlarge " << pPfoToEnlarge << "\n";

          if(pfosGroupedAlongZ[i].size() == 1) continue; // extra coutious check
                                                         //
          for (size_t j = 1; j < pfosGroupedAlongZ[i].size(); ++j)
          {
            const ParticleFlowObject* pPfoToDelete = pfosGroupedAlongZ[i][j];

            pAlgorithm->StitchPfos(pPfoToEnlarge, pPfoToDelete, pfoToLArTPCMap);

            std::cout << "pfo " << pPfoToDelete << " deleted\n";
          }
          isZgroupStitched[i] = true;
        }
      }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float StitchingCosmicRayMergingTool::GetSignedX0(const ParticleFlowObject *const pPfoToShift,
    const ParticleFlowObject *const pMatchedPfo, const float x0, const PfoToLArTPCMap &pfoToLArTPCMap,
    const PfoToPointingVertexMatrix &pfoToPointingVertexMatrix) const
{
  // get stitching vertex for the pfo to be shifted
    const PfoToPointingVertexMatrix::const_iterator pfoToPointingVertexMatrixIter(pfoToPointingVertexMatrix.find(pPfoToShift));
    const LArPointingCluster::Vertex stitchingVertex(pfoToPointingVertexMatrixIter->second.at(pMatchedPfo));

    const LArTPC *const pShiftLArTPC(pfoToLArTPCMap.at(pPfoToShift));
    const LArTPC *const pMatchedLArTPC(pfoToLArTPCMap.at(pMatchedPfo));

    // determine shift sign from the relative position of stitching vertex and the relevant TPC boundary position
    const float tpcBoundaryCenterX(LArStitchingHelper::GetTPCBoundaryCenterX(*pShiftLArTPC, *pMatchedLArTPC));
    float tpcBoundaryX(0.f);

    if (pShiftLArTPC->GetCenterX() < tpcBoundaryCenterX)
    {
        tpcBoundaryX = pShiftLArTPC->GetCenterX() + (0.5f * pShiftLArTPC->GetWidthX());
    }
    else
    {
        tpcBoundaryX = pShiftLArTPC->GetCenterX() - (0.5f * pShiftLArTPC->GetWidthX());
    }
    // std:: cout << "test\n";
    // std::cout<< LArStitchingHelper::GetTPCBoundaryCenterX(*pShiftLArTPC, *pMatchedLArTPC) << "\n";
    // std::cout<< LArStitchingHelper::GetTPCBoundaryCenterX(*pMatchedLArTPC, *pShiftLArTPC) << "\n";
    // const float test1 = LArStitchingHelper::GetTPCBoundaryCenterX(*pShiftLArTPC, *pMatchedLArTPC); 
    // const float test2 = LArStitchingHelper::GetTPCBoundaryCenterX(*pMatchedLArTPC, *pShiftLArTPC); 
    // std::cout << test1 <<"\n";
    // std::cout << test2 <<"\n";
    // std:: cout << "test\n";
    //
    const float signX0 = stitchingVertex.GetPosition().GetX() < tpcBoundaryX ? 1.f : -1.f;
    // std::cout << "\n";
    // std::cout << "tpcBoundaryX : " << tpcBoundaryX 
    //           << ", TPC with pfo to shift : " << pShiftLArTPC->GetLArTPCVolumeId() << ", with X center : " <<pShiftLArTPC->GetCenterX()
    //           << ", other : " << pMatchedLArTPC->GetLArTPCVolumeId() << ", with X center : " <<pMatchedLArTPC->GetCenterX()
    //           << ", TPC boundaryX " << tpcBoundaryX
    //           << ", tpcBoundaryCenterX " << tpcBoundaryCenterX 
    //           << ", vertexX " << stitchingVertex.GetPosition().GetX()
    //           << ", signedX0 : "<< signX0
    //           <<"\n";
    //
    return std::fabs(x0) * signX0; 
}
//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::ShiftPfo(const MasterAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfoToShift, const float signedX0, const PfoToLArTPCMap &pfoToLArTPCMap) const
{
    // ATTN: No CPA/APA sign needed since x0 calculation corresponds to an APA
    object_creation::ParticleFlowObject::Metadata metadata;
    metadata.m_propertiesToAdd["X0"] = std::fabs(signedX0);

    // ATTN: Set the X0 shift for all particles in hierarchy
    PfoList downstreamPfoList;
    LArPfoHelper::GetAllDownstreamPfos(pPfoToShift, downstreamPfoList);

    for (const ParticleFlowObject *const pHierarchyPfo : downstreamPfoList)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pHierarchyPfo, metadata));

    pAlgorithm->ShiftPfoHierarchy(pPfoToShift, pfoToLArTPCMap, signedX0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool StitchingCosmicRayMergingTool::CalculateX0(const PfoToLArTPCMap &pfoToLArTPCMap, const ThreeDPointingClusterMap &pointingClusterMap,
    const PfoVector &pfoVector, float &x0, PfoToPointingVertexMatrix &pfoToPointingVertexMatrix) const
{
    float sumX(0.f), sumN(0.f);

    for (PfoVector::const_iterator iter1 = pfoVector.begin(), iterEnd = pfoVector.end(); iter1 != iterEnd; ++iter1)
    {
        const ParticleFlowObject *const pPfo1(*iter1);
        PfoToLArTPCMap::const_iterator tpcIter1 = pfoToLArTPCMap.find(pPfo1);
        ThreeDPointingClusterMap::const_iterator pointingIter1 = pointingClusterMap.find(pPfo1);

        if (pfoToLArTPCMap.end() == tpcIter1 || pointingClusterMap.end() == pointingIter1)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const LArTPC *const pLArTPC1(tpcIter1->second);
        const LArPointingCluster &pointingCluster1(pointingIter1->second);

        for (PfoVector::const_iterator iter2 = std::next(iter1); iter2 != iterEnd; ++iter2)
        {
            const ParticleFlowObject *const pPfo2(*iter2);
            PfoToLArTPCMap::const_iterator tpcIter2 = pfoToLArTPCMap.find(pPfo2);
            ThreeDPointingClusterMap::const_iterator pointingIter2 = pointingClusterMap.find(pPfo2);

            if (pfoToLArTPCMap.end() == tpcIter2 || pointingClusterMap.end() == pointingIter2)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const LArTPC *const pLArTPC2(tpcIter2->second);
            const LArPointingCluster &pointingCluster2(pointingIter2->second);

            if (!LArStitchingHelper::CanTPCsBeStitched(*pLArTPC1, *pLArTPC2))
                continue;

            // Calculate X0 for the closest pair of vertices
            LArPointingCluster::Vertex pointingVertex1, pointingVertex2;
            try
            {
                LArStitchingHelper::GetClosestVertices(*pLArTPC1, *pLArTPC2, pointingCluster1, pointingCluster2, pointingVertex1, pointingVertex2);

                // Record pfo1 stitching vertex for pfo1<->pfo2 match, used to determine shifting direction in later step
                const PfoToPointingVertexMatrix::iterator pfoToPointingVertexMatrixIter1(pfoToPointingVertexMatrix.find(pPfo1));
                if (pfoToPointingVertexMatrixIter1 == pfoToPointingVertexMatrix.end())
                {
                    // No matches present in map, add pfo1<->pfo2 match
                    const PfoToPointingVertexMap pfoToPointingVertexMap({{pPfo2, pointingVertex1}});
                    (void)pfoToPointingVertexMatrix.insert(PfoToPointingVertexMatrix::value_type(pPfo1, pfoToPointingVertexMap));
                }
                else
                {
                    // ATTN: another match for a different TPC boundary may be present, add pfo1<->pfo2 match
                    PfoToPointingVertexMap &pfoToPointingVertexMap(pfoToPointingVertexMatrixIter1->second);
                    const PfoToPointingVertexMap::iterator pfoToPointingVertexMapIter(pfoToPointingVertexMap.find(pPfo2));
                    if (pfoToPointingVertexMapIter == pfoToPointingVertexMap.end())
                    {
                        (void)pfoToPointingVertexMap.insert(PfoToPointingVertexMap::value_type(pPfo2, pointingVertex1));
                    }
                    else
                    {
                        if ((pfoToPointingVertexMapIter->second.GetPosition() - pointingVertex1.GetPosition()).GetMagnitude() >
                            std::numeric_limits<float>::epsilon())
                            throw StatusCodeException(STATUS_CODE_FAILURE);
                    }
                }

                // Record pfo2 stitching vertex for pfo1<->pfo2 match, used to determine shifting direction in later step
                const PfoToPointingVertexMatrix::iterator pfoToPointingVertexMatrixIter2(pfoToPointingVertexMatrix.find(pPfo2));
                if (pfoToPointingVertexMatrixIter2 == pfoToPointingVertexMatrix.end())
                {
                    // No matches present in map, add pfo1<->pfo2 match
                    const PfoToPointingVertexMap pfoToPointingVertexMap({{pPfo1, pointingVertex2}});
                    (void)pfoToPointingVertexMatrix.insert(PfoToPointingVertexMatrix::value_type(pPfo2, pfoToPointingVertexMap));
                }
                else
                {
                    // ATTN: another match for a different TPC boundary may be present, add pfo1<->pfo2 match
                    PfoToPointingVertexMap &pfoToPointingVertexMap(pfoToPointingVertexMatrixIter2->second);
                    const PfoToPointingVertexMap::iterator pfoToPointingVertexMapIter(pfoToPointingVertexMap.find(pPfo1));
                    if (pfoToPointingVertexMapIter == pfoToPointingVertexMap.end())
                    {
                        (void)pfoToPointingVertexMap.insert(PfoToPointingVertexMap::value_type(pPfo1, pointingVertex2));
                    }
                    else
                    {
                        if ((pfoToPointingVertexMapIter->second.GetPosition() - pointingVertex2.GetPosition()).GetMagnitude() >
                            std::numeric_limits<float>::epsilon())
                            throw StatusCodeException(STATUS_CODE_FAILURE);
                    }
                }

                const float tpcBoundaryCenterX(LArStitchingHelper::GetTPCBoundaryCenterX(*pLArTPC1, *pLArTPC2));
                const bool isCPAStitch(pLArTPC1->GetCenterX() < tpcBoundaryCenterX ? !pLArTPC1->IsDriftInPositiveX() : !pLArTPC2->IsDriftInPositiveX());
                float thisX0(LArStitchingHelper::CalculateX0(*pLArTPC1, *pLArTPC2, pointingVertex1, pointingVertex2));

                thisX0 *= isCPAStitch ? -1.f : 1.f;

                // If multiple boundaries identified, check if stitching contribution is consistent
                if ((sumN > std::numeric_limits<float>::epsilon()) && (sumX > std::numeric_limits<float>::epsilon()))
                {
                    const float fractionalDiff(std::fabs((sumX - (thisX0 * sumN)) / sumX));

                    if ((fractionalDiff > m_maxX0FractionalDeviation) && (std::fabs(sumX / sumN) > m_boundaryToleranceWidth))
                        return false;
                }

                sumX += thisX0;
                sumN += 1.f;
            }
            catch (const pandora::StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    std::cout << "StitchingCosmicRayMergingTool: Attempting to stitch a pfo with multiple vertices for the same match" << std::endl;
            }
        }
    }

    if ((sumN < std::numeric_limits<float>::epsilon()) || (std::fabs(sumX) < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    x0 = (sumX / sumN);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StitchingCosmicRayMergingTool::PfoAssociation::PfoAssociation(const VertexType parent, const VertexType daughter, const float fom) :
    m_parent(parent),
    m_daughter(daughter),
    m_fom(fom)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StitchingCosmicRayMergingTool::PfoAssociation::VertexType StitchingCosmicRayMergingTool::PfoAssociation::GetParent() const
{
    return m_parent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StitchingCosmicRayMergingTool::PfoAssociation::VertexType StitchingCosmicRayMergingTool::PfoAssociation::GetDaughter() const
{
    return m_daughter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float StitchingCosmicRayMergingTool::PfoAssociation::GetFigureOfMerit() const
{
    return m_fom;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingCosmicRayMergingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ThreeDStitchingMode", m_useXcoordinate));

PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ApplyZStitching", m_stitchPfosZMatches));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "AlwaysApplyT0Calculation", m_alwaysApplyT0Calculation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HalfWindowLayers", m_halfWindowLayers));

    float minLength(std::sqrt(m_minLengthSquared));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPfoLength", minLength));
    m_minLengthSquared = minLength * minLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCosRelativeAngle", m_minCosRelativeAngle));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxXYdisplacementPfoZmatching", m_maxXYdisplacementPfoZmatching));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxZdistancePfoZmatching", m_maxZdistancePfoZmatching));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "RelaxMinLongitudinalDisplacement", m_relaxMinLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxLongitudinalDisplacementX", m_maxLongitudinalDisplacementX));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "LooserMinCosRelativeAngle", m_relaxCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "LooserMaxTransverseDisplacement", m_relaxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinNCaloHits3D", m_minNCaloHits3D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxX0FractionalDeviation", m_maxX0FractionalDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BoundaryToleranceWidth", m_boundaryToleranceWidth));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
