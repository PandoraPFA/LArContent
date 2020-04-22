/**
 *  @file   LArContent/src/LArControlFlow/StitchingCosmicRayMergingTool.cc
 *
 *  @brief  Implementation of the stitching cosmic ray merging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

#include "larpandoracontent/LArControlFlow/StitchingCosmicRayMergingTool.h"

using namespace pandora;

namespace lar_content
{

StitchingCosmicRayMergingTool::StitchingCosmicRayMergingTool() :
    m_useXcoordinate(false),
    m_alwaysApplyT0Calculation(true),
    m_halfWindowLayers(30),
    m_minLengthSquared(50.f),
    m_minCosRelativeAngle(0.966),
    m_maxLongitudinalDisplacementX(15.f),
    m_maxTransverseDisplacement(5.f),
    m_relaxCosRelativeAngle(0.906),
    m_relaxTransverseDisplacement(2.5f),
    m_minNCaloHits3D(0),
    m_writeToTree(false),
    m_fileName("ImpactParameters.root"),
    m_treeName("ImpactParameterTree"),
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StitchingCosmicRayMergingTool::~StitchingCosmicRayMergingTool() 
{
    if(m_writeToTree) {
        try {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
            std::cout << "TREE SAVED" << std::endl;
        }
        catch (const StatusCodeException &) {
            std::cout << "UNABLE TO WRITE TREE" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::Run(const MasterAlgorithm *const pAlgorithm, const PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map)
{
    
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    m_eventNumber++;
    
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
    
    PfoAssociationMatrix pfoAssociationMatrix;
    this->CreatePfoMatches(larTPCToPfoMap, pointingClusterMap, pfoAssociationMatrix);
    
    PfoMergeMap pfoSelectedMatches;
    this->SelectPfoMatches(pfoAssociationMatrix, pfoSelectedMatches);
    
    PfoMergeMap pfoSelectedMerges;
    this->SelectPfoMerges(pfoSelectedMatches, pfoSelectedMerges);
    
    PfoMergeMap pfoOrderedMerges;
    this->OrderPfoMerges(pfoToLArTPCMap, pointingClusterMap, pfoSelectedMerges, pfoOrderedMerges);

    this->StitchPfos(pAlgorithm, pointingClusterMap, pfoOrderedMerges, pfoToLArTPCMap, stitchedPfosToX0Map);

    std::cout << "HERE" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SelectPrimaryPfos(const PfoList *pInputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, PfoList &outputPfoList) const
{
    for (const ParticleFlowObject *const pPfo : *pInputPfoList)
    {
        if (!LArPfoHelper::IsFinalState(pPfo) || !LArPfoHelper::IsTrack(pPfo))
            continue;

        if (!pfoToLArTPCMap.count(pPfo))
            continue;

        outputPfoList.push_back(pPfo);
    }

    outputPfoList.sort(LArPfoHelper::SortByNHits);

    /*
    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &outputPfoList, "TO MATCH PFOS", RED);

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::BuildPointingClusterMaps(const PfoList &inputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, ThreeDPointingClusterMap &pointingClusterMap) const
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
            (void) pointingClusterMap.insert(ThreeDPointingClusterMap::value_type(pPfo, LArPointingCluster(slidingFitResult)));
        }
        catch (const StatusCodeException &) {}
    }
    /*
    for (const ParticleFlowObject *const pPfo : inputPfoList)
    {
        const auto iter(pfoToLArTPCMap.find(pPfo));
        if (iter == pfoToLArTPCMap.end())
            continue;
        PfoList aList;
        aList.push_back(pPfo);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &aList, "CLUSTERED", RED);
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
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

void StitchingCosmicRayMergingTool::CreatePfoMatches(const LArTPCToPfoMap &larTPCToPfoMap, const ThreeDPointingClusterMap &pointingClusterMap,
    PfoAssociationMatrix &pfoAssociationMatrix) const
{
    LArTPCVector larTPCVector;
    for (const auto &mapEntry : larTPCToPfoMap) larTPCVector.push_back(mapEntry.first);
    std::sort(larTPCVector.begin(), larTPCVector.end(), LArStitchingHelper::SortTPCs);

    /*
    for (const auto &lartpc : larTPCVector)
    {
        const CartesianVector centre(lartpc->GetCenterX(), lartpc->GetCenterY(), lartpc->GetCenterZ());
        const CartesianVector origin(0,0,0);

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &centre, "CENTRE", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &origin, "CENTRE", RED, 2);
        PandoraMonitoringApi::Pause(this->GetPandora());
    }
    
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
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
                {
                    //////////
                    /*
                    PfoList pfo1List, pfo2List;
                    pfo1List.push_back(pPfo1);
                    pfo2List.push_back(pPfo2);
                    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfo1List, "PFO1", RED);
                    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfo2List, "PFO2", VIOLET);
                    */
                    //////////
                    
                    this->CreatePfoMatches(*pLArTPC1, *pLArTPC2, pPfo1, pPfo2, pointingClusterMap, pfoAssociationMatrix);

                    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                }
            }
        }
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CreatePfoMatches(const LArTPC &larTPC1, const LArTPC &larTPC2,
    const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2,
    const ThreeDPointingClusterMap &pointingClusterMap, PfoAssociationMatrix &pfoAssociationMatrix) const
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

    //std::cout << "A" << std::endl;
    
    const LArPointingCluster &pointingCluster1(iter1->second);
    const LArPointingCluster &pointingCluster2(iter2->second);

    // Check length of pointing clusters
    if (pointingCluster1.GetLengthSquared() < m_minLengthSquared || pointingCluster2.GetLengthSquared() < m_minLengthSquared)
        return;

    //std::cout << "B" << std::endl;



    // Get number of 3D hits in each of the pfos
    CaloHitList caloHitList3D1;
    LArPfoHelper::GetCaloHits(pPfo1, TPC_3D, caloHitList3D1);

    CaloHitList caloHitList3D2;
    LArPfoHelper::GetCaloHits(pPfo2, TPC_3D, caloHitList3D2);

    // Check number of 3D hits in each of the pfos
    if (caloHitList3D1.size() < m_minNCaloHits3D || caloHitList3D2.size() < m_minNCaloHits3D)
        return;

    //std::cout << "C" << std::endl;

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


    //////////
    /*
    CartesianVector position1(pointingVertex1.GetPosition());
    CartesianVector position2(pointingVertex2.GetPosition());
    
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position1, "VERTEX 1", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position2, "VERTEX 2", VIOLET, 2);
    */
    //////////
    
    //std::cout << "D" << std::endl;
    

    // Pointing clusters must have a parallel direction
    const float cosRelativeAngle(-pointingVertex1.GetDirection().GetDotProduct(pointingVertex2.GetDirection()));

    //std::cout << "COS RELATIVE ANGLE: " << cosRelativeAngle << std::endl;
    
    if (cosRelativeAngle < m_relaxCosRelativeAngle)
        return;

    //std::cout << "E" << std::endl;

    // Pointing clusters must have a non-zero X direction (so that they point across drift volume boundary)
    const float pX1(std::fabs(pointingVertex1.GetDirection().GetX()));
    const float pX2(std::fabs(pointingVertex2.GetDirection().GetX()));

    if (pX1 < std::numeric_limits<float>::epsilon() || pX2 < std::numeric_limits<float>::epsilon())
        return;

    //std::cout << "F" << std::endl;
    
    // Pointing clusters must intersect at a drift volume boundary
    const float intersectX(0.5 * (pointingVertex1.GetPosition().GetX() + pointingVertex2.GetPosition().GetX()));

    if (std::fabs(intersectX - boundaryCenterX) > maxLongitudinalDisplacementX)
        return;

    //std::cout << "G" << std::endl;

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

    //std::cout << "H" << std::endl;

    // Selection cuts on longitudinal impact parameters
    const float minL(-1.f);
    const float dXdL1(m_useXcoordinate ? pX1 :
        (1.f - pX1 * pX1 > std::numeric_limits<float>::epsilon()) ? pX1 / std::sqrt(1.f - pX1 * pX1) : minL);
    const float dXdL2(m_useXcoordinate ? pX2 :
        (1.f - pX2 * pX2 > std::numeric_limits<float>::epsilon()) ? pX2 / std::sqrt(1.f - pX2 * pX2) : minL);
    const float maxL1(maxLongitudinalDisplacementX / dXdL1);
    const float maxL2(maxLongitudinalDisplacementX / dXdL2);

    /*
    std::cout << "minL: " << minL << std::endl;
    std::cout << "maxL1: " << maxL1 << std::endl;
    std::cout << "maxL2: " << maxL2 << std::endl;
    std::cout << "rL1: " << rL1 << std::endl;
    std::cout << "rL2: " << rL2 << std::endl;

    std::cout << "rT1" << rT1 << std::endl;
    std::cout << "rT2" << rT2 << std::endl;
    */


    if (rL1 < minL || rL1 > maxL1 || rL2 < minL || rL2 > maxL2)
        return;
    
    //std::cout << "I" << std::endl;
    
    // Selection cuts on transverse impact parameters
    const bool minPass(std::min(rT1, rT2) < m_relaxTransverseDisplacement && cosRelativeAngle > m_relaxCosRelativeAngle);
    const bool maxPass(std::max(rT1, rT2) < m_maxTransverseDisplacement && cosRelativeAngle > m_minCosRelativeAngle);

    if (!minPass && !maxPass)
        return;

    //std::cout << "J" << std::endl;

    // Store this association
    const PfoAssociation::VertexType vertexType1(pointingVertex1.IsInnerVertex() ? PfoAssociation::INNER : PfoAssociation::OUTER);
    const PfoAssociation::VertexType vertexType2(pointingVertex2.IsInnerVertex() ? PfoAssociation::INNER : PfoAssociation::OUTER);

    const float particleLength1(pointingCluster1.GetLengthSquared());
    const float particleLength2(pointingCluster2.GetLengthSquared());

    pfoAssociationMatrix[pPfo1].insert(PfoAssociationMap::value_type(pPfo2, PfoAssociation(vertexType1, vertexType2, particleLength2)));
    pfoAssociationMatrix[pPfo2].insert(PfoAssociationMap::value_type(pPfo1, PfoAssociation(vertexType2, vertexType1, particleLength1)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SelectPfoMatches(const PfoAssociationMatrix &pfoAssociationMatrix, PfoMergeMap &pfoMatches) const
{
    // First step: loop over association matrix and find best associations A -> X and B -> Y
    // =====================================================================================
    PfoAssociationMatrix bestAssociationMatrix;

    PfoVector pfoVector1;
    for (const auto &mapEntry : pfoAssociationMatrix) pfoVector1.push_back(mapEntry.first);
    std::sort(pfoVector1.begin(), pfoVector1.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo1 : pfoVector1)
    {
        const PfoAssociationMap &pfoAssociationMap(pfoAssociationMatrix.at(pPfo1));

        const ParticleFlowObject *pBestPfoInner(nullptr);
        PfoAssociation bestAssociationInner(PfoAssociation::UNDEFINED, PfoAssociation::UNDEFINED, 0.f);

        const ParticleFlowObject *pBestPfoOuter(nullptr);
        PfoAssociation bestAssociationOuter(PfoAssociation::UNDEFINED, PfoAssociation::UNDEFINED, 0.f);

        PfoVector pfoVector2;
        for (const auto &mapEntry : pfoAssociationMap) pfoVector2.push_back(mapEntry.first);
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
            (void) bestAssociationMatrix[pPfo1].insert(PfoAssociationMap::value_type(pBestPfoInner, bestAssociationInner));

        if (pBestPfoOuter)
            (void) bestAssociationMatrix[pPfo1].insert(PfoAssociationMap::value_type(pBestPfoOuter, bestAssociationOuter));
    }

    // Second step: make the merge if A -> X and B -> Y is in fact A -> B and B -> A
    // =============================================================================
    PfoVector pfoVector3;
    for (const auto &mapEntry : bestAssociationMatrix) pfoVector3.push_back(mapEntry.first);
    std::sort(pfoVector3.begin(), pfoVector3.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pParentPfo : pfoVector3)
    {
        const PfoAssociationMap &parentAssociationMap(bestAssociationMatrix.at(pParentPfo));

        PfoVector pfoVector4;
        for (const auto &mapEntry : parentAssociationMap) pfoVector4.push_back(mapEntry.first);
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
    for (const auto &mapEntry : pfoMatches) inputPfoVector.push_back(mapEntry.first);
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

void StitchingCosmicRayMergingTool::CollectAssociatedPfos(const ParticleFlowObject *const pSeedPfo, const ParticleFlowObject *const pCurrentPfo,
    const PfoMergeMap &pfoMergeMap, const PfoSet &vetoSet, PfoList &associatedList) const
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
    for (const auto &mapEntry : inputPfoMerges) inputPfoVector.push_back(mapEntry.first);
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

                    const LArPointingCluster::Vertex &farVertex1(nearVertex1.IsInnerVertex() ? pointingCluster1.GetOuterVertex() : pointingCluster1.GetInnerVertex());
                    const LArPointingCluster::Vertex &farVertex2(nearVertex2.IsInnerVertex() ? pointingCluster2.GetOuterVertex() : pointingCluster2.GetInnerVertex());
                    const float deltaY(farVertex1.GetPosition().GetY() - farVertex2.GetPosition().GetY());

                    if (std::fabs(deltaY) < std::numeric_limits<float>::epsilon())
                         throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                    pVertexPfo = ((deltaY > 0.f) ? pPfo1 : pPfo2);
                }
                catch (const pandora::StatusCodeException &) {}
            }
        }

        if (pVertexPfo)
            outputPfoMerges[pVertexPfo].insert(outputPfoMerges[pVertexPfo].begin(), pfoList.begin(), pfoList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::StitchPfos(const MasterAlgorithm *const pAlgorithm, const ThreeDPointingClusterMap &pointingClusterMap,
    const PfoMergeMap &pfoMerges, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map) const
{
    PfoVector pfoVectorToEnlarge;
    for (const auto &mapEntry : pfoMerges) pfoVectorToEnlarge.push_back(mapEntry.first);
    std::sort(pfoVectorToEnlarge.begin(), pfoVectorToEnlarge.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfoToEnlarge : pfoVectorToEnlarge)
    {
        const PfoList &pfoList(pfoMerges.at(pPfoToEnlarge));
        const PfoVector pfoVector(pfoList.begin(), pfoList.end());
        
        float x0(0.f);
        PfoToPointingVertexMatrix pfoToPointingVertexMatrix;
        if (!m_useXcoordinate || m_alwaysApplyT0Calculation)
        {
            try
            {
                this->CalculateX0(pfoToLArTPCMap, pointingClusterMap, pfoVector, x0, pfoToPointingVertexMatrix);
                //std::cout << "XO: " << x0 << std::endl;
            }
            catch (const pandora::StatusCodeException &)
            {
                continue;
            }
        }

        ///////////////
        /*
        PfoList thePfoToEnlarge;
        thePfoToEnlarge.push_back(pPfoToEnlarge);
        std::cout << "SIZE OF MATCHES VECTOR: " << pfoVector.size() << std::endl;
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &thePfoToEnlarge, "PFO TO ENLARGE", RED);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "MERGE PFOS BEFORE SHIFTED", BLUE);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        ///////////////

        if (m_writeToTree)
        {
            int matchesSize(pfoVector.size());
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "APA_X0", x0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Matches", matchesSize));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "EventNumber", m_eventNumber));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
        }
        

        // first shift the pfos
        PfoList shiftedPfos;
        for (PfoVector::const_iterator iterI = pfoVector.begin(); iterI != pfoVector.end(); ++iterI)
        {
            const ParticleFlowObject *const pPfoI(*iterI);
            const LArTPC *const pLArTPCI(pfoToLArTPCMap.at(pPfoI));

            for (PfoVector::const_iterator iterJ = iterI; iterJ != pfoVector.end(); ++iterJ)
            {
                if (iterI == iterJ)
                    continue;

                const ParticleFlowObject *const pPfoJ(*iterJ);
                const LArTPC *const pLArTPCJ(pfoToLArTPCMap.at(pPfoJ));

                if (!LArStitchingHelper::CanTPCsBeStitched(*pLArTPCI, *pLArTPCJ))
                    continue;

                PfoList pfosToShift({pPfoI, pPfoJ});

                ///////////////
                /*
                PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfosToShift, "MERGE PAIR", BLACK);
                std::cout << "XO FOR APA: " << x0 << std::endl;
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                */
                ///////////////
                
                for (const ParticleFlowObject *const pPfoToShift : pfosToShift)
                {
                    // if pfo has already been shifted
                    if (std::find(shiftedPfos.begin(), shiftedPfos.end(), pPfoToShift) != shiftedPfos.end())
                        continue;

                    if (!m_useXcoordinate || m_alwaysApplyT0Calculation)
                    {                    
                        const LArTPC *const pLArTPCShift(pfoToLArTPCMap.at(pPfoToShift));
                        const float tpcBoundaryCenterX(LArStitchingHelper::GetTPCBoundaryCenterX(*pLArTPCI, *pLArTPCJ));
                        float tpcBoundaryX(0.f);

                        if (pLArTPCShift->GetCenterX() < tpcBoundaryCenterX)
                        {
                            tpcBoundaryX = pLArTPCShift->GetCenterX() + (pLArTPCShift->GetWidthX()/2);
                        }
                        else
                        {
                            tpcBoundaryX = pLArTPCShift->GetCenterX() - (pLArTPCShift->GetWidthX()/2);
                        }

                        const PfoToPointingVertexMatrix::iterator pfoToPointingVertexMatrixIter(pfoToPointingVertexMatrix.find(pPfoToShift));
                        LArPointingCluster::Vertex stitchingVertex;
                        for (const ParticleFlowObject *const pJam : pfosToShift)
                        {
                            if (pJam == pPfoToShift)
                                continue;

                            stitchingVertex = pfoToPointingVertexMatrixIter->second.at(pJam);
                        }
                    
                    float positionShiftSign(0.f);
                    positionShiftSign = stitchingVertex.GetPosition().GetX() < tpcBoundaryX ? 1.f : -1.f;

                //ISOBEL: WORK OUT WHAT TO DO HERE
                //const float t0Sign(isCPAStitch ? -1.f : 1.f);
                //object_creation::ParticleFlowObject::Metadata metadata;
                //metadata.m_propertiesToAdd["X0"] = x0 * t0Sign;
            

                // ATTN: Set the X0 shift for all particles in hierarchy
                //PfoList downstreamPfoList;
                //LArPfoHelper::GetAllDownstreamPfos(pPfoToShift, downstreamPfoList);

                //for (const ParticleFlowObject *const pHierarchyPfo : downstreamPfoList)
                    //PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pHierarchyPfo, metadata));
                     
                    
                        const float signedX0(std::fabs(x0) * positionShiftSign);

                        ///////////////
                        /*
                        PfoList theParticles;
                        theParticles.push_back(pPfoToShift);
                        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &theParticles, "WILL SHIFT PFO", DARKGREEN);
                        std::cout << "SIGNED XO: " << signedX0 << std::endl;
                        PandoraMonitoringApi::Pause(this->GetPandora());
                        */
                        ///////////////
                        
                        pAlgorithm->ShiftPfoHierarchy(pPfoToShift, pfoToLArTPCMap, signedX0);
                    }
                    

                    shiftedPfos.push_back(pPfoToShift);
                }
                
            }
        }

        ///////////////
        //PfoList shiftList(shiftedPfos.begin(), shiftedPfos.end());
        //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &shiftList, "MERGE PFOS AFTER SHIFTED", RED);
        ///////////////

        // now merge
        for (const ParticleFlowObject *const pPfoToDelete : shiftedPfos)
        {
            if (pPfoToDelete == pPfoToEnlarge)
                continue;

            pAlgorithm->StitchPfos(pPfoToEnlarge, pPfoToDelete, pfoToLArTPCMap);
<<<<<<< HEAD
            stitchedPfosToX0Map.insert(PfoToFloatMap::value_type(pPfoToEnlarge, x0)); // i think this x0 is okay to ignore the sign of ISOBEL PLEASE CHECK
=======
>>>>>>> 25a640b2... X0 Distribtution Tree
        }

        stitchedPfosToX0Map.insert(PfoToFloatMap::value_type(pPfoToEnlarge, x0));

        ///////////////
        /*
        PfoList mergedList;
        mergedList.push_back(pPfoToEnlarge);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &mergedList, "MERGED PFOS", VIOLET);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        ///////////////
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *StitchingCosmicRayMergingTool::ReduceToLongestStitch(const PfoVector &pfoVector, const ParticleFlowObject *const pPfoToEnlarge,
    const PfoToLArTPCMap &pfoToLArTPCMap, PfoVector &reducedPfoVector) const
{
    if (pfoVector.size() > 2)
    {
        this->SelectLongestStitch(pfoVector, pfoToLArTPCMap, reducedPfoVector);
        std::sort(reducedPfoVector.begin(), reducedPfoVector.end(), LArPfoHelper::SortByNHits);

        return this->GetClosestPfo(pPfoToEnlarge, reducedPfoVector);
    }
    else
    {
        reducedPfoVector = pfoVector;

        return pPfoToEnlarge;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SelectLongestStitch(const PfoVector &pfoVector, const PfoToLArTPCMap &pfoToLArTPCMap, PfoVector &reducedPfoVector) const
{
    //ATTN this is the initialization, in case the rest of logic fails
    reducedPfoVector = pfoVector;
    reducedPfoVector.resize(2);

    unsigned int totalHitsStitched(0);

    for (PfoVector::const_iterator iterPfo1 = pfoVector.begin(), iterPfoEnd = pfoVector.end(); iterPfo1 != iterPfoEnd; ++iterPfo1)
    {
        for (PfoVector::const_iterator iterPfo2 =  iterPfo1; iterPfo2 != iterPfoEnd; ++iterPfo2)
        {
            const ParticleFlowObject *const pPfoToShift1(*iterPfo1);
            const ParticleFlowObject *const pPfoToShift2(*iterPfo2);

            if (pPfoToShift2 == pPfoToShift1)
                continue;

            const unsigned int twoDHitsPfo1(LArPfoHelper::GetNumberOfTwoDHits(pPfoToShift1)), twoDHitsPfo2(LArPfoHelper::GetNumberOfTwoDHits(pPfoToShift2));

            PfoToLArTPCMap::const_iterator iter1(pfoToLArTPCMap.find(pPfoToShift1));
            PfoToLArTPCMap::const_iterator iter2(pfoToLArTPCMap.find(pPfoToShift2));

            if ((iter1 != pfoToLArTPCMap.end()) && (iter2 != pfoToLArTPCMap.end()))
            {
                const LArTPC *const pLArTPC1(iter1->second);
                const LArTPC *const pLArTPC2(iter2->second);

                if (((twoDHitsPfo1 + twoDHitsPfo2) > totalHitsStitched) && LArStitchingHelper::CanTPCsBeStitched(*pLArTPC1, *pLArTPC2))
                {
                    totalHitsStitched = twoDHitsPfo1 + twoDHitsPfo2;
                    reducedPfoVector.clear();
                    reducedPfoVector.push_back(pPfoToShift1);
                    reducedPfoVector.push_back(pPfoToShift2);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *StitchingCosmicRayMergingTool::GetClosestPfo(const ParticleFlowObject *const pPfoToEnlarge, const PfoVector &pfoVector) const
{
    const ParticleFlowObject *pClosestPfo(*pfoVector.begin());

    float minDistance(std::numeric_limits<float>::max());
    for (const ParticleFlowObject *const pPfoToShift : pfoVector)
    {
        if (pPfoToShift == pPfoToEnlarge)
            return pPfoToEnlarge;

        const float thisDistance(LArPfoHelper::GetTwoDSeparation(pPfoToShift, pPfoToEnlarge));
        if (thisDistance < minDistance)
        {
            minDistance = thisDistance;
            pClosestPfo = pPfoToShift;
        }
    }

    return pClosestPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::FindStitchedLArTPCs(const PfoVector &pfoVector, const PfoToLArTPCMap &pfoToLArTPCMap,
    LArTPCPair &stitchedLArTPCs) const
{
    for (const ParticleFlowObject *const pPfoToShift : pfoVector)
    {
        PfoToLArTPCMap::const_iterator iter(pfoToLArTPCMap.find(pPfoToShift));

        if (iter != pfoToLArTPCMap.end())
        {
            const LArTPC *const pLArTPC(iter->second);

            if (pLArTPC && (stitchedLArTPCs.first == pLArTPC || stitchedLArTPCs.second == pLArTPC))
                continue;

            if (!stitchedLArTPCs.first)
            {
                stitchedLArTPCs.first = pLArTPC;
            }
            else if (!stitchedLArTPCs.second)
            {
                stitchedLArTPCs.second = pLArTPC;
            }
            else
            {
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }
        }
        else
        {
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CalculateX0(const PfoToLArTPCMap &pfoToLArTPCMap, const ThreeDPointingClusterMap &pointingClusterMap,
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

        for (PfoVector::const_iterator iter2 = iter1; iter2 != iterEnd; ++iter2)
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
            unsigned int count(0);
            try
            {
                ++count;
                
                LArStitchingHelper::GetClosestVertices(*pLArTPC1, *pLArTPC2, pointingCluster1, pointingCluster2,
                    pointingVertex1, pointingVertex2);
                
                PfoToPointingVertexMatrix::iterator pfoToPointingVertexMatrixIter1(pfoToPointingVertexMatrix.find(pPfo1));
                if (pfoToPointingVertexMatrixIter1 == pfoToPointingVertexMatrix.end())
                {
                    PfoToPointingVertexMap pfoToPointingVertexMap({{pPfo2, pointingVertex1}});
                    (void) pfoToPointingVertexMatrix.insert(PfoToPointingVertexMatrix::value_type(pPfo1, pfoToPointingVertexMap));
                }
                else
                {
                    PfoToPointingVertexMap pfoToPointingVertexMap(pfoToPointingVertexMatrixIter1->second);
                    PfoToPointingVertexMap::iterator pfoToPointingVertexMapIter(pfoToPointingVertexMap.find(pPfo2));
                    if(pfoToPointingVertexMapIter == pfoToPointingVertexMap.end())
                    {
                        (void) pfoToPointingVertexMap.insert(PfoToPointingVertexMap::value_type(pPfo2, pointingVertex1));
                    }
                    else
                    {
                        if ((pfoToPointingVertexMapIter->second.GetPosition() - pointingVertex1.GetPosition()).GetMagnitude() > std::numeric_limits<float>::epsilon())
                            throw StatusCodeException(STATUS_CODE_FAILURE);;
                    }
                }

                PfoToPointingVertexMatrix::iterator pfoToPointingVertexMatrixIter2(pfoToPointingVertexMatrix.find(pPfo2));
                if (pfoToPointingVertexMatrixIter2 == pfoToPointingVertexMatrix.end())
                {
                    PfoToPointingVertexMap pfoToPointingVertexMap({{pPfo1, pointingVertex2}});
                    (void) pfoToPointingVertexMatrix.insert(PfoToPointingVertexMatrix::value_type(pPfo2, pfoToPointingVertexMap));
                }
                else
                {
                    PfoToPointingVertexMap pfoToPointingVertexMap(pfoToPointingVertexMatrixIter2->second);
                    PfoToPointingVertexMap::iterator pfoToPointingVertexMapIter(pfoToPointingVertexMap.find(pPfo1));
                    if(pfoToPointingVertexMapIter == pfoToPointingVertexMap.end())
                    {
                        (void) pfoToPointingVertexMap.insert(PfoToPointingVertexMap::value_type(pPfo1, pointingVertex2));
                    }
                    else
                    {
                        if ((pfoToPointingVertexMapIter->second.GetPosition() - pointingVertex2.GetPosition()).GetMagnitude() > std::numeric_limits<float>::epsilon())
                            throw StatusCodeException(STATUS_CODE_FAILURE);;
                    }
                }

                float thisX0(LArStitchingHelper::CalculateX0(*pLArTPC1, *pLArTPC2, pointingVertex1, pointingVertex2));
                thisX0 *= (count % 2 == 0) ? 1.f : -1.f; 
                sumX += thisX0; sumN += 1.f;
            }
            catch (const pandora::StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    std::cout << "StitchingCosmicRayMergingTool: Attempting to stitch a pfo using multiple cluster pointing vertices" << std::endl;
            }
        }
    }

    if ((sumN < std::numeric_limits<float>::epsilon()) || (std::fabs(sumX) < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    x0 = (sumX / sumN);
}

//------------------------------------------------------------------------------------------------------------------------------------------
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThreeDStitchingMode", m_useXcoordinate));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AlwaysApplyT0Calculation", m_alwaysApplyT0Calculation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HalfWindowLayers", m_halfWindowLayers));

    float minLength(std::sqrt(m_minLengthSquared));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPfoLength", minLength));
    m_minLengthSquared = minLength * minLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacementX", m_maxLongitudinalDisplacementX));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LooserMinCosRelativeAngle", m_relaxCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LooserMaxTransverseDisplacement", m_relaxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNCaloHits3D", m_minNCaloHits3D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_fileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
