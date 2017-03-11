/**
 *  @file   LArContent/src/LArStitching/StitchingCosmicRayMergingTool.cc
 *
 *  @brief  Implementation of the stitching pfo merging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"
#include "larpandoracontent/LArStitching/LArStitchingHelper.h"
#include "larpandoracontent/LArStitching/StitchingCosmicRayMergingTool.h"

using namespace pandora;

namespace lar_content
{

StitchingCosmicRayMergingTool::StitchingCosmicRayMergingTool() :
    m_useXcoordinate(false),
    m_halfWindowLayers(30),
    m_minLengthSquared(50.f),
    m_minCosRelativeAngle(0.966),
    m_maxLongitudinalDisplacementX(15.f),
    m_maxTransverseDisplacement(5.f),
    m_relaxCosRelativeAngle(0.906),
    m_relaxTransverseDisplacement(2.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Ensure that at least one Pfo knows its drift volume
    const PfoToVolumeIdMap &pfoToVolumeIdMap(stitchingInfo.m_pfoToVolumeIdMap);
    if (pfoToVolumeIdMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Get the complete list of Pfos (from all drift volumes)
    const PfoList *pMultiPfoList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMultiPfoList));

    // Build map of primary Pfos (these must be in the top level of the Pfo hierarchy and flagged as track-like)
    PfoList primaryPfos;
    this->SelectPrimaryPfos(pMultiPfoList, primaryPfos);

    // Build 3D pointing clusters from primary Pfos
    ThreeDPointingClusterMap pointingClusterMap;
    this->BuildPointingClusterMaps(primaryPfos, pointingClusterMap);

    // Separate Pfos by their volume ID
    VolumeIdToPfoMap volumeIdToPfoMap;
    this->BuildVolumeMaps(primaryPfos, pfoToVolumeIdMap, volumeIdToPfoMap);

    // Use pointing clusters to form associations between Pfos
    PfoAssociationMatrix pfoAssociationMatrix;
    this->CreatePfoMatches(volumeIdToPfoMap, pointingClusterMap, pfoAssociationMatrix);

    // Select the best associations between Pfos
    PfoMergeMap pfoSelectedMatches;
    this->SelectPfoMatches(pfoAssociationMatrix, pfoSelectedMatches);

    // Create an initial map of Pfo merges to be made
    PfoMergeMap pfoSelectedMerges;
    this->SelectPfoMerges(pfoSelectedMatches, pfoSelectedMerges);

    // Identify the vertex Pfo and re-order the map of Pfo merges accordingly
    PfoMergeMap pfoOrderedMerges;
    this->OrderPfoMerges(pfoToVolumeIdMap, pointingClusterMap, pfoSelectedMerges, pfoOrderedMerges);

    // Apply X0 corrections, and then stitch together Pfos
    this->StitchPfos(pAlgorithm, pointingClusterMap, pfoOrderedMerges, stitchingInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SelectPrimaryPfos(const PfoList *pInputPfoList, PfoList &outputPfoList) const
{
    for (const ParticleFlowObject *const pPfo : *pInputPfoList)
    {
        if (!LArPfoHelper::IsFinalState(pPfo) || !LArPfoHelper::IsTrack(pPfo))
            continue;

        outputPfoList.push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::BuildVolumeMaps(const PfoList &inputPfoList, const PfoToVolumeIdMap &pfoToVolumeIdMap,
    VolumeIdToPfoMap &volumeIdToPfoMap) const
{
    for (const ParticleFlowObject *const pPfo : inputPfoList)
    {
        PfoToVolumeIdMap::const_iterator iter = pfoToVolumeIdMap.find(pPfo);

        if (pfoToVolumeIdMap.end() == iter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const int volumeID(iter->second);

        try
        {
            const VolumeInfo &checkVolume(MultiPandoraApi::GetVolumeInfo(volumeID));

            if (volumeID != checkVolume.GetIdNumber())
                throw StatusCodeException(STATUS_CODE_FAILURE);

            (void) volumeIdToPfoMap[volumeID].push_back(pPfo);
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }

/*
// Display volume maps here
PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f));
for (VolumeIdToPfoMap::const_iterator iter = volumeIdToPfoMap.begin(), iterEnd = volumeIdToPfoMap.end(); iter != iterEnd; ++iter)
{
const PfoList &pfoList(iter->second);
PfoList myPfoList(pfoList.begin(), pfoList.end());
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &myPfoList, "PrimaryPfos", BLUE));
}
PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::BuildPointingClusterMaps(const PfoList &inputPfoList, ThreeDPointingClusterMap &pointingClusterMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (const ParticleFlowObject *const pPfo : inputPfoList)
    {
        try
        {
            ClusterList clusterList;
            LArPfoHelper::GetThreeDClusterList(pPfo, clusterList);

            if (1 != clusterList.size())
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            const Cluster *pCluster = *(clusterList.begin());
            const ThreeDSlidingFitResult slidingFitResult(pCluster, m_halfWindowLayers, slidingFitPitch);
            const LArPointingCluster pointingCluster(slidingFitResult);

            (void) pointingClusterMap.insert(ThreeDPointingClusterMap::value_type(pPfo, pointingCluster));
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CreatePfoMatches(const VolumeIdToPfoMap &volumeIdToPfoMap, const ThreeDPointingClusterMap &pointingClusterMap,
    PfoAssociationMatrix &pfoAssociationMatrix) const
{
    for (VolumeIdToPfoMap::const_iterator iter1A = volumeIdToPfoMap.begin(), iterEnd1A = volumeIdToPfoMap.end(); iter1A != iterEnd1A; ++iter1A)
    {
        const VolumeInfo &pfoVolume1(MultiPandoraApi::GetVolumeInfo(iter1A->first));
        const PfoList &pfoList1(iter1A->second);

        for (VolumeIdToPfoMap::const_iterator iter2A = iter1A, iterEnd2A = volumeIdToPfoMap.end(); iter2A != iterEnd2A; ++iter2A)
        {
            const VolumeInfo &pfoVolume2(MultiPandoraApi::GetVolumeInfo(iter2A->first));
            const PfoList &pfoList2(iter2A->second);

            if (!LArStitchingHelper::CanVolumesBeStitched(pfoVolume1, pfoVolume2))
                continue;

            for (PfoList::const_iterator iter1 = pfoList1.begin(), iterEnd1 = pfoList1.end(); iter1 != iterEnd1; ++iter1)
            {
                const ParticleFlowObject *const pPfo1 = *iter1;

                for (PfoList::const_iterator iter2 = pfoList2.begin(), iterEnd2 = pfoList2.end(); iter2 != iterEnd2; ++iter2)
                {
                    const ParticleFlowObject *const pPfo2 = *iter2;

                    this->CreatePfoMatches(pfoVolume1, pfoVolume2, pPfo1, pPfo2, pointingClusterMap, pfoAssociationMatrix);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CreatePfoMatches(const VolumeInfo &pfoVolume1, const VolumeInfo &pfoVolume2,
    const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2,
    const ThreeDPointingClusterMap &pointingClusterMap, PfoAssociationMatrix &pfoAssociationMatrix) const
{
    // Get the pointing cluster corresponding to each of these Pfos
    ThreeDPointingClusterMap::const_iterator iter1 = pointingClusterMap.find(pPfo1);
    ThreeDPointingClusterMap::const_iterator iter2 = pointingClusterMap.find(pPfo2);

    if (pointingClusterMap.end() == iter1 || pointingClusterMap.end() == iter2)
        return;

    const LArPointingCluster pointingCluster1(iter1->second);
    const LArPointingCluster pointingCluster2(iter2->second);

    // Check length of pointing clusters
    if (pointingCluster1.GetLengthSquared() < m_minLengthSquared || pointingCluster2.GetLengthSquared() < m_minLengthSquared)
        return;

    // Get closest pair of vertices
    LArPointingCluster::Vertex pointingVertex1, pointingVertex2;

    try
    {
        LArStitchingHelper::GetClosestVertices(pfoVolume1, pfoVolume2, pointingCluster1, pointingCluster2, pointingVertex1, pointingVertex2);
    }
    catch (pandora::StatusCodeException& )
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
    catch (pandora::StatusCodeException& )
    {
        return;
    }

    // Selection cuts on longitudinal impact parameters
    const float minL(-1.f);
    const float dXdL1(m_useXcoordinate ? pX1 :
        (1.f - pX1 * pX1 > std::numeric_limits<float>::epsilon()) ? pX1 / std::sqrt(1.f - pX1 * pX1) : minL);
    const float dXdL2(m_useXcoordinate ? pX2 :
        (1.f - pX2 * pX2 > std::numeric_limits<float>::epsilon()) ? pX2 / std::sqrt(1.f - pX2 * pX2) : minL);
    const float maxL1(m_maxLongitudinalDisplacementX / dXdL1);
    const float maxL2(m_maxLongitudinalDisplacementX / dXdL2);

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


/*
//Display matches here
std::cout << " *** StitchingCosmicRayMergingTool::MatchParticles(...) *** " << std::endl
          << "     pPfo1=" << pPfo1 << " pPfo2=" << pPfo2 << std::endl
          << "     cosTheta=" << cosRelativeAngle << std::endl
          << "     rL1=" << rL1 << ", rT1=" << rT1 << ", rL2=" << rL2 << ", rT2=" << rT2 << std::endl
          << "     Length1=" << std::sqrt(particleLength1) << " Length2=" << std::sqrt(particleLength2) << std::endl;
PfoList myPfoList1, myPfoList2;
myPfoList1.push_back(pPfo1);
myPfoList2.push_back(pPfo2);
PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f));
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &myPfoList1, "Pfo1", BLUE));
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &myPfoList2, "Pfo2", GREEN));
PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pointingVertex1.GetPosition(), "Vertex1", RED, 3.5f));
PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pointingVertex2.GetPosition(), "Vertex2", RED, 3.5f));
PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::SelectPfoMatches(const PfoAssociationMatrix &pfoAssociationMatrix, PfoMergeMap &pfoMatches) const
{
    // First step: loop over association matrix and find best associations A -> X and B -> Y
    // =====================================================================================
    PfoAssociationMatrix bestAssociationMatrix;

    for (PfoAssociationMatrix::const_iterator iter1 = pfoAssociationMatrix.begin(), iterEnd1 = pfoAssociationMatrix.end();
        iter1 != iterEnd1; ++iter1)
    {
        const ParticleFlowObject *const pPfo1 = iter1->first;
        const PfoAssociationMap &pfoAssociationMap(iter1->second);

        const ParticleFlowObject *pBestPfoInner = NULL;
        PfoAssociation bestAssociationInner(PfoAssociation::UNDEFINED, PfoAssociation::UNDEFINED, 0.f);

        const ParticleFlowObject *pBestPfoOuter = NULL;
        PfoAssociation bestAssociationOuter(PfoAssociation::UNDEFINED, PfoAssociation::UNDEFINED, 0.f);

        for (PfoAssociationMap::const_iterator iter2 = pfoAssociationMap.begin(), iterEnd2 = pfoAssociationMap.end();
            iter2 != iterEnd2; ++iter2)
       {
            const ParticleFlowObject *const pPfo2 = iter2->first;
            const PfoAssociation &pfoAssociation(iter2->second);

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
    for (PfoAssociationMatrix::const_iterator iter3 = bestAssociationMatrix.begin(), iterEnd3 = bestAssociationMatrix.end();
        iter3 != iterEnd3; ++iter3)
    {
        const ParticleFlowObject *const pParentPfo = iter3->first;
        const PfoAssociationMap &parentAssociationMap(iter3->second);

        for (PfoAssociationMap::const_iterator iter4 = parentAssociationMap.begin(), iterEnd4 = parentAssociationMap.end();
            iter4 != iterEnd4; ++iter4)
        {
            const ParticleFlowObject *const pDaughterPfo = iter4->first;
            const PfoAssociation &parentToDaughterAssociation(iter4->second);

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
    PfoList vetoList;

    for (PfoMergeMap::const_iterator iter = pfoMatches.begin(), iterEnd = pfoMatches.end(); iter != iterEnd; ++iter)
    {
        const PfoList &pfoList(iter->second);

        for (const ParticleFlowObject *const pSeedPfo : pfoList)
        {
            if (vetoList.end() != std::find(vetoList.begin(), vetoList.end(), pSeedPfo))
                continue;

            PfoList mergeList;
            this->CollectAssociatedPfos(pSeedPfo, pSeedPfo, pfoMatches, vetoList, mergeList);

            vetoList.push_back(pSeedPfo);
            pfoMerges[pSeedPfo].push_back(pSeedPfo);

            for (const ParticleFlowObject *const pAssociatedPfo : mergeList)
            {
                // Throw an exception if this particle has already been counted
                if (vetoList.end() != std::find(vetoList.begin(), vetoList.end(), pAssociatedPfo) ||
                    pfoMerges[pSeedPfo].end() != std::find(pfoMerges[pSeedPfo].begin(), pfoMerges[pSeedPfo].end(), pAssociatedPfo))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                vetoList.push_back(pAssociatedPfo);
                pfoMerges[pSeedPfo].push_back(pAssociatedPfo);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CollectAssociatedPfos(const ParticleFlowObject *const pSeedPfo, const ParticleFlowObject *const pCurrentPfo,
    const PfoMergeMap &pfoMergeMap, const PfoList &vetoList, PfoList &associatedList) const
{
    if (vetoList.end() != std::find(vetoList.begin(), vetoList.end(), pCurrentPfo))
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

        this->CollectAssociatedPfos(pSeedPfo, pAssociatedPfo, pfoMergeMap, vetoList, associatedList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::OrderPfoMerges(const PfoToVolumeIdMap &pfoToVolumeIdMap, const ThreeDPointingClusterMap &pointingClusterMap,
    const PfoMergeMap &inputPfoMerges, PfoMergeMap &outputPfoMerges) const
{
    for (PfoMergeMap::const_iterator iter = inputPfoMerges.begin(), iterEnd = inputPfoMerges.end(); iter != iterEnd; ++iter)
    {
        float bestLength(0.f);
        const ParticleFlowObject *pVertexPfo = NULL;
        const PfoList &pfoList = iter->second;

        for (PfoList::const_iterator iter1 = pfoList.begin(), iterEnd1 = pfoList.end(); iter1 != iterEnd1; ++iter1)
        {
            const ParticleFlowObject *const pPfo1 = *iter1;

            PfoToVolumeIdMap::const_iterator iter1A = pfoToVolumeIdMap.find(pPfo1);
            ThreeDPointingClusterMap::const_iterator iter1B = pointingClusterMap.find(pPfo1);

            if (pfoToVolumeIdMap.end() == iter1A || pointingClusterMap.end() == iter1B)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const VolumeInfo &pfoVolume1(MultiPandoraApi::GetVolumeInfo(iter1A->second));
            const LArPointingCluster &pointingCluster1(iter1B->second);

            for (PfoList::const_iterator iter2 = iter1, iterEnd2 = pfoList.end(); iter2 != iterEnd2; ++iter2)
            {
                const ParticleFlowObject *const pPfo2 = *iter2;

                PfoToVolumeIdMap::const_iterator iter2A = pfoToVolumeIdMap.find(pPfo2);
                ThreeDPointingClusterMap::const_iterator iter2B = pointingClusterMap.find(pPfo2);

                if (pfoToVolumeIdMap.end() == iter2A || pointingClusterMap.end() == iter2B)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const VolumeInfo &pfoVolume2(MultiPandoraApi::GetVolumeInfo(iter2A->second));
                const LArPointingCluster &pointingCluster2(iter2B->second);

                if (pfoVolume1.GetIdNumber() == pfoVolume2.GetIdNumber())
                    continue;

                const float thisLength(LArStitchingHelper::GetVolumeDisplacement(pfoVolume1, pfoVolume2));

                if (thisLength < bestLength)
                    continue;

                bestLength = thisLength;

                try
                {
                    pVertexPfo = NULL;

                    LArPointingCluster::Vertex nearVertex1, nearVertex2;
                    LArStitchingHelper::GetClosestVertices(pfoVolume1, pfoVolume2, pointingCluster1, pointingCluster2,
                        nearVertex1, nearVertex2);

                    const LArPointingCluster::Vertex &farVertex1(nearVertex1.IsInnerVertex() ? pointingCluster1.GetOuterVertex() : pointingCluster1.GetInnerVertex());
                    const LArPointingCluster::Vertex &farVertex2(nearVertex2.IsInnerVertex() ? pointingCluster2.GetOuterVertex() : pointingCluster2.GetInnerVertex());
                    const float deltaY(farVertex1.GetPosition().GetY() - farVertex2.GetPosition().GetY());

                    if (std::fabs(deltaY) < std::numeric_limits<float>::epsilon())
                         throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                    pVertexPfo = ((deltaY > 0.f) ? pPfo1 : pPfo2);
                }
                catch (pandora::StatusCodeException& )
                {
                }
            }
        }

        if (pVertexPfo)
            outputPfoMerges[pVertexPfo].insert(outputPfoMerges[pVertexPfo].begin(), pfoList.begin(), pfoList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::StitchPfos(const StitchingAlgorithm *const pAlgorithm, const ThreeDPointingClusterMap &pointingClusterMap,
    const PfoMergeMap &pfoMerges, StitchingAlgorithm::StitchingInfo &stitchingInfo) const
{
    const StitchingAlgorithm::PfoToVolumeIdMap &pfoToVolumeIdMap(stitchingInfo.m_pfoToVolumeIdMap);

    for (PfoMergeMap::const_iterator iter1 = pfoMerges.begin(), iterEnd1 = pfoMerges.end(); iter1 != iterEnd1; ++iter1)
    {
        const ParticleFlowObject *const pPfoToEnlarge = iter1->first;
        const PfoList &pfoList = iter1->second;

        if (!m_useXcoordinate)
        {
            float x0(0.f);

            try
            {
                this->CalculateX0(pfoToVolumeIdMap, pointingClusterMap, pfoList, x0);
            }
            catch (pandora::StatusCodeException& )
            {
                continue;
            }
/*
// Display X0 information
PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f));
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "PfoList", BLUE));
for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
{
const ParticleFlowObject *const pPfo = *iter;
PfoToVolumeIdMap::const_iterator iterA = pfoToVolumeIdMap.find(pPfo);
ThreeDPointingClusterMap::const_iterator iterB = pointingClusterMap.find(pPfo);
const VolumeInfo &pfoVolume(MultiPandoraApi::GetVolumeInfo(iterA->second));
const LArPointingCluster &pointingCluster(iterB->second);
const CartesianVector xOffset(pfoVolume.IsDriftInPositiveX() ? -x0 : x0, 0.f, 0.f);
const CartesianVector innerVertex((pointingCluster.GetInnerVertex().GetPosition() + xOffset));
const CartesianVector outerVertex((pointingCluster.GetOuterVertex().GetPosition() + xOffset));
PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &innerVertex, "Vertex1", RED, 3.5f));
PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &outerVertex, "Vertex2", BLUE, 3.5f));
}
PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
*/
            for (PfoList::const_iterator iter2 = pfoList.begin(), iterEnd2 = pfoList.end(); iter2 != iterEnd2; ++iter2)
            {
                const ParticleFlowObject *const pPfoToShift = *iter2;

                pAlgorithm->ShiftPfoHierarchy(pPfoToShift, stitchingInfo, x0);
            }
        }

        for (PfoList::const_iterator iter2 = pfoList.begin(), iterEnd2 = pfoList.end(); iter2 != iterEnd2; ++iter2)
        {
            const ParticleFlowObject *const pPfoToDelete = *iter2;

            if (pPfoToEnlarge == pPfoToDelete)
                continue;
/*
//Display Pfos Before Merge
PfoList enlargePfoList, deletePfoList;
enlargePfoList.push_back(pPfoToEnlarge);
deletePfoList.push_back(pPfoToDelete);
PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f));
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &enlargePfoList, "enlargePfoList", BLUE));
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &deletePfoList, "deletePfoList", GREEN));
PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
*/
            pAlgorithm->StitchPfos(pPfoToEnlarge, pPfoToDelete, stitchingInfo);
/*
//Display Pfos After Merge
PfoList mergedPfoList;
mergedPfoList.push_back(pPfoToEnlarge);
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &mergedPfoList, "mergedPfoList", RED));
PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
*/
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StitchingCosmicRayMergingTool::CalculateX0(const PfoToVolumeIdMap &pfoToVolumeIdMap, const ThreeDPointingClusterMap &pointingClusterMap,
    const PfoList &pfoList, float &x0) const
{
    float SumX(0.f), SumN(0.f);

    for (PfoList::const_iterator iter1 = pfoList.begin(), iterEnd1 = pfoList.end(); iter1 != iterEnd1; ++iter1)
    {
        const ParticleFlowObject *const pPfo1 = *iter1;

        PfoToVolumeIdMap::const_iterator iter1A = pfoToVolumeIdMap.find(pPfo1);
        ThreeDPointingClusterMap::const_iterator iter1B = pointingClusterMap.find(pPfo1);

        if (pfoToVolumeIdMap.end() == iter1A || pointingClusterMap.end() == iter1B)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        for (PfoList::const_iterator iter2 = iter1, iterEnd2 = pfoList.end(); iter2 != iterEnd2; ++iter2)
        {
            const ParticleFlowObject *const pPfo2 = *iter2;

            PfoToVolumeIdMap::const_iterator iter2A = pfoToVolumeIdMap.find(pPfo2);
            ThreeDPointingClusterMap::const_iterator iter2B = pointingClusterMap.find(pPfo2);

            if (pfoToVolumeIdMap.end() == iter2A || pointingClusterMap.end() == iter2B)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            // Check that this pair of Pfos can be stitched
            const VolumeInfo &pfoVolume1(MultiPandoraApi::GetVolumeInfo(iter1A->second));
            const VolumeInfo &pfoVolume2(MultiPandoraApi::GetVolumeInfo(iter2A->second));

            if (!LArStitchingHelper::CanVolumesBeStitched(pfoVolume1, pfoVolume2))
                continue;

            // Calculate X0 for the closest pair of vertices
            const LArPointingCluster &pointingCluster1(iter1B->second);
            const LArPointingCluster &pointingCluster2(iter2B->second);

            LArPointingCluster::Vertex pointingVertex1, pointingVertex2;

            try
            {
                LArStitchingHelper::GetClosestVertices(pfoVolume1, pfoVolume2, pointingCluster1, pointingCluster2,
                    pointingVertex1, pointingVertex2);

                const float thisX0(LArStitchingHelper::CalculateX0(pfoVolume1, pfoVolume2, pointingVertex1, pointingVertex2));

                SumX += thisX0; SumN += 1.f;
            }
            catch (pandora::StatusCodeException& )
            {
            }
        }
    }

    if (SumN < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    x0 = (SumX / SumN);
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

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
