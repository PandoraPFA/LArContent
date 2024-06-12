/**
 *  @file   larpandoracontent/LArControlFlow/TestBeamCosmicRayTaggingTool.cc
 *
 *  @brief  Implementation of the cosmic-ray tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/TestBeamCosmicRayTaggingTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

TestBeamCosmicRayTaggingTool::TestBeamCosmicRayTaggingTool() :
    m_angularUncertainty(5.f),
    m_positionalUncertainty(3.f),
    m_maxAssociationDist(3.f * 18.f),
    m_minimumHits(15),
    m_inTimeMargin(5.f),
    m_inTimeMaxX0(1.f),
    m_marginY(20.f),
    m_tagTopEntering(true),
    m_tagTopToBottom(true),
    m_tagOutOfTime(true),
    m_tagInVetoedTPCs(false),
    m_face_Xa(0.f),
    m_face_Xc(0.f),
    m_face_Yb(0.f),
    m_face_Yt(0.f),
    m_face_Zu(0.f),
    m_face_Zd(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamCosmicRayTaggingTool::Initialize()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos, const MasterAlgorithm *const /*pAlgorithm*/)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // TODO First time only, TODO Refactor with master algorithm
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);

    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    m_face_Xa = parentMinX;
    m_face_Xc = parentMaxX;
    m_face_Yb = parentMinY;
    m_face_Yt = parentMaxY;
    m_face_Zu = parentMinZ;
    m_face_Zd = parentMaxZ;

    PfoToPfoListMap pfoAssociationMap;
    this->GetPfoAssociations(parentCosmicRayPfos, pfoAssociationMap);

    PfoToSliceIdMap pfoToSliceIdMap;
    this->SliceEvent(parentCosmicRayPfos, pfoAssociationMap, pfoToSliceIdMap);

    CRCandidateList candidates;
    this->GetCRCandidates(parentCosmicRayPfos, pfoToSliceIdMap, candidates);

    // Create a map that will be used to find all cosmic ray pfos
    PfoToBoolMap pfoToIsCosmicRayMap;
    for (const CRCandidate &candidate : candidates)
    {
        pfoToIsCosmicRayMap[candidate.m_pPfo] = false;
    }

    // Tag all particles appearing to enter from the top
    if (m_tagTopEntering)
        this->CheckIfTopEntering(candidates, pfoToIsCosmicRayMap);

    // Tag particles going from top to bottom
    if (m_tagTopToBottom)
        this->CheckIfTopToBottom(candidates, pfoToIsCosmicRayMap);

    // Tag particles with t0 != 0 or appears outside the detector with t0 = 0
    if (m_tagOutOfTime)
        this->CheckIfOutOfTime(candidates, pfoToIsCosmicRayMap);

    // Tag particles with their highest point in a given list of TPCs
    if (m_tagInVetoedTPCs)
        this->CheckIfInVetoedTPC(candidates, pfoToIsCosmicRayMap);

    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        if (!pfoToIsCosmicRayMap.at(pPfo))
            ambiguousPfos.push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TestBeamCosmicRayTaggingTool::GetValid3DCluster(const ParticleFlowObject *const pPfo, const Cluster *&pCluster3D) const
{
    pCluster3D = nullptr;
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

    // ATTN Only uses first cluster with hits of type TPC_3D
    if (clusters3D.empty() || (clusters3D.front()->GetNCaloHits() < m_minimumHits))
        return false;

    pCluster3D = clusters3D.front();
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::GetPfoAssociations(const PfoList &parentCosmicRayPfos, PfoToPfoListMap &pfoAssociationMap) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pFirstLArTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float layerPitch(pFirstLArTPC->GetWirePitchW());

    PfoToSlidingFitsMap pfoToSlidingFitsMap;

    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        const pandora::Cluster *pCluster(nullptr);
        if (!this->GetValid3DCluster(pPfo, pCluster) || !pCluster)
            continue;

        (void)pfoToSlidingFitsMap.insert(PfoToSlidingFitsMap::value_type(pPfo,
            std::make_pair(ThreeDSlidingFitResult(pCluster, 5, layerPitch), ThreeDSlidingFitResult(pCluster, 100, layerPitch)))); // TODO Configurable
    }

    for (const ParticleFlowObject *const pPfo1 : parentCosmicRayPfos)
    {
        PfoToSlidingFitsMap::const_iterator iter1(pfoToSlidingFitsMap.find(pPfo1));
        if (pfoToSlidingFitsMap.end() == iter1)
            continue;

        const ThreeDSlidingFitResult &fitPos1(iter1->second.first), &fitDir1(iter1->second.second);

        for (const ParticleFlowObject *const pPfo2 : parentCosmicRayPfos)
        {
            if (pPfo1 == pPfo2)
                continue;

            PfoToSlidingFitsMap::const_iterator iter2(pfoToSlidingFitsMap.find(pPfo2));
            if (pfoToSlidingFitsMap.end() == iter2)
                continue;

            const ThreeDSlidingFitResult &fitPos2(iter2->second.first), &fitDir2(iter2->second.second);

            // TODO Use existing LArPointingClusters and IsEmission/IsNode logic, for consistency
            if (!(this->CheckAssociation(fitPos1.GetGlobalMinLayerPosition(), fitDir1.GetGlobalMinLayerDirection() * -1.f,
                      fitPos2.GetGlobalMinLayerPosition(), fitDir2.GetGlobalMinLayerDirection() * -1.f) ||
                    this->CheckAssociation(fitPos1.GetGlobalMinLayerPosition(), fitDir1.GetGlobalMinLayerDirection() * -1.f,
                        fitPos2.GetGlobalMaxLayerPosition(), fitDir2.GetGlobalMaxLayerDirection()) ||
                    this->CheckAssociation(fitPos1.GetGlobalMaxLayerPosition(), fitDir1.GetGlobalMaxLayerDirection(),
                        fitPos2.GetGlobalMinLayerPosition(), fitDir2.GetGlobalMinLayerDirection() * -1.f) ||
                    this->CheckAssociation(fitPos1.GetGlobalMaxLayerPosition(), fitDir1.GetGlobalMaxLayerDirection(),
                        fitPos2.GetGlobalMaxLayerPosition(), fitDir2.GetGlobalMaxLayerDirection())))
            {
                continue;
            }

            PfoList &pfoList1(pfoAssociationMap[pPfo1]), &pfoList2(pfoAssociationMap[pPfo2]);

            if (pfoList1.end() == std::find(pfoList1.begin(), pfoList1.end(), pPfo2))
                pfoList1.push_back(pPfo2);

            if (pfoList2.end() == std::find(pfoList2.begin(), pfoList2.end(), pPfo1))
                pfoList2.push_back(pPfo1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TestBeamCosmicRayTaggingTool::CheckAssociation(
    const CartesianVector &endPoint1, const CartesianVector &endDir1, const CartesianVector &endPoint2, const CartesianVector &endDir2) const
{
    // TODO This function needs significant tidying and checks (variable names must be self describing, check deltaTheta, etc.)
    const CartesianVector n(endDir1.GetUnitVector());
    const CartesianVector m(endDir2.GetUnitVector());
    const CartesianVector a(endPoint2 - endPoint1);
    const float b(n.GetDotProduct(m));

    // Parallel lines never meet
    if (std::fabs(b - 1.f) < std::numeric_limits<float>::epsilon())
        return false;

    // Distance from endPoint1 along endDir1 to the point of closest approach
    const float lambda((n - m * b).GetDotProduct(a) / (1.f - b * b));

    // Distance from endPoint2 along endDir2 to the point of closest approach
    const float mu((n * b - m).GetDotProduct(a) / (1.f - b * b));

    // Calculate the maximum "vertex uncertainty"
    const float deltaTheta(m_angularUncertainty * M_PI / 180.f);
    const float maxVertexUncertainty(m_maxAssociationDist * std::sin(deltaTheta) + m_positionalUncertainty);

    // Ensure that the distances to the point of closest approch are within the limits
    if ((lambda < -maxVertexUncertainty) || (mu < -maxVertexUncertainty) || (lambda > m_maxAssociationDist + maxVertexUncertainty) ||
        (mu > m_maxAssociationDist + maxVertexUncertainty))
    {
        return false;
    }

    // Ensure point of closest approach is within detector
    const CartesianVector impactPosition((endPoint1 + n * lambda + endPoint2 + m * mu) * 0.5f);

    if ((impactPosition.GetX() < m_face_Xa - maxVertexUncertainty) || (impactPosition.GetX() > m_face_Xc + maxVertexUncertainty) ||
        (impactPosition.GetY() < m_face_Yb - maxVertexUncertainty) || (impactPosition.GetY() > m_face_Yt + maxVertexUncertainty) ||
        (impactPosition.GetZ() < m_face_Zu - maxVertexUncertainty) || (impactPosition.GetZ() > m_face_Zd + maxVertexUncertainty))
    {
        return false;
    }

    // Compare distance of closest approach and max uncertainty in impact parameter
    const float maxImpactDist(std::sin(deltaTheta) * (std::fabs(mu) + std::fabs(lambda)) + m_positionalUncertainty);
    const CartesianVector d(a - n * lambda + m * mu);

    return (d.GetMagnitude() < maxImpactDist);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::SliceEvent(const PfoList &parentCosmicRayPfos, const PfoToPfoListMap &pfoAssociationMap, PfoToSliceIdMap &pfoToSliceIdMap) const
{
    SliceList sliceList;

    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        bool isAlreadyInSlice(false);

        for (const PfoList &slice : sliceList)
        {
            if (std::find(slice.begin(), slice.end(), pPfo) != slice.end())
            {
                isAlreadyInSlice = true;
                break;
            }
        }

        if (!isAlreadyInSlice)
        {
            sliceList.push_back(PfoList());
            this->FillSlice(pPfo, pfoAssociationMap, sliceList.back());
        }
    }

    unsigned int sliceId(0);
    for (const PfoList &slice : sliceList)
    {
        for (const ParticleFlowObject *const pPfo : slice)
        {
            if (!pfoToSliceIdMap.insert(PfoToSliceIdMap::value_type(pPfo, sliceId)).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }

        ++sliceId;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::FillSlice(const ParticleFlowObject *const pPfo, const PfoToPfoListMap &pfoAssociationMap, PfoList &slice) const
{
    if (std::find(slice.begin(), slice.end(), pPfo) != slice.end())
        return;

    slice.push_back(pPfo);

    PfoToPfoListMap::const_iterator iter(pfoAssociationMap.find(pPfo));

    if (pfoAssociationMap.end() != iter)
    {
        for (const ParticleFlowObject *const pAssociatedPfo : iter->second)
            this->FillSlice(pAssociatedPfo, pfoAssociationMap, slice);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::GetCRCandidates(const PfoList &parentCosmicRayPfos, const PfoToSliceIdMap &pfoToSliceIdMap, CRCandidateList &candidates) const
{
    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        if (!LArPfoHelper::IsFinalState(pPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        candidates.push_back(CRCandidate(this->GetPandora(), pPfo, pfoToSliceIdMap.at(pPfo)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::CheckIfOutOfTime(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsCosmicRayMap) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());

    for (const CRCandidate &candidate : candidates)
    {
        // Move on if this candidate has already been tagged
        if (pfoToIsCosmicRayMap.at(candidate.m_pPfo))
            continue;
        // Cosmic-ray muons extending outside of (single) physical volume if given t0 is that of the beam particle
        float minX(std::numeric_limits<float>::max()), maxX(-std::numeric_limits<float>::max());

        if (candidate.m_canFit)
        {
            minX = ((candidate.m_endPoint1.GetX() < candidate.m_endPoint2.GetX()) ? candidate.m_endPoint1.GetX() : candidate.m_endPoint2.GetX());
            maxX = ((candidate.m_endPoint1.GetX() > candidate.m_endPoint2.GetX()) ? candidate.m_endPoint1.GetX() : candidate.m_endPoint2.GetX());
        }
        else
        {
            // Handle any particles with small numbers of 3D hits, for which no 3D sliding fit information is available
            for (const Cluster *const pCluster : candidate.m_pPfo->GetClusterList())
            {
                float clusterMinX(std::numeric_limits<float>::max()), clusterMaxX(-std::numeric_limits<float>::max());
                pCluster->GetClusterSpanX(clusterMinX, clusterMaxX);
                minX = std::min(clusterMinX, minX);
                maxX = std::max(clusterMaxX, maxX);
            }
        }

        bool isInTime((minX > m_face_Xa - m_inTimeMargin) && (maxX < m_face_Xc + m_inTimeMargin));

        // Cosmic-ray muons that have been shifted and stitched across mid plane between volumes
        if (isInTime)
        {
            try
            {
                if (std::fabs(LArPfoHelper::GetVertex(candidate.m_pPfo)->GetX0()) > m_inTimeMaxX0)
                    isInTime = false;
            }
            catch (const StatusCodeException &)
            {
            }
        }

        // Cosmic-ray muons extending outside of (any individual) physical volume if given t0 is that of the beam particle
        if (isInTime)
        {
            CaloHitList caloHitList;
            LArPfoHelper::GetCaloHits(candidate.m_pPfo, TPC_VIEW_U, caloHitList);
            LArPfoHelper::GetCaloHits(candidate.m_pPfo, TPC_VIEW_V, caloHitList);
            LArPfoHelper::GetCaloHits(candidate.m_pPfo, TPC_VIEW_W, caloHitList);

            bool isFirstHit(true);
            bool isInSingleVolume(true);
            unsigned int volumeId(std::numeric_limits<unsigned int>::max());

            for (const CaloHit *const pCaloHit : caloHitList)
            {
                const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));

                if (!pLArCaloHit)
                    continue;

                if (isFirstHit)
                {
                    isFirstHit = false;
                    volumeId = pLArCaloHit->GetLArTPCVolumeId();
                }
                else if (volumeId != pLArCaloHit->GetLArTPCVolumeId())
                {
                    isInSingleVolume = false;
                    break;
                }
            }

            LArTPCMap::const_iterator tpcIter(larTPCMap.find(volumeId));

            if (isInSingleVolume && (larTPCMap.end() != tpcIter))
            {
                const float thisFaceXLow(tpcIter->second->GetCenterX() - 0.5f * tpcIter->second->GetWidthX());
                const float thisFaceXHigh(tpcIter->second->GetCenterX() + 0.5f * tpcIter->second->GetWidthX());

                if (!((minX > thisFaceXLow - m_inTimeMargin) && (maxX < thisFaceXHigh + m_inTimeMargin)))
                    isInTime = false;
            }
        }

        pfoToIsCosmicRayMap[candidate.m_pPfo] = !isInTime;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::CheckIfTopToBottom(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsCosmicRayMap) const
{
    for (const CRCandidate &candidate : candidates)
    {
        // Move on if this candidate has already been tagged
        if (pfoToIsCosmicRayMap.at(candidate.m_pPfo))
            continue;
        const float upperY(
            (candidate.m_endPoint1.GetY() > candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetY() : candidate.m_endPoint2.GetY());
        const float lowerY(
            (candidate.m_endPoint1.GetY() < candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetY() : candidate.m_endPoint2.GetY());

        pfoToIsCosmicRayMap[candidate.m_pPfo] = ((upperY > m_face_Yt - m_marginY) && (lowerY < m_face_Yb + m_marginY));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::CheckIfTopEntering(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsCosmicRayMap) const
{
    for (const CRCandidate &candidate : candidates)
    {
        // Move on if this candidate has already been tagged
        if (pfoToIsCosmicRayMap.at(candidate.m_pPfo))
            continue;
        const float upperY(
            (candidate.m_endPoint1.GetY() > candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetY() : candidate.m_endPoint2.GetY());
        pfoToIsCosmicRayMap[candidate.m_pPfo] = (upperY > m_face_Yt - m_marginY);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::CheckIfInVetoedTPC(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsCosmicRayMap) const
{
    if (m_vetoedTPCs.empty())
        return;

    for (const CRCandidate &candidate : candidates)
    {
        // Move on if this candidate has already been tagged
        if (pfoToIsCosmicRayMap.at(candidate.m_pPfo))
            continue;

        // We have to go via the hit volume id to get the TPC
        CaloHitList caloHitList3D;
        LArPfoHelper::GetCaloHits(candidate.m_pPfo, TPC_3D, caloHitList3D);
        
        float maxYPosition{-1.f * std::numeric_limits<float>::max()};
        float maxYXPosition{-1.f * std::numeric_limits<float>::max()};
        unsigned int maxYTPCId{std::numeric_limits<unsigned int>::max()};
        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            const float hitYPos(pCaloHit->GetPositionVector().GetY());
            const float hitXPos(pCaloHit->GetPositionVector().GetX());
            if (maxYPosition < hitYPos)
            {
                maxYPosition = hitYPos;
                maxYXPosition = hitXPos;
                const unsigned int maxYDriftVol = reinterpret_cast<LArCaloHit*>(const_cast<void*>(pCaloHit->GetParentAddress()))->GetLArTPCVolumeId();
                const unsigned int maxYAPA = reinterpret_cast<LArCaloHit*>(const_cast<void*>(pCaloHit->GetParentAddress()))->GetDaughterVolumeId();
                maxYTPCId = maxYDriftVol + maxYAPA * 4;
            }
        }
        std::cout << "Pfo " << candidate.m_pPfo << " starts in TPC " << maxYTPCId << " X " << maxYXPosition << std::endl;
        if (std::find(m_vetoedTPCs.begin(), m_vetoedTPCs.end(), maxYTPCId) != m_vetoedTPCs.end())
            pfoToIsCosmicRayMap[candidate.m_pPfo] = true;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TestBeamCosmicRayTaggingTool::CRCandidate::CRCandidate(const Pandora &pandora, const ParticleFlowObject *const pPfo, const unsigned int sliceId) :
    m_pPfo(pPfo),
    m_sliceId(sliceId),
    m_canFit(false),
    m_endPoint1(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
    m_endPoint2(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
    m_length(std::numeric_limits<float>::max()),
    m_curvature(std::numeric_limits<float>::max()),
    m_theta(std::numeric_limits<float>::max())
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

    if (!clusters3D.empty() && (clusters3D.front()->GetNCaloHits() > 15)) // TODO Configurable
    {
        m_canFit = true;
        const LArTPC *const pFirstLArTPC(pandora.GetGeometry()->GetLArTPCMap().begin()->second);
        const ThreeDSlidingFitResult slidingFitResult(clusters3D.front(), 5, pFirstLArTPC->GetWirePitchW()); // TODO Configurable
        this->CalculateFitVariables(slidingFitResult);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestBeamCosmicRayTaggingTool::CRCandidate::CalculateFitVariables(const ThreeDSlidingFitResult &slidingFitResult)
{
    m_endPoint1 = slidingFitResult.GetGlobalMinLayerPosition();
    m_endPoint2 = slidingFitResult.GetGlobalMaxLayerPosition();
    m_length = (m_endPoint2 - m_endPoint1).GetMagnitude();

    if (std::fabs(m_length) > std::numeric_limits<float>::epsilon())
        m_theta = std::fabs(m_endPoint2.GetY() - m_endPoint1.GetY()) / m_length;

    const float layerPitch(slidingFitResult.GetFirstFitResult().GetLayerPitch());

    CartesianPointVector directionList;
    for (int i = slidingFitResult.GetMinLayer(); i < slidingFitResult.GetMaxLayer(); ++i)
    {
        CartesianVector direction(0.f, 0.f, 0.f);
        if (STATUS_CODE_SUCCESS == slidingFitResult.GetGlobalFitDirection(static_cast<float>(i) * layerPitch, direction))
            directionList.push_back(direction);
    }

    CartesianVector meanDirection(0.f, 0.f, 0.f);
    for (const CartesianVector &direction : directionList)
        meanDirection += direction;

    if (!directionList.empty() > 0)
        meanDirection *= 1.f / static_cast<float>(directionList.size());

    m_curvature = 0.f;
    for (const CartesianVector &direction : directionList)
        m_curvature += (direction - meanDirection).GetMagnitude();

    if (!directionList.empty() > 0)
        m_curvature *= 1.f / static_cast<float>(directionList.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TestBeamCosmicRayTaggingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AngularUncertainty", m_angularUncertainty));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PositionalUncertainty", m_positionalUncertainty));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAssociationDist", m_maxAssociationDist));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitThreshold", m_minimumHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InTimeMargin", m_inTimeMargin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InTimeMaxX0", m_inTimeMaxX0));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MarginY", m_marginY));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TagTopEntering", m_tagTopEntering));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TagTopToBottom", m_tagTopToBottom));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TagOutOfTime", m_tagOutOfTime));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TagInVetoedTPCs", m_tagInVetoedTPCs));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "VetoedTPCs", m_vetoedTPCs));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
