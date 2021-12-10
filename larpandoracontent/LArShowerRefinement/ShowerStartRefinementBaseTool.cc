/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.cc
 *
 *  @brief  Implementation of the shower start refinement base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

using namespace pandora;

namespace lar_content
{

ShowerStartRefinementBaseTool::ShowerStartRefinementBaseTool() : 
    m_maxDistanceForConnection(5.f),
    m_growingFitInitialLength(10.f),
    m_macroSlidingFitWindow(1000),
    m_growingFitSegmentLength(5.0f),
    m_distanceToLine(1.f),
    m_initialFitDistanceToLine(0.5f),
    m_maxFittingHits(20),
    m_longitudinalCoordinateBinSize(1.f),
    m_hitConnectionDistance(1.f),
    m_minInitialHitsFound(10)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementBaseTool::HasPathToNuVertex(const ParticleFlowObject *const pShowerPfo, const CartesianVector &neutrinoVertex) const
{
    // this algorithm has to run after we have 3D hits but not after shower vertex creation
    // therefore, find the closest space point to the neutrino vertex

    ClusterList clusterList3D;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, clusterList3D);

    for (const Cluster *const pCluster3D : clusterList3D)
    {
        CaloHitList caloHitList3D;
        pCluster3D->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

        for (const CaloHit *const caloHit3D : caloHitList3D)
        {
            const CartesianVector &hitPosition(caloHit3D->GetPositionVector());
            const float separationSquared((neutrinoVertex - hitPosition).GetMagnitudeSquared());

            if (separationSquared < (m_maxDistanceForConnection * m_maxDistanceForConnection))
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::FindShowerSpine(const ShowerStartRefinementAlgorithm *pAlgorithm, const CaloHitList &viewShowerHitList, 
    const CartesianVector &projectedNuVertexPosition, const CartesianVector &initialDirection, CaloHitList &unavailableHitList, CaloHitList &showerSpineHitList)
{
    // Construct initial fit
    float highestL(0.f);
    CartesianPointVector runningFitPositionVector;

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector &displacementVector(hitPosition - projectedNuVertexPosition);
        const float l(initialDirection.GetDotProduct(displacementVector));
        const float t(initialDirection.GetCrossProduct(displacementVector).GetMagnitude());

        if ((l < m_growingFitInitialLength) && (l > 0.f) && (t < m_initialFitDistanceToLine))
        {
            if (l > highestL)
                highestL = l;

            if (std::find(unavailableHitList.begin(), unavailableHitList.end(), pCaloHit) == unavailableHitList.end())
            {
                showerSpineHitList.push_back(pCaloHit);
                unavailableHitList.push_back(pCaloHit);
            }

            runningFitPositionVector.push_back(hitPosition);
        }
    }

    // Require significant number of initial hits
    if (runningFitPositionVector.size() < m_minInitialHitsFound)
    {
        std::cout << "Not enough initial hits to form fit" << std::endl;
        showerSpineHitList.clear();
        return;
    }

    //////////////////////////////////
    /*
    for (const CartesianVector position : runningFitPositionVector)
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &position, "initial hit", VIOLET, 2);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////////////////////

    // Perform a running fit to follow pathway
    const bool isEndDownstream(initialDirection.GetZ() > 0.f);
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    CartesianVector extrapolatedStartPosition(projectedNuVertexPosition), extrapolatedDirection(initialDirection);
    CartesianVector extrapolatedEndPosition(extrapolatedStartPosition + (extrapolatedDirection * highestL));

    unsigned int count(0);
    bool hitsCollected(true);

    // make this a boolean
    while (hitsCollected)
    {
        ++count;

        try
        {
            const int excessHitsInFit(runningFitPositionVector.size() - m_maxFittingHits);

            if (excessHitsInFit > 0)
            {
                // Remove furthest away hits
                std::sort(runningFitPositionVector.begin(), runningFitPositionVector.end(), SortByDistanceToPoint(extrapolatedEndPosition));

                for (int i = 0; i < excessHitsInFit; ++i)
                    runningFitPositionVector.erase(std::prev(runningFitPositionVector.end()));
            }

            const TwoDSlidingFitResult extrapolatedFit(&runningFitPositionVector, m_macroSlidingFitWindow, slidingFitPitch);            

            extrapolatedStartPosition = count == 1 ? extrapolatedEndPosition : isEndDownstream ? extrapolatedFit.GetGlobalMaxLayerPosition() : extrapolatedFit.GetGlobalMinLayerPosition();
            extrapolatedDirection = isEndDownstream ? extrapolatedFit.GetGlobalMaxLayerDirection() : extrapolatedFit.GetGlobalMinLayerDirection() * (-1.f);
            extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);

            hitsCollected = this->CollectSubsectionHits(pAlgorithm, extrapolatedFit, extrapolatedStartPosition, extrapolatedEndPosition, extrapolatedDirection,
                isEndDownstream, viewShowerHitList, runningFitPositionVector, unavailableHitList, showerSpineHitList);

            // As a final effort, reduce the sliding fit window
            if (!hitsCollected)
            {
                const TwoDSlidingFitResult microExtrapolatedFit(&runningFitPositionVector, 5.f, slidingFitPitch);
            
                extrapolatedStartPosition = count == 1 ? extrapolatedStartPosition : isEndDownstream ? microExtrapolatedFit.GetGlobalMaxLayerPosition() : microExtrapolatedFit.GetGlobalMinLayerPosition();
                extrapolatedDirection = isEndDownstream ? microExtrapolatedFit.GetGlobalMaxLayerDirection() : microExtrapolatedFit.GetGlobalMinLayerDirection() * (-1.f);
                extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);

                hitsCollected = this->CollectSubsectionHits(pAlgorithm, extrapolatedFit, extrapolatedStartPosition, extrapolatedEndPosition, extrapolatedDirection,
                    isEndDownstream, viewShowerHitList, runningFitPositionVector, unavailableHitList, showerSpineHitList);
            }
        }
        catch (const StatusCodeException &)
        {
            std::cout << "couldn't fit runningFitPositionVector" << std::endl;
            return;
        }
    }

    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementBaseTool::CollectSubsectionHits(const ShowerStartRefinementAlgorithm *pAlgorithm, const TwoDSlidingFitResult &extrapolatedFit, 
    const CartesianVector &extrapolatedStartPosition, const CartesianVector &extrapolatedEndPosition, const CartesianVector &extrapolatedDirection, 
    const bool isEndDownstream, const CaloHitList &viewShowerHitList, CartesianPointVector &runningFitPositionVector, CaloHitList &unavailableHitList, 
    CaloHitList &showerSpineHitList)
{ 
    ////////////////////////////
    //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &extrapolatedStartPosition, "start", GREEN, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &extrapolatedEndPosition, "end", GREEN, 2);
    ////////////////////////////

    float extrapolatedStartL(0.f), extrapolatedStartT(0.f);
    extrapolatedFit.GetLocalPosition(extrapolatedStartPosition, extrapolatedStartL, extrapolatedStartT);

    float extrapolatedEndL(0.f), extrapolatedEndT(0.f);
    extrapolatedFit.GetLocalPosition(extrapolatedEndPosition, extrapolatedEndL, extrapolatedEndT);

    CaloHitList collectedHits;

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        if (std::find(showerSpineHitList.begin(), showerSpineHitList.end(), pCaloHit) != showerSpineHitList.end())
            continue;

        if (std::find(unavailableHitList.begin(), unavailableHitList.end(), pCaloHit) != unavailableHitList.end())
            continue;

        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        float hitL(0.f), hitT(0.f);
        extrapolatedFit.GetLocalPosition(hitPosition, hitL, hitT);

        // Assess whether hit is within section boundaries
        if (isEndDownstream && ((hitL < extrapolatedStartL) || (hitL > extrapolatedEndL)))
            continue;

        if (!isEndDownstream && ((hitL > extrapolatedStartL) || (hitL < extrapolatedEndL)))
            continue;

        // Assess whether hit is close to connecting line - taking account hit width if necessary
        if (this->IsCloseToLine(hitPosition, extrapolatedStartPosition, extrapolatedDirection, m_distanceToLine))
        {
            collectedHits.push_back(pCaloHit);
        }
        else
        {
            const CartesianVector closestPointInHit(LArHitWidthHelper::GetClosestPointToLine2D(extrapolatedStartPosition, extrapolatedDirection, pCaloHit));

            if (this->IsCloseToLine(closestPointInHit, extrapolatedStartPosition, extrapolatedDirection, m_distanceToLine))
                collectedHits.push_back(pCaloHit);
        }
    }

    const int nInitialHits(showerSpineHitList.size());
    this->CollectConnectedHits(pAlgorithm, collectedHits, extrapolatedStartPosition, extrapolatedDirection, runningFitPositionVector, unavailableHitList, showerSpineHitList);
    const int nFinalHits(showerSpineHitList.size());

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return (nFinalHits != nInitialHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::CollectConnectedHits(const ShowerStartRefinementAlgorithm *pAlgorithm, const CaloHitList &collectedHits, 
    const CartesianVector &extrapolatedStartPosition, const CartesianVector &extrapolatedDirection, CartesianPointVector &runningFitPositionVector, 
    CaloHitList &unavailableHitList, CaloHitList &showerSpineHitList)
{
    // Now add connected hits - taking into account hit width if necessary
    bool found = true;

    while (found)
    {
        found = false;

        for (const CaloHit *const pCaloHit : collectedHits)
        {
            if (std::find(showerSpineHitList.begin(), showerSpineHitList.end(), pCaloHit) != showerSpineHitList.end())
                continue;

            CartesianVector hitPosition(pCaloHit->GetPositionVector());

            if (this->GetClosestDistance(hitPosition, runningFitPositionVector) > m_hitConnectionDistance)
            {
                hitPosition = LArHitWidthHelper::GetClosestPointToLine2D(extrapolatedStartPosition, extrapolatedDirection, pCaloHit);
                        
                if (LArHitWidthHelper::GetClosestDistance(pCaloHit, showerSpineHitList) > m_hitConnectionDistance)
                    continue;
            }

            found = true;

            runningFitPositionVector.push_back(hitPosition);
            showerSpineHitList.push_back(pCaloHit);
            unavailableHitList.push_back(pCaloHit);
            ////////////////////////////////
            //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "added hit", BLUE, 2);
            ////////////////////////////////
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::GetHitsInBoundingBox(const CartesianVector &firstCorner, const CartesianVector &secondCorner, const CaloHitList &inputHitList,
    const float distanceToLine, CaloHitList &hitsInBoundingBox) const 
{
    const float minX(std::min(firstCorner.GetX(), secondCorner.GetX())), maxX(std::max(firstCorner.GetX(), secondCorner.GetX()));
    const float minZ(std::min(firstCorner.GetZ(), secondCorner.GetZ())), maxZ(std::max(firstCorner.GetZ(), secondCorner.GetZ()));

    CartesianVector connectingLineDirection(firstCorner - secondCorner);
    connectingLineDirection = connectingLineDirection.GetUnitVector();

    for (const CaloHit *const pCaloHit : inputHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if(!this->IsInBoundingBox(minX, maxX, minZ, maxZ, hitPosition))
            continue;

        if (distanceToLine > 0.f)
        {
            if (!this->IsCloseToLine(hitPosition, firstCorner, connectingLineDirection, distanceToLine))
                continue;
        }

        hitsInBoundingBox.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementBaseTool::IsInBoundingBox(const float minX, const float maxX, const float minZ, const float maxZ, 
    const CartesianVector &hitPosition) const
{
    if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < minZ) || (hitPosition.GetZ() > maxZ))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementBaseTool::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart, 
    const CartesianVector &lineDirection, const float distanceToLine) const
{
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());
    
    if (transverseDistanceFromLine > distanceToLine)
        return false;

    return true;
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementBaseTool::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
{
    const float segmentBoundaryGradient = (-1.f) * (upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ()) / segmentBoundaryGradient + upperBoundary.GetX());
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ()) / segmentBoundaryGradient + lowerBoundary.GetX());
    
    if (std::fabs(xPointOnUpperLine - point.GetX()) < std::numeric_limits<float>::epsilon())
        return true;

    if (std::fabs(xPointOnLowerLine - point.GetX()) < std::numeric_limits<float>::epsilon())
        return true;
    
    if ((point.GetX() > xPointOnUpperLine) && (point.GetX() > xPointOnLowerLine))
        return false;

    if ((point.GetX() < xPointOnUpperLine) && (point.GetX() < xPointOnLowerLine))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::BuildProtoShowers(const ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const
{
    // cheat this to find all shower starts of significant MC particles in the pfo
    this->FindShowerCores(pShowerPfo, protoShowerVector);

    if (protoShowerVector.empty())
    {
        std::cout << "NO SHOWER CORES FOUND" << std::endl;
        return;
    }

    // bail if it can't be done? or take out of vector?
    this->FindShowerStartPositions(pShowerPfo, protoShowerVector);

    this->FindConnectionPathways(pShowerPfo, protoShowerVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShowerStartRefinementBaseTool::GetClosestDistance(const CartesianVector &position, const CartesianPointVector &testPositions) const
{
    float closestDistanceSqaured(std::numeric_limits<float>::max());

    for (const CartesianVector &testPosition : testPositions)
    {
        const float separationSquared((testPosition - position).GetMagnitudeSquared());

        if (separationSquared < closestDistanceSqaured)
            closestDistanceSqaured = separationSquared;
    }

    return std::sqrt(closestDistanceSqaured);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::FindShowerCores(const ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const
{
    // do fancy logic
    std::cout << pShowerPfo << protoShowerVector.size() << std::endl;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::FindShowerStartPositions(const ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const
{
    // do fancy logic
    std::cout << pShowerPfo << protoShowerVector.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::FindConnectionPathways(const ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const
{
    // do fancy logic
    std::cout << pShowerPfo << protoShowerVector.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementBaseTool::IsElectronPathway(const ProtoShower &protoShower)
{
    std::cout << protoShower.m_showerCore.m_startPosition.GetX() << std::endl;
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    bool SortByDistanceToPoint::operator()(const CartesianVector &lhs, const CartesianVector &rhs)
    {
        return (m_referencePoint.GetDistanceSquared(lhs) < m_referencePoint.GetDistanceSquared(rhs));
    }

//------------------------------------------------------------------------------------------------------------------------------------------


StatusCode ShowerStartRefinementBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDistanceForConnection", m_maxDistanceForConnection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "GrowingFitInitialLength", m_growingFitInitialLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MacroSlidingFitWindow", m_macroSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "GrowingFitSegmentLength", m_growingFitSegmentLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "DistanceToLine", m_distanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "InitialFitDistanceToLine", m_initialFitDistanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxFittingsHits", m_maxFittingHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "LongitudinalCoordinateBinSize", m_longitudinalCoordinateBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitConnectionDistance", m_hitConnectionDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitInitialHitsFound", m_minInitialHitsFound));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
