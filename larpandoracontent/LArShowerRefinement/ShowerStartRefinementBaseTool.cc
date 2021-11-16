/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.cc
 *
 *  @brief  Implementation of the shower start refinement base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

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
    m_distanceToLine(1.0f),
    m_maxFittingHits(20)
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

void ShowerStartRefinementBaseTool::FindShowerSpine(const ShowerStartRefinementAlgorithm *pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &neutrinoVertex, 
    const CartesianVector &initialDirection, const HitType hitType, CaloHitList &showerSpineHitList, LongitudinalPositionMap &longitudinalPositionMap)
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, caloHitList);

    // Construct initial fit
    const bool isEndDownstream(initialDirection.GetZ() > 0.f);
    CartesianVector extrapolatedStartPosition(neutrinoVertex), extrapolatedDirection(initialDirection);
    const CartesianVector clusterSubsetBoundary(extrapolatedStartPosition + (extrapolatedDirection * m_growingFitInitialLength));

    CaloHitList subsetHitList;
    this->GetHitsInBoundingBox(extrapolatedStartPosition, clusterSubsetBoundary, caloHitList, m_distanceToLine, subsetHitList);

    if (subsetHitList.empty())
        return;

    showerSpineHitList.insert(showerSpineHitList.end(), subsetHitList.begin(), subsetHitList.end());
    
    CartesianPointVector runningFitPositionVector;
    for (const CaloHit *const pCaloHit : subsetHitList)
    {
        std::cout << "extrapolatedDirection: " << extrapolatedDirection << std::endl;

        longitudinalPositionMap[pCaloHit] = extrapolatedDirection.GetDotProduct(pCaloHit->GetPositionVector() - extrapolatedStartPosition);
        runningFitPositionVector.push_back(pCaloHit->GetPositionVector());
    }

    float longitudinalDistance(m_growingFitInitialLength);

    // Collect extrapolated hits by performing a running fit
    unsigned int count(0);
    CartesianVector extrapolatedEndPosition(0.f, 0.f, 0.f);
    unsigned int hitsCollected(1);
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    std::cout << "INITIAL HITS IN FIT " << std::endl;

    for (const CartesianVector position : runningFitPositionVector)
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &position, "initial hit", VIOLET, 2);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    while (hitsCollected)
    {
        hitsCollected = 0;

        try
        {
            int excess(runningFitPositionVector.size() - m_maxFittingHits);

            if (excess > 0)
            {
                std::sort(runningFitPositionVector.begin(), runningFitPositionVector.end(), SortByDistanceToPoint(extrapolatedEndPosition));

                for (int i = 0; i < excess; ++i)
                    runningFitPositionVector.erase(runningFitPositionVector.end() - 1);
            }

            std::cout << "runningFitPositionVector.size: " << runningFitPositionVector.size() << std::endl;

            const TwoDSlidingFitResult extrapolatedFit(&runningFitPositionVector, m_macroSlidingFitWindow, slidingFitPitch);
            
            extrapolatedStartPosition = isEndDownstream ? extrapolatedFit.GetGlobalMaxLayerPosition() : extrapolatedFit.GetGlobalMinLayerPosition();
            extrapolatedDirection = isEndDownstream ? extrapolatedFit.GetGlobalMaxLayerDirection() : extrapolatedFit.GetGlobalMinLayerDirection() * (-1.f);
            extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);

            for (const CaloHit *const pCaloHit : caloHitList)
            {
                if (std::find(subsetHitList.begin(), subsetHitList.end(), pCaloHit) != subsetHitList.end())
                    continue;

                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                    
                if (!this->IsInLineSegment(extrapolatedStartPosition, extrapolatedEndPosition, hitPosition))
                    continue;

                if (!this->IsCloseToLine(hitPosition, extrapolatedStartPosition, extrapolatedDirection, m_distanceToLine))
                    continue;

                ++hitsCollected;

                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &extrapolatedStartPosition, "start", GREEN, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "added hit", BLUE, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &extrapolatedEndPosition, "end", GREEN, 2);
                    
                longitudinalPositionMap[pCaloHit] = longitudinalDistance + extrapolatedDirection.GetDotProduct(hitPosition - extrapolatedStartPosition);

                subsetHitList.push_back(pCaloHit);
                showerSpineHitList.push_back(pCaloHit);
                runningFitPositionVector.push_back(hitPosition);
            }
        }
        catch (const StatusCodeException &)
        {
            return;
        }

        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        ++count;
        longitudinalDistance += m_growingFitSegmentLength;
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
        XmlHelper::ReadValue(xmlHandle, "MaxFittingsHits", m_maxFittingHits));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
