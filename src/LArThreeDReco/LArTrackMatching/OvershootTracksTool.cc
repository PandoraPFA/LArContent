/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/OvershootTracksTool.cc
 * 
 *  @brief  Implementation of the overshoot tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArTrackMatching/OvershootTracksTool.h"

using namespace pandora;

namespace lar
{

OvershootTracksTool::OvershootTracksTool() :
    ThreeDKinkBaseTool(1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootTracksTool::GetIteratorListModifications(const IteratorList &iteratorList, ModificationList &modificationList) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = iIter1; iIter2 != iIter1End; ++iIter2)
        {
            if (iIter1 == iIter2)
                continue;

            try
            {
                const unsigned int nMatchedSamplingPoints1((*iIter1)->GetOverlapResult().GetNMatchedSamplingPoints());
                const unsigned int nMatchedSamplingPoints2((*iIter2)->GetOverlapResult().GetNMatchedSamplingPoints());
                IteratorList::const_iterator iIterA((nMatchedSamplingPoints1 >= nMatchedSamplingPoints2) ? iIter1 : iIter2);
                IteratorList::const_iterator iIterB((nMatchedSamplingPoints1 >= nMatchedSamplingPoints2) ? iIter2 : iIter1);

                SplitParticle splitParticle(*(*iIterA), *(*iIterB));
                const LArPointingCluster pointingClusterA1(splitParticle.m_pClusterA1);
                const LArPointingCluster pointingClusterB1(splitParticle.m_pClusterB1);
                const LArPointingCluster pointingClusterA2(splitParticle.m_pClusterA2);
                const LArPointingCluster pointingClusterB2(splitParticle.m_pClusterB2);

                LArPointingCluster::Vertex vertexA1, vertexB1, vertexA2, vertexB2;
                LArPointingClusterHelper::GetClosestVerticesInX(pointingClusterA1, pointingClusterB1, vertexA1, vertexB1);
                LArPointingClusterHelper::GetClosestVerticesInX(pointingClusterA2, pointingClusterB2, vertexA2, vertexB2);

                if (!this->PassesVertexCuts(vertexA1, vertexB1) || !this->PassesVertexCuts(vertexA2, vertexB2))
                    continue;

                this->SetSplitPosition(vertexA1, vertexA2, vertexB1, vertexB2, splitParticle);

                Modification modification;
                modification.m_splitPositionMap[splitParticle.m_pCommonCluster].push_back(splitParticle.m_splitPosition);

                modification.m_affectedClusters.insert(splitParticle.m_pCommonCluster);
                modification.m_affectedClusters.insert(splitParticle.m_pClusterA1);
                modification.m_affectedClusters.insert(splitParticle.m_pClusterA2);
                modification.m_affectedClusters.insert(splitParticle.m_pClusterB1);
                modification.m_affectedClusters.insert(splitParticle.m_pClusterB2);

                modificationList.push_back(modification);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;

                continue;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OvershootTracksTool::PassesVertexCuts(const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB) const
{
    float longitudinalAB(-std::numeric_limits<float>::max()), transverseAB(std::numeric_limits<float>::max());
    LArPointingClusterHelper::GetImpactParameters(vertexA, vertexB, longitudinalAB, transverseAB);

    float longitudinalBA(-std::numeric_limits<float>::max()), transverseBA(std::numeric_limits<float>::max());
    LArPointingClusterHelper::GetImpactParameters(vertexB, vertexA, longitudinalBA, transverseBA);

    if (std::min(longitudinalAB, longitudinalBA) < m_minLongitudinalImpactParameter)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootTracksTool::SetSplitPosition(const LArPointingCluster::Vertex &vertexA1, const LArPointingCluster::Vertex &vertexA2,
    const LArPointingCluster::Vertex &vertexB1, const LArPointingCluster::Vertex &vertexB2, SplitParticle &splitParticle) const
{
    bool splitAtElementA(false), splitAtElementB(false);

    if (std::fabs(vertexA1.GetPosition().GetX() - vertexA2.GetPosition().GetX()) < m_maxVertexXSeparation)
    {
        splitAtElementA = true;
    }
    else if (std::fabs(vertexB1.GetPosition().GetX() - vertexB2.GetPosition().GetX()) < m_maxVertexXSeparation)
    {
        splitAtElementB = true;
    }

    if (!splitAtElementA && !splitAtElementB)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianVector splitPosition(0.f, 0.f, 0.f);
    float chiSquared(std::numeric_limits<float>::max());

    LArGeometryHelper::MergeTwoPositions(LArThreeDHelper::GetClusterHitType(splitParticle.m_pClusterA1),
        LArThreeDHelper::GetClusterHitType(splitParticle.m_pClusterA2), splitAtElementA ? vertexA1.GetPosition() : vertexB1.GetPosition(),
        splitAtElementA ? vertexA2.GetPosition() : vertexB2.GetPosition(), splitPosition, chiSquared);

    splitParticle.m_splitPosition = splitPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

OvershootTracksTool::SplitParticle::SplitParticle(const TensorType::Element &elementA, const TensorType::Element &elementB) :
    m_splitPosition(0.f, 0.f, 0.f)
{
    const HitType commonView((elementA.GetClusterU() == elementB.GetClusterU()) ? VIEW_U :
        (elementA.GetClusterV() == elementB.GetClusterV()) ? VIEW_V :
        (elementA.GetClusterW() == elementB.GetClusterW()) ? VIEW_W :
        throw StatusCodeException(STATUS_CODE_FAILURE));

    m_pCommonCluster = (VIEW_U == commonView) ? elementA.GetClusterU() : (VIEW_V == commonView) ? elementA.GetClusterV() : elementA.GetClusterW();
    m_pClusterA1 = (VIEW_U == commonView) ? elementA.GetClusterV() : (VIEW_V == commonView) ? elementA.GetClusterU() : elementA.GetClusterU();
    m_pClusterA2 = (VIEW_U == commonView) ? elementA.GetClusterW() : (VIEW_V == commonView) ? elementA.GetClusterW() : elementA.GetClusterV();
    m_pClusterB1 = (VIEW_U == commonView) ? elementB.GetClusterV() : (VIEW_V == commonView) ? elementB.GetClusterU() : elementB.GetClusterU();
    m_pClusterB2 = (VIEW_U == commonView) ? elementB.GetClusterW() : (VIEW_V == commonView) ? elementB.GetClusterW() : elementB.GetClusterV();

    if ((m_pClusterA1 == m_pClusterB1) || (m_pClusterA2 == m_pClusterB2))
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OvershootTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_maxVertexXSeparation = 2.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexXSeparation", m_maxVertexXSeparation));

    return ThreeDKinkBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar
