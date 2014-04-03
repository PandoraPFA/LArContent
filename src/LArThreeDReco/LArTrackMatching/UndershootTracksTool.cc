/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/UndershootTracksTool.cc
 * 
 *  @brief  Implementation of the undershoot tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArTrackMatching/UndershootTracksTool.h"

using namespace pandora;

namespace lar
{

UndershootTracksTool::UndershootTracksTool() :
    ThreeDKinkBaseTool(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UndershootTracksTool::GetIteratorListModifications(ThreeDTransverseTracksAlgorithm *pAlgorithm, const IteratorList &iteratorList,
    ModificationList &modificationList) const
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

                Particle particle(*(*iIterA), *(*iIterB));
                const LArPointingCluster pointingClusterA(particle.m_pClusterA);
                const LArPointingCluster pointingClusterB(particle.m_pClusterB);

                LArPointingCluster::Vertex vertexA, vertexB;
                LArPointingClusterHelper::GetClosestVertices(pointingClusterA, pointingClusterB, vertexA, vertexB);

                float transverseAB(std::numeric_limits<float>::max()), transverseBA(std::numeric_limits<float>::max());
                float longitudinalAB(-std::numeric_limits<float>::max()), longitudinalBA(-std::numeric_limits<float>::max());

                LArPointingClusterHelper::GetImpactParameters(vertexA, vertexB, longitudinalAB, transverseAB);
                LArPointingClusterHelper::GetImpactParameters(vertexB, vertexA, longitudinalBA, transverseBA);

                if (std::min(longitudinalAB, longitudinalBA) < m_minLongitudinalImpactParameter)
                    continue;

                if (std::min(transverseAB, transverseBA) > m_maxTransverseImpactParameter)
                    continue;

                const float cosTheta(-vertexA.GetDirection().GetCosOpeningAngle(vertexB.GetDirection()));

                if (cosTheta < m_minImpactParameterCosTheta)
                    continue;

                const bool isThreeDKink(false); // TODO;

                if (isThreeDKink != m_splitMode)
                    continue;

                // Construct the modification object
                Modification modification;

                if (m_splitMode)
                {
                    const float splitX((vertexA.GetPosition().GetX() + vertexB.GetPosition().GetX()) * 0.5f);

                    CartesianVector splitPosition1(0.f, 0.f, 0.f);
                    pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster1).GetGlobalFitPosition(splitX, true, splitPosition1);
                    modification.m_splitPositionMap[particle.m_pCommonCluster1].push_back(splitPosition1);

                    CartesianVector splitPosition2(0.f, 0.f, 0.f);
                    pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster2).GetGlobalFitPosition(splitX, true, splitPosition2);
                    modification.m_splitPositionMap[particle.m_pCommonCluster2].push_back(splitPosition2);
                }
                else
                {
                    const bool vertexAIsLowX(vertexA.GetPosition().GetX() < vertexB.GetPosition().GetX());
                    Cluster *pLowXCluster(vertexAIsLowX ? particle.m_pClusterA : particle.m_pClusterB);
                    Cluster *pHighXCluster(vertexAIsLowX ? particle.m_pClusterB : particle.m_pClusterA);
                    modification.m_clusterMergeMap[pLowXCluster].insert(pHighXCluster);
                }

                modification.m_affectedClusters.insert(particle.m_pClusterA);
                modification.m_affectedClusters.insert(particle.m_pClusterB);
                modification.m_affectedClusters.insert(particle.m_pCommonCluster1);
                modification.m_affectedClusters.insert(particle.m_pCommonCluster2);

                modificationList.push_back(modification);
            }
            catch (StatusCodeException &)
            {
                continue;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

UndershootTracksTool::Particle::Particle(const TensorType::Element &elementA, const TensorType::Element &elementB)
{
    m_pClusterA = (elementA.GetClusterU() != elementB.GetClusterU()) ? elementA.GetClusterU() :
        (elementA.GetClusterV() != elementB.GetClusterV()) ? elementA.GetClusterV() : elementA.GetClusterW();
    m_pClusterB = (elementA.GetClusterU() != elementB.GetClusterU()) ? elementB.GetClusterU() :
        (elementA.GetClusterV() != elementB.GetClusterV()) ? elementB.GetClusterV() : elementB.GetClusterW();
    m_pCommonCluster1 = (elementA.GetClusterU() == elementB.GetClusterU()) ? elementA.GetClusterU() :
        (elementA.GetClusterV() == elementB.GetClusterV()) ? elementA.GetClusterV() : elementA.GetClusterW();
    m_pCommonCluster2 = ((m_pClusterA != elementA.GetClusterU()) && (m_pCommonCluster1 != elementA.GetClusterU())) ? elementA.GetClusterU() :
        ((m_pClusterA != elementA.GetClusterV()) && (m_pCommonCluster1 != elementA.GetClusterV())) ? elementA.GetClusterV() : elementA.GetClusterW();

    if ((m_pClusterA == m_pClusterB) || (m_pCommonCluster1 == m_pCommonCluster2))
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode UndershootTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_splitMode = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SplitMode", m_splitMode));

    m_minLongitudinalImpactParameter = -1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLongitudinalImpactParameter", m_minLongitudinalImpactParameter));

    m_maxTransverseImpactParameter = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseImpactParameter", m_maxTransverseImpactParameter));

    m_minImpactParameterCosTheta = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinImpactParameterCosTheta", m_minImpactParameterCosTheta));

    return ThreeDKinkBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar
