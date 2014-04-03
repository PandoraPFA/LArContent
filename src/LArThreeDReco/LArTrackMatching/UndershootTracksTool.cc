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

void UndershootTracksTool::GetIteratorListModifications(const IteratorList &iteratorList, ModificationList &modificationList) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = iIter1; iIter2 != iIter1End; ++iIter2)
        {
            if (iIter1 == iIter2)
                continue;

            Cluster *pCluster1(((*iIter1)->GetClusterU() != (*iIter2)->GetClusterU()) ? (*iIter1)->GetClusterU() :
                ((*iIter1)->GetClusterV() != (*iIter2)->GetClusterV()) ? (*iIter1)->GetClusterV() : (*iIter1)->GetClusterW());

            Cluster *pCluster2(((*iIter1)->GetClusterU() != (*iIter2)->GetClusterU()) ? (*iIter2)->GetClusterU() :
                ((*iIter1)->GetClusterV() != (*iIter2)->GetClusterV()) ? (*iIter2)->GetClusterV() : (*iIter2)->GetClusterW());

            if (pCluster1 == pCluster2)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            try
            {
                const LArPointingCluster pointingCluster1(pCluster1);
                const LArPointingCluster pointingCluster2(pCluster2);

                LArPointingCluster::Vertex vertex1, vertex2;
                LArPointingClusterHelper::GetClosestVertices(pointingCluster1, pointingCluster2, vertex1, vertex2);

                float transverse12(std::numeric_limits<float>::max()), transverse21(std::numeric_limits<float>::max());
                float longitudinal12(-std::numeric_limits<float>::max()), longitudinal21(-std::numeric_limits<float>::max());

                LArPointingClusterHelper::GetImpactParameters(vertex1, vertex2, longitudinal12, transverse12);
                LArPointingClusterHelper::GetImpactParameters(vertex2, vertex1, longitudinal21, transverse21);

                if (std::min(longitudinal12, longitudinal21) < m_minLongitudinalImpactParameter)
                    continue;

                if (std::min(transverse12, transverse21) > m_maxTransverseImpactParameter)
                    continue;

                const float cosTheta(-vertex1.GetDirection().GetCosOpeningAngle(vertex2.GetDirection()));

                if (cosTheta < m_minImpactParameterCosTheta)
                    continue;

                const bool vertex1IsLowX(vertex1.GetPosition().GetX() < vertex2.GetPosition().GetX());
                Cluster *pLowXCluster(vertex1IsLowX ? pCluster1 : pCluster2);
                Cluster *pHighXCluster(vertex1IsLowX ? pCluster2 : pCluster1);

                Modification modification;
                modification.m_clusterMergeMap[pLowXCluster].insert(pHighXCluster);

                modification.m_affectedClusters.insert((*iIter1)->GetClusterU());
                modification.m_affectedClusters.insert((*iIter1)->GetClusterV());
                modification.m_affectedClusters.insert((*iIter1)->GetClusterW());
                modification.m_affectedClusters.insert((*iIter2)->GetClusterU());
                modification.m_affectedClusters.insert((*iIter2)->GetClusterV());
                modification.m_affectedClusters.insert((*iIter2)->GetClusterW());

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

StatusCode UndershootTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
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
