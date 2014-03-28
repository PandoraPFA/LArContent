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

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArTrackMatching/OvershootTracksTool.h"

using namespace pandora;

namespace lar
{

StatusCode OvershootTracksTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindOvershootTracks(overlapTensor, protoParticleVector);
    pAlgorithm->CreateThreeDParticles(protoParticleVector);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootTracksTool::FindOvershootTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        if (nU * nV * nW < 2)
            continue;

        std::sort(elementList.begin(), elementList.end(), ThreeDTransverseTracksAlgorithm::SortByNMatchedSamplingPoints);

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (!this->PassesElementCuts(eIter, usedClusters))
                continue;

            IteratorList iteratorList;
            this->SelectOvershootElements(eIter, elementList, usedClusters, iteratorList);

            if (iteratorList.size() < 2)
                continue;

            ProtoParticleVector localProtoParticleVector;
            this->BuildProtoParticle(iteratorList, localProtoParticleVector);

            if (localProtoParticleVector.empty())
                continue;

            protoParticleVector.insert(protoParticleVector.end(), localProtoParticleVector.begin(), localProtoParticleVector.end());

            for (ProtoParticleVector::const_iterator pIter = localProtoParticleVector.begin(); pIter != localProtoParticleVector.end(); ++pIter)
            {
                usedClusters.insert(pIter->m_clusterListU.begin(), pIter->m_clusterListU.end());
                usedClusters.insert(pIter->m_clusterListV.begin(), pIter->m_clusterListV.end());
                usedClusters.insert(pIter->m_clusterListW.begin(), pIter->m_clusterListW.end());
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootTracksTool::SelectOvershootElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
    const ClusterList &usedClusters, IteratorList &iteratorList) const
{
    iteratorList.push_back(eIter);

    for (TensorType::ElementList::const_iterator eIter2 = elementList.begin(); eIter2 != elementList.end(); ++eIter2)
    {
        if (eIter == eIter2)
            continue;

        if (!this->PassesElementCuts(eIter2, usedClusters))
            continue;

        for (IteratorList::const_iterator iIter = iteratorList.begin(); iIter != iteratorList.end(); ++iIter)
        {
            if ((*iIter) == eIter2)
                continue;

            unsigned int nMatchedClusters(0);

            if ((*iIter)->GetClusterU() == eIter2->GetClusterU())
                ++nMatchedClusters;

            if ((*iIter)->GetClusterV() == eIter2->GetClusterV())
                ++nMatchedClusters;

            if ((*iIter)->GetClusterW() == eIter2->GetClusterW())
                ++nMatchedClusters;

            if (1 == nMatchedClusters)
            {
                iteratorList.push_back(eIter2);
                return;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootTracksTool::BuildProtoParticle(const IteratorList &iteratorList, ProtoParticleVector &protoParticleVector) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = iteratorList.begin(), iIter2End = iteratorList.end(); iIter2 != iIter2End; ++iIter2)
        {
            if (iIter1 == iIter2)
                continue;

            const unsigned int nMatchedSamplingPoints1((*iIter1)->GetOverlapResult().GetNMatchedSamplingPoints());
            const unsigned int nMatchedSamplingPoints2((*iIter2)->GetOverlapResult().GetNMatchedSamplingPoints());
            IteratorList::const_iterator iIterA((nMatchedSamplingPoints1 >= nMatchedSamplingPoints2) ? iIter1 : iIter2);
            IteratorList::const_iterator iIterB((nMatchedSamplingPoints1 >= nMatchedSamplingPoints2) ? iIter2 : iIter1);

            const HitType commonView(((*iIterA)->GetClusterU() == (*iIterB)->GetClusterU()) ? VIEW_U :
                 ((*iIterA)->GetClusterV() == (*iIterB)->GetClusterV()) ? VIEW_V : VIEW_W);

            if ((VIEW_W == commonView) && ((*iIterA)->GetClusterW() != (*iIterB)->GetClusterW()))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            HitType hitType1(CUSTOM), hitType2(CUSTOM);
            Cluster *pCommonCluster(NULL), *pClusterA1(NULL), *pClusterA2(NULL), *pClusterB1(NULL), *pClusterB2(NULL);

            if (VIEW_U == commonView)
            {
                pCommonCluster = (*iIterA)->GetClusterU();
                pClusterA1 = (*iIterA)->GetClusterV();
                pClusterA2 = (*iIterA)->GetClusterW();
                pClusterB1 = (*iIterB)->GetClusterV();
                pClusterB2 = (*iIterB)->GetClusterW();
                hitType1 = VIEW_V;
                hitType2 = VIEW_W;
            }
            else if (VIEW_V == commonView)
            {
                pCommonCluster = (*iIterA)->GetClusterV();
                pClusterA1 = (*iIterA)->GetClusterU();
                pClusterA2 = (*iIterA)->GetClusterW();
                pClusterB1 = (*iIterB)->GetClusterU();
                pClusterB2 = (*iIterB)->GetClusterW();
                hitType1 = VIEW_U;
                hitType2 = VIEW_W;
            }
            else
            {
                pCommonCluster = (*iIterA)->GetClusterW();
                pClusterA1 = (*iIterA)->GetClusterU();
                pClusterA2 = (*iIterA)->GetClusterV();
                pClusterB1 = (*iIterB)->GetClusterU();
                pClusterB2 = (*iIterB)->GetClusterV();
                hitType1 = VIEW_U;
                hitType2 = VIEW_V;
            }

            if ((NULL == pCommonCluster) || (NULL == pClusterA1) || (NULL == pClusterA2) || (NULL == pClusterB1) || (NULL == pClusterB2))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            if ((pClusterA1 == pClusterB1) || (pClusterA2 == pClusterB2))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            try
            {
                // TODO Replace with something that returns "internal" vertices
                LArPointingCluster::Vertex vertexA1, vertexB1;
                LArPointingClusterHelper::GetClosestVertices(LArPointingCluster(pClusterA1), LArPointingCluster(pClusterB1), vertexA1, vertexB1);

                LArPointingCluster::Vertex vertexA2, vertexB2;
                LArPointingClusterHelper::GetClosestVertices(LArPointingCluster(pClusterA2), LArPointingCluster(pClusterB2), vertexA2, vertexB2);

                float transverseAB1(std::numeric_limits<float>::max()), transverseAB2(std::numeric_limits<float>::max());
                float longitudinalAB1(-std::numeric_limits<float>::max()), longitudinalAB2(-std::numeric_limits<float>::max());
                LArPointingClusterHelper::GetImpactParameters(vertexA1, vertexB1, longitudinalAB1, transverseAB1);
                LArPointingClusterHelper::GetImpactParameters(vertexA2, vertexB2, longitudinalAB2, transverseAB2);

                if (std::min(longitudinalAB1, longitudinalAB2) < m_minLongitudinalImpactParameter)
                    continue;

                float transverseBA1(std::numeric_limits<float>::max()), transverseBA2(std::numeric_limits<float>::max());
                float longitudinalBA1(-std::numeric_limits<float>::max()), longitudinalBA2(-std::numeric_limits<float>::max());
                LArPointingClusterHelper::GetImpactParameters(vertexB1, vertexA1, longitudinalBA1, transverseBA1);
                LArPointingClusterHelper::GetImpactParameters(vertexB2, vertexA2, longitudinalBA2, transverseBA2);

                if (std::min(longitudinalBA1, longitudinalBA2) < m_minLongitudinalImpactParameter)
                    continue;

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
                    continue;

                CartesianVector projectedPosition(0.f, 0.f, 0.f);
                float chiSquared(std::numeric_limits<float>::max());

                LArGeometryHelper::MergeTwoPositions(hitType1, hitType2, splitAtElementA ? vertexA1.GetPosition() : vertexB1.GetPosition(),
                     splitAtElementA ? vertexA2.GetPosition() : vertexB2.GetPosition(), projectedPosition, chiSquared);

                // TODO
            }
            catch (StatusCodeException &)
            {
                continue;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OvershootTracksTool::PassesElementCuts(TensorType::ElementList::const_iterator eIter, const ClusterList &usedClusters) const
{
    if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
        return false;

    if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
        return false;

    if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OvershootTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minLongitudinalImpactParameter = -1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLongitudinalImpactParameter", m_minLongitudinalImpactParameter));

    m_maxVertexXSeparation = 2.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexXSeparation", m_maxVertexXSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
