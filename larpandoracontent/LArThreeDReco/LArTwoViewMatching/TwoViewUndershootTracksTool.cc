/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewUndershootTracksTool.cc
 *
 *  @brief  Implementation of the undershoot tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewUndershootTracksTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewUndershootTracksTool::TwoViewUndershootTracksTool() :
    //TwoViewThreeDKinkBaseTool(2),
    TwoViewThreeDKinkBaseTool(1),
    m_splitMode(false),
    m_maxTransverseImpactParameter(5.f),
    m_minImpactParameterCosTheta(0.5f),
    m_cosThetaCutForKinkSearch(0.75f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewUndershootTracksTool::GetIteratorListModifications(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList,
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
		//Adatta queste condizioniusando le variabili del calo matching (GetLocallyMatchedFraction(), GetMatchingScore(), GetNMatchedReUpsampledSamplingPoints(), GetTwoViewXOverlap().GetXOverlapFraction0() and 1)
                const unsigned int nMatchedReUpsampledSamplingPoints1((*iIter1)->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints());
                const unsigned int nMatchedReUpsampledSamplingPoints2((*iIter2)->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints());

                IteratorList::const_iterator iIterA((nMatchedReUpsampledSamplingPoints1 >= nMatchedReUpsampledSamplingPoints2) ? iIter1 : iIter2);
                IteratorList::const_iterator iIterB((nMatchedReUpsampledSamplingPoints1 >= nMatchedReUpsampledSamplingPoints2) ? iIter2 : iIter1);

                Particle particle(*(*iIterA), *(*iIterB));

                ClusterList commonClusters, clustersA, clustersB;
                commonClusters.push_back(particle.m_pCommonCluster);
                clustersA.push_back(particle.m_pClusterA);
                clustersB.push_back(particle.m_pClusterB);
	        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
	        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&commonClusters, "CommonClusters", RED);
	        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&clustersA, "ClustersA", BLUE);
	        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&clustersB, "ClustersB", GREEN);
	        PandoraMonitoringApi::ViewEvent(this->GetPandora());


                const LArPointingCluster pointingClusterA(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA));
                const LArPointingCluster pointingClusterB(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB));

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

                const bool isALowestInX(this->IsALowestInX(pointingClusterA, pointingClusterB));
                const CartesianVector splitPosition((vertexA.GetPosition() + vertexB.GetPosition()) * 0.5f);
                //m_majorityRulesMode = Whether to run in majority rules mode (always split overshoots, always merge undershoots)
                const bool isThreeDKink(m_majorityRulesMode ? false : this->IsThreeDKink(pAlgorithm, particle, splitPosition, isALowestInX));

                if (isThreeDKink != m_splitMode)
                    continue;

                // Construct the modification object
                Modification modification;

                if (m_splitMode)
                {
                    const TwoDSlidingFitResult &fitResult1(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster));
                    //const TwoDSlidingFitResult &fitResult2(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster2));

                    CartesianVector splitPosition1(0.f, 0.f, 0.f)/*, splitPosition2(0.f, 0.f, 0.f)*/;
                    if ((STATUS_CODE_SUCCESS != fitResult1.GetGlobalFitPositionAtX(splitPosition.GetX(), splitPosition1))/* ||
                        (STATUS_CODE_SUCCESS != fitResult2.GetGlobalFitPositionAtX(splitPosition.GetX(), splitPosition2))*/)
                    {
                        continue;
                    }

                    modification.m_splitPositionMap[particle.m_pCommonCluster].push_back(splitPosition1);
                    //modification.m_splitPositionMap[particle.m_pCommonCluster2].push_back(splitPosition2);
                }
                else
                {
                    const bool vertexAIsLowX(vertexA.GetPosition().GetX() < vertexB.GetPosition().GetX());
                    const Cluster *const pLowXCluster(vertexAIsLowX ? particle.m_pClusterA : particle.m_pClusterB);
                    const Cluster *const pHighXCluster(vertexAIsLowX ? particle.m_pClusterB : particle.m_pClusterA);
                    modification.m_clusterMergeMap[pLowXCluster].push_back(pHighXCluster);
                }

                modification.m_affectedClusters.push_back(particle.m_pClusterA);
                modification.m_affectedClusters.push_back(particle.m_pClusterB);
                modification.m_affectedClusters.push_back(particle.m_pCommonCluster);
                //modification.m_affectedClusters.push_back(particle.m_pCommonCluster2);

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

bool TwoViewUndershootTracksTool::IsThreeDKink(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const Particle &particle, const CartesianVector &splitPosition,
    const bool isALowestInX) const
{
    try
    {
        const TwoDSlidingFitResult &fitResultCommon(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster));
        //const TwoDSlidingFitResult &fitResultCommon2(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster2));
        const TwoDSlidingFitResult &lowXFitResult(isALowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA) : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB));
        const TwoDSlidingFitResult &highXFitResult(isALowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB) : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA));

	///QUESTO DEVE ESSERE ADATTATO!!! Non so ancora come
        const float minusX(this->GetXSamplingPoint(splitPosition, false, lowXFitResult, fitResultCommon));
        const float plusX(this->GetXSamplingPoint(splitPosition, true, highXFitResult, fitResultCommon));
        const float splitX(splitPosition.GetX());

        CartesianVector minus1(0.f, 0.f, 0.f), split1(0.f, 0.f, 0.f), plus1(0.f, 0.f, 0.f);
        //CartesianVector minus2(0.f, 0.f, 0.f), split2(0.f, 0.f, 0.f), plus2(0.f, 0.f, 0.f);
        CartesianVector minus2(0.f, 0.f, 0.f), split2(splitPosition), plus2(0.f, 0.f, 0.f);

        if ((STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(minusX, minus1)) ||
            (STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(splitX, split1)) ||
            (STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(plusX, plus1)) ||
//            (STATUS_CODE_SUCCESS != fitResultCommon2.GetGlobalFitPositionAtX(minusX, minus2)) ||
//            (STATUS_CODE_SUCCESS != fitResultCommon2.GetGlobalFitPositionAtX(splitX, split2)) ||
//            (STATUS_CODE_SUCCESS != fitResultCommon2.GetGlobalFitPositionAtX(plusX, plus2)) ||
            (STATUS_CODE_SUCCESS != lowXFitResult.GetGlobalFitPositionAtX(minusX, minus2)) ||
            (STATUS_CODE_SUCCESS != highXFitResult.GetGlobalFitPositionAtX(plusX, plus2)))
        {
            return false; // majority rules, by default
        }

        // Extract results
        const HitType hitType1(LArClusterHelper::GetClusterHitType(particle.m_pCommonCluster));
        //const HitType hitType2(LArClusterHelper::GetClusterHitType(particle.m_pCommonCluster2));
        const HitType hitType2(LArClusterHelper::GetClusterHitType(particle.m_pClusterA));

        CartesianVector minus(0.f, 0.f, 0.f), split(0.f, 0.f, 0.f), plus(0.f, 0.f, 0.f);
        float chi2Minus(std::numeric_limits<float>::max()), chi2Split(std::numeric_limits<float>::max()), chi2Plus(std::numeric_limits<float>::max());
        /*LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType1, hitType2, hitType3, minus1, minus2, minus3, minus, chi2Minus);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType1, hitType2, hitType3, split1, split2, split3, split, chi2Split);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType1, hitType2, hitType3, plus1, plus2, plus3, plus, chi2Plus);*/
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, minus1, minus2,  minus, chi2Minus);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, split1, split2,  split, chi2Split);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, plus1, plus2, plus, chi2Plus);

        // Apply final cuts
        const CartesianVector minusToSplit((split - minus).GetUnitVector());
        const CartesianVector splitToPlus((plus - split).GetUnitVector());
        const float dotProduct(minusToSplit.GetDotProduct(splitToPlus));

        if (dotProduct < m_cosThetaCutForKinkSearch)
            return true;
    }
    catch (StatusCodeException &s)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewUndershootTracksTool::Particle::Particle(const MatrixType::Element &elementA, const MatrixType::Element &elementB)
{
    m_pClusterA = (elementA.GetCluster1() != elementB.GetCluster1()) ? elementA.GetCluster1() : elementA.GetCluster2();
    m_pClusterB = (elementA.GetCluster1() != elementB.GetCluster1()) ? elementB.GetCluster1() : elementB.GetCluster2();
    m_pCommonCluster = (elementA.GetCluster1() == elementB.GetCluster1()) ? elementA.GetCluster1() : elementA.GetCluster2();
    //m_pCommonCluster2 = ((m_pClusterA != elementA.GetCluster1()) && (m_pCommonCluster != elementA.GetCluster1())) ? elementA.GetCluster1() :
        //((m_pClusterA != elementA.GetCluster2()) && (m_pCommonCluster != elementA.GetCluster2())) ? elementA.GetCluster2() : elementA.GetClusterW();

    //if ((m_pClusterA == m_pClusterB) || (m_pCommonCluster == m_pCommonCluster2))
    if (m_pClusterA == m_pClusterB)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewUndershootTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SplitMode", m_splitMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseImpactParameter", m_maxTransverseImpactParameter));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinImpactParameterCosTheta", m_minImpactParameterCosTheta));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CosThetaCutForKinkSearch", m_cosThetaCutForKinkSearch));

    return TwoViewThreeDKinkBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
