/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.cc
 *
 *  @brief  Implementation of the two view transverse tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArDiscreteProbabilityHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoViewTransverseTracksAlgorithm::TwoViewTransverseTracksAlgorithm() :
    m_nMaxMatrixToolRepeats(1000),
    m_downsampleFactor(5),
    m_minSamples(11),
    m_nPermutations(1000),
    m_localMatchingScoreThreshold(0.99f),
    m_maxDotProduct(0.998f),
    m_minOverallMatchingScore(0.1f),
    m_minOverallLocallyMatchedFraction(0.1f),
    m_randomNumberGenerator(static_cast<std::mt19937::result_type>(0))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const)
{
    m_randomNumberGenerator.seed(
        static_cast<std::mt19937::result_type>(pCluster1->GetOrderedCaloHitList().size() + pCluster2->GetOrderedCaloHitList().size()));

    TwoViewTransverseOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pCluster1, pCluster2, overlapResult));

    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapMatrix().SetOverlapResult(pCluster1, pCluster2, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TwoViewTransverseTracksAlgorithm::CalculateOverlapResult(
    const Cluster *const pCluster1, const Cluster *const pCluster2, TwoViewTransverseOverlapResult &overlapResult)
{
    UIntSet daughterVolumeIntersection;
    LArGeometryHelper::GetCommonDaughterVolumes(pCluster1, pCluster2, daughterVolumeIntersection);

    if (daughterVolumeIntersection.empty())
        return STATUS_CODE_NOT_FOUND;

    if (this->GetPrimaryAxisDotDriftAxis(pCluster1) > m_maxDotProduct || this->GetPrimaryAxisDotDriftAxis(pCluster2) > m_maxDotProduct)
        return STATUS_CODE_NOT_FOUND;

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);

    const TwoViewXOverlap twoViewXOverlap(xMin1, xMax1, xMin2, xMax2);
    if (twoViewXOverlap.GetXSpan0() < std::numeric_limits<float>::epsilon() || twoViewXOverlap.GetXSpan1() < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (twoViewXOverlap.GetTwoViewXOverlapSpan() < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;

    float xOverlapMin(twoViewXOverlap.GetTwoViewXOverlapMin());
    float xOverlapMax(twoViewXOverlap.GetTwoViewXOverlapMax());
    float zMin1(0.f), zMax1(0.f);
    float zMin2(0.f), zMax2(0.f);
    pCluster1->GetClusterSpanZ(xMin1, xMax1, zMin1, zMax1);
    pCluster2->GetClusterSpanZ(xMin2, xMax2, zMin2, zMax2);
    const CartesianVector boundingBoxMin1(xOverlapMin, 0.f, zMin1), boundingBoxMax1(xOverlapMax, 0.f, zMax1);
    const CartesianVector boundingBoxMin2(xOverlapMin, 0.f, zMin2), boundingBoxMax2(xOverlapMax, 0.f, zMax2);

    pandora::CaloHitList overlapHits1, overlapHits2;
    LArClusterHelper::GetCaloHitListInBoundingBox(pCluster1, boundingBoxMin1, boundingBoxMax1, overlapHits1);
    LArClusterHelper::GetCaloHitListInBoundingBox(pCluster2, boundingBoxMin2, boundingBoxMax2, overlapHits2);

    if (m_minSamples > std::min(overlapHits1.size(), overlapHits2.size()))
        return STATUS_CODE_NOT_FOUND;

    if (1 > m_downsampleFactor)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    const unsigned int nSamples(
        std::max(m_minSamples, static_cast<unsigned int>(std::min(overlapHits1.size(), overlapHits2.size())) / m_downsampleFactor));

    DiscreteProbabilityVector::AllFloatInputData inputData1;
    for (const pandora::CaloHit *const pCaloHit : overlapHits1)
        inputData1.emplace_back(pCaloHit->GetPositionVector().GetX(), pCaloHit->GetInputEnergy());

    DiscreteProbabilityVector::AllFloatInputData inputData2;
    for (const pandora::CaloHit *const pCaloHit : overlapHits2)
        inputData2.emplace_back(pCaloHit->GetPositionVector().GetX(), pCaloHit->GetInputEnergy());

    const DiscreteProbabilityVector discreteProbabilityVector1(inputData1, xOverlapMax, false);
    const DiscreteProbabilityVector discreteProbabilityVector2(inputData2, xOverlapMax, false);

    DiscreteProbabilityVector::ResamplingPoints resamplingPointsX;
    for (unsigned int iSample = 0; iSample < nSamples; ++iSample)
    {
        resamplingPointsX.emplace_back(
            (xOverlapMin + (xOverlapMax - xOverlapMin) * static_cast<float>(iSample + 1) / static_cast<float>(nSamples + 1)));
    }

    const DiscreteProbabilityVector resampledDiscreteProbabilityVector1(discreteProbabilityVector1, resamplingPointsX);
    const DiscreteProbabilityVector resampledDiscreteProbabilityVector2(discreteProbabilityVector2, resamplingPointsX);

    const float correlation(
        LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(resampledDiscreteProbabilityVector1, resampledDiscreteProbabilityVector2));

    const float pvalue(LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromPermutationTest(
        resampledDiscreteProbabilityVector1, resampledDiscreteProbabilityVector2, m_randomNumberGenerator, m_nPermutations));

    const float matchingScore(1.f - pvalue);
    if (matchingScore < m_minOverallMatchingScore)
        return STATUS_CODE_NOT_FOUND;

    const unsigned int nLocallyMatchedSamplingPoints(this->CalculateNumberOfLocallyMatchingSamplingPoints(
        resampledDiscreteProbabilityVector1, resampledDiscreteProbabilityVector2, m_randomNumberGenerator));
    const int nComparisons(static_cast<int>(resampledDiscreteProbabilityVector1.GetSize()) - (static_cast<int>(m_minSamples) - 1));
    if (1 > nComparisons)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float locallyMatchedFraction(static_cast<float>(nLocallyMatchedSamplingPoints) / static_cast<float>(nComparisons));
    if (locallyMatchedFraction < m_minOverallLocallyMatchedFraction)
        return STATUS_CODE_NOT_FOUND;

    overlapResult = TwoViewTransverseOverlapResult(
        matchingScore, m_downsampleFactor, nComparisons, nLocallyMatchedSamplingPoints, correlation, twoViewXOverlap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int TwoViewTransverseTracksAlgorithm::CalculateNumberOfLocallyMatchingSamplingPoints(const DiscreteProbabilityVector &discreteProbabilityVector1,
    const DiscreteProbabilityVector &discreteProbabilityVector2, std::mt19937 &randomNumberGenerator)
{
    if (discreteProbabilityVector1.GetSize() != discreteProbabilityVector2.GetSize() ||
        0 == discreteProbabilityVector1.GetSize() * discreteProbabilityVector2.GetSize())
        throw STATUS_CODE_INVALID_PARAMETER;

    if (m_minSamples > discreteProbabilityVector1.GetSize())
        throw STATUS_CODE_INVALID_PARAMETER;

    pandora::FloatVector localValues1, localValues2;
    unsigned int nMatchedComparisons(0);

    for (unsigned int iValue = 0; iValue < discreteProbabilityVector1.GetSize(); ++iValue)
    {
        localValues1.emplace_back(discreteProbabilityVector1.GetProbability(iValue));
        localValues2.emplace_back(discreteProbabilityVector2.GetProbability(iValue));
        if (localValues1.size() == m_minSamples)
        {
            float localPValue(0);
            try
            {
                localPValue = LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromPermutationTest(
                    localValues1, localValues2, randomNumberGenerator, m_nPermutations);
            }
            catch (const StatusCodeException &)
            {
                std::cout << "TwoViewTransverseTracksAlgorithm: failed to calculate correlation coefficient p-value for these numbers" << std::endl;
                ;
                std::cout << "----view 0: ";
                for (unsigned int iElement = 0; iElement < localValues1.size(); ++iElement)
                    std::cout << localValues1.at(iElement) << " ";
                std::cout << std::endl;
                std::cout << "----view 1: ";
                for (unsigned int iElement = 0; iElement < localValues2.size(); ++iElement)
                    std::cout << localValues2.at(iElement) << " ";
                std::cout << std::endl;
            }

            if ((1.f - localPValue) - m_localMatchingScoreThreshold > std::numeric_limits<float>::epsilon())
                nMatchedComparisons++;

            localValues1.erase(localValues1.begin());
            localValues2.erase(localValues2.begin());
        }
    }
    return nMatchedComparisons;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoViewTransverseTracksAlgorithm::GetPrimaryAxisDotDriftAxis(const pandora::Cluster *const pCluster)
{
    pandora::CartesianPointVector pointVector;
    LArClusterHelper::GetCoordinateVector(pCluster, pointVector);

    pandora::CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(pointVector, centroid, eigenValues, eigenVecs);

    const pandora::CartesianVector primaryAxis(eigenVecs.at(0));
    const pandora::CartesianVector driftAxis(1.f, 0.f, 0.f);
    return primaryAxis.GetDotProduct(driftAxis);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);
    for (MatrixToolVector::const_iterator iter = m_algorithmToolVector.begin(), iterEnd = m_algorithmToolVector.end(); iter != iterEnd;)
    {
        if ((*iter)->Run(this, this->GetMatchingControl().GetOverlapMatrix()))
        {
            iter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxMatrixToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewTransverseTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "TrackTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        TransverseMatrixTool *const pTransverseMatrixTool(dynamic_cast<TransverseMatrixTool *>(*iter));

        if (!pTransverseMatrixTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pTransverseMatrixTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NMaxMatrixToolRepeats", m_nMaxMatrixToolRepeats));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DownsampleFactor", m_downsampleFactor));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSamples", m_minSamples));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NPermutations", m_nPermutations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "LocalMatchingScoreThreshold", m_localMatchingScoreThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDotProduct", m_maxDotProduct));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinOverallMatchingScore", m_minOverallMatchingScore));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinOverallLocallyMatchedFraction", m_minOverallLocallyMatchedFraction));

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
