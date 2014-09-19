/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/KinkSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the kink splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/KinkSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

KinkSplittingAlgorithm::KinkSplittingAlgorithm() :
    m_maxScatterRms(0.2f),
    m_maxScatterCosTheta(0.905f),
    m_maxSlidingCosTheta(0.985f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkSplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, CartesianVector& splitPosition) const
{
    // Search for scatters by scanning over the layers in the sliding fit result
    const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());
    const int minLayer(layerFitResultMap.begin()->first), maxLayer(layerFitResultMap.rbegin()->first);

    const int nLayersHalfWindow(slidingFitResult.GetLayerFitHalfWindow());
    const int nLayersSpanned(1 + maxLayer - minLayer);

    if (nLayersSpanned <= 2 * nLayersHalfWindow)
        return STATUS_CODE_NOT_FOUND;

    bool foundSplit(false);

    float bestCosTheta(1.f);

    for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        const int iLayer(iter->first);

        const float rL(slidingFitResult.GetL(iLayer));
        const float rL1(slidingFitResult.GetL(iLayer - nLayersHalfWindow));
        const float rL2(slidingFitResult.GetL(iLayer + nLayersHalfWindow));

        try
        {
            CartesianVector centralPosition(0.f,0.f,0.f);
            CartesianVector firstDirection(0.f,0.f,0.f);
            CartesianVector secondDirection(0.f,0.f,0.f);

            slidingFitResult.GetGlobalFitPosition(rL, centralPosition);
            slidingFitResult.GetGlobalFitDirection(rL1, firstDirection);
            slidingFitResult.GetGlobalFitDirection(rL2, secondDirection);

            const float cosTheta(firstDirection.GetDotProduct(secondDirection));
            const float rms1(slidingFitResult.GetFitRms(rL1));
            const float rms2(slidingFitResult.GetFitRms(rL2));
            const float rms(std::max(rms1, rms2));

            float rmsCut(m_maxScatterRms);

            if (cosTheta > m_maxScatterCosTheta)
            {
                rmsCut *= ((m_maxSlidingCosTheta > cosTheta) ? (m_maxSlidingCosTheta - cosTheta) /
                        (m_maxSlidingCosTheta - m_maxScatterCosTheta) : 0.f);
            }

            if (rms < rmsCut && cosTheta < bestCosTheta)
            {
                bestCosTheta = cosTheta;
                splitPosition = centralPosition;
                foundSplit = true;
            }
        }
        catch (StatusCodeException &)
        {
        }
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxScatterRms", m_maxScatterRms));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxScatterCosTheta", m_maxScatterCosTheta));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSlidingCosTeta", m_maxSlidingCosTheta));

    return TwoDSlidingFitSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
