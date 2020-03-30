/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.cc
 *
 *  @brief  Implementation of the two view transverse tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoViewTransverseTracksAlgorithm::TwoViewTransverseTracksAlgorithm() :
    m_nMaxMatrixToolRepeats(1000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const)
{
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xSpan1(xMax1 - xMin1), xSpan2(xMax2 - xMin2);

    if ((xSpan1 < std::numeric_limits<float>::epsilon()) || (xSpan2 < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // TODO Create a fully-fledged overlap result object and fill it with useful information
    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));

    if (xOverlap > std::numeric_limits<float>::epsilon())
        this->GetMatchingControl().GetOverlapMatrix().SetOverlapResult(pCluster1, pCluster2, xOverlap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (MatrixToolVector::const_iterator iter = m_algorithmToolVector.begin(), iterEnd = m_algorithmToolVector.end(); iter != iterEnd; )
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
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewTransverseTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        TransverseMatrixTool *const pTransverseMatrixTool(dynamic_cast<TransverseMatrixTool*>(*iter));

        if (!pTransverseMatrixTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pTransverseMatrixTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxMatrixToolRepeats", m_nMaxMatrixToolRepeats));

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
