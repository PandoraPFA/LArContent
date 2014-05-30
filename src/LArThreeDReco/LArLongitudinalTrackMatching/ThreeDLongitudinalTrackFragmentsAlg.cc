/**
 *  @file   LArContent/src/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTrackFragmentsAlg.cc
 *
 *  @brief  Implementation of the three dimensional longitudinal track fragments algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTrackFragmentsAlg.h"

using namespace pandora;

namespace lar
{

void ThreeDLongitudinalTrackFragmentsAlg::CalculateOverlapResult(Cluster */*pClusterU*/, Cluster */*pClusterV*/, Cluster */*pClusterW*/)
{
    try
    {
//        LongitudinalOverlapResult overlapResult;
//        this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
//
//        if (overlapResult.IsInitialized())
//            m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode() || 
              STATUS_CODE_NOT_INITIALIZED == statusCodeException.GetStatusCode()))
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLongitudinalTrackFragmentsAlg::ExamineTensor()
{
    unsigned int repeatCounter(0);

    for (TensorToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, m_overlapTensor))
        {
            iter = m_algorithmToolList.begin();

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDLongitudinalTrackFragmentsAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        LongitudinalFragmentTensorTool *pTensorManipulationTool(dynamic_cast<LongitudinalFragmentTensorTool*>(*iter));

        if (NULL == pTensorManipulationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pTensorManipulationTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    return ThreeDTracksBaseAlgorithm<float>::ReadSettings(xmlHandle);
}

} // namespace lar
