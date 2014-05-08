/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/ThreeDRemnantTracksAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional longitudinal tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArTrackMatching/ThreeDRemnantTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDRemnantTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDRemnantTracksAlgorithm::ExamineTensor()
{
    unsigned int repeatCounter(0);

    for (RemnantTensorToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; )
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

StatusCode ThreeDRemnantTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        RemnantTensorTool *pTensorManipulationTool(dynamic_cast<RemnantTensorTool*>(*iter));

        if (NULL == pTensorManipulationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pTensorManipulationTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    return ThreeDBaseAlgorithm<float>::ReadSettings(xmlHandle);
}

} // namespace lar
