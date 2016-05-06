/**
 *  @file   LArContent/src/LArStitching/StitchingAlgorithm.cc
 * 
 *  @brief  Implementation of the Stitching algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArStitching/StitchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode StitchingAlgorithm::Run()
{
    StitchingInfo stitchingInfo;

    for (StitchingTool *const pStitchingTool : m_algorithmToolList)
    {
        pStitchingTool->Run(this, stitchingInfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "StitchingTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        StitchingTool *const pStitchingTool(dynamic_cast<StitchingTool*>(*iter));

        if (NULL == pStitchingTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pStitchingTool);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
