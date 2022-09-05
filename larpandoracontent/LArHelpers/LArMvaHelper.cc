/**
 *  @file   larpandoracontent/include/LArHelpers/LArMvaHelper.cc
 *
 *  @brief  Implementation for the lar mva helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

using namespace pandora;

namespace lar_content
{

StatusCode LArMvaHelper::ProcessAlgorithmToolListToMap(const Algorithm &algorithm, const TiXmlHandle &xmlHandle,
    const std::string &listName, StringVector &algorithmToolNameVector, AlgorithmToolMap &algorithmToolMap)
{
    // Fill a vector with names in the desired run order as well as the map

    if ("algorithm" != xmlHandle.ToNode()->ValueStr())
        return STATUS_CODE_NOT_ALLOWED;

    const TiXmlHandle algorithmListHandle = TiXmlHandle(xmlHandle.FirstChild(listName).Element());

    for (TiXmlElement *pXmlElement = algorithmListHandle.FirstChild("tool").Element(); nullptr != pXmlElement;
         pXmlElement = pXmlElement->NextSiblingElement("tool"))
    {
        AlgorithmTool *pAlgorithmTool(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateAlgorithmTool(algorithm, pXmlElement, pAlgorithmTool));
        std::string toolName = pXmlElement->Attribute("type");
        // If already exists, then make the second have the instance name attached
        if (algorithmToolMap.find(toolName) != algorithmToolMap.end())
            toolName += "_" + pAlgorithmTool->GetInstanceName();
        algorithmToolMap[toolName] = pAlgorithmTool;
        algorithmToolNameVector.push_back(toolName);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
