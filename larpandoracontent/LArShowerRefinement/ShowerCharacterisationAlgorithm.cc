/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the shower characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArShowerRefinement/ShowerCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerCharacterisationAlgorithm::ShowerCharacterisationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerCharacterisationAlgorithm::Run()
{
    /*
    PfoList electronPfos, gammaPfos;

    for (const std::string pfoListName : m_pfoListNames)
    {
        std::cout << "Pfo list name: " << pfoListName < std::endl;

        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ShowerCharacterisationAlgorithm: unable to find pfo list " << pfoListName << std::endl;

            continue;
        }

        std::cout << "Pfo list obtained" << std::endl;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            if (this->isElectron(pPfo))
                
        }


    }
    */
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    /*
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "OutputElectronListName", m_outputElectronListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "OutputGammaListName", m_outputGammaListName));
    */
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
