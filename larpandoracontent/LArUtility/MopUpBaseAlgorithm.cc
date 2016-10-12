/**
 *  @file   larpandoracontent/LArUtility/MopUpBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the mop up algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/MopUpBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
const std::string MopUpBaseAlgorithm::GetListName(const T *const pT) const
{
    std::string currentListName;
    const MANAGED_CONTAINER<const T*> *pCurrentList(nullptr);
    (void) PandoraContentApi::GetCurrentList(*this, pCurrentList, currentListName);

    if (pCurrentList && (pCurrentList->end() != std::find(pCurrentList->begin(), pCurrentList->end(), pT)))
        return currentListName;

    for (const std::string &listName : m_daughterListNames)
    {
        const MANAGED_CONTAINER<const T*> *pList(nullptr);
        (void) PandoraContentApi::GetList(*this, listName, pList);

        if (pList && (pList->end() != std::find(pList->begin(), pList->end(), pT)))
            return listName;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MopUpBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "DaughterListNames", m_daughterListNames));

    if (m_daughterListNames.empty())
    {
        std::cout << "MopUpBaseAlgorithm::ReadSettings - Must provide names of daughter object lists for use in mop up." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template const std::string MopUpBaseAlgorithm::GetListName(const CaloHit *const) const;
template const std::string MopUpBaseAlgorithm::GetListName(const MCParticle *const) const;
template const std::string MopUpBaseAlgorithm::GetListName(const Track *const) const;
template const std::string MopUpBaseAlgorithm::GetListName(const Cluster *const) const;
template const std::string MopUpBaseAlgorithm::GetListName(const ParticleFlowObject *const) const;
template const std::string MopUpBaseAlgorithm::GetListName(const Vertex *const) const;

} // namespace lar_content
