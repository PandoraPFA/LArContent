/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;

    if (m_inputClusterListName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pClusterList));
    }

    /*
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster* pCluster = *iter;

        // TODO: 'SetIsMipTrackFlag' for tracks and 'SetIsFixedElectronFlag' for showers
        
    }
    */

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_inputClusterListName.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputClusterListName", m_inputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
