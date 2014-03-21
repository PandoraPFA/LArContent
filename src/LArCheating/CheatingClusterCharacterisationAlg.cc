/**
 *  @file   LArContent/src/LArCheating/CheatingClusterCharacterisationAlg.cc
 * 
 *  @brief  Implementation of the cheater for the cluster characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArCheating/CheatingClusterCharacterisationAlg.h"

using namespace pandora;

namespace lar
{

StatusCode CheatingClusterCharacterisationAlg::Run()
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

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
            const int pdgCode(pMCParticle->GetParticleId());

            // Muons, pions and protons are track-like
            if ((std::abs(pdgCode) == PROTON) || (std::abs(pdgCode) == MU_MINUS) || (std::abs(pdgCode) == PI_PLUS))
            {
                pCluster->SetIsMipTrackFlag(true);
            }

            // Electrons and photons are shower-like
            else if ((std::abs(pdgCode) == E_MINUS) || (std::abs(pdgCode) == PHOTON))
            {
                pCluster->SetIsFixedElectronFlag(true);
            }
        }
        catch (StatusCodeException &)
        {
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterCharacterisationAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_inputClusterListName.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputClusterListName", m_inputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
