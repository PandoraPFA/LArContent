/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/CheatingClusterCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating cluster characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCheating/CheatingClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode CheatingClusterCharacterisationAlgorithm::Run()
{
    for (StringVector::const_iterator listIter = m_inputClusterListNames.begin(), listIterEnd = m_inputClusterListNames.end(); listIter != listIterEnd; ++listIter)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *listIter, pClusterList));

        for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        {
            Cluster *const pCluster(*iter);

            if (this->IsClearTrack(pCluster))
            {
                pCluster->SetIsFixedMuonFlag(true);
            }
            else
            {
                pCluster->SetIsFixedMuonFlag(false);
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    try
    {
        // ATTN Slightly curious definition of a clear track, but this is most-likely what is needed for shower-growing
        const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

        if ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())))
            return true;
    }
    catch (StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
