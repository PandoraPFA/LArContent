/**
 *  @file   LArContent/src/LArCheating/CheatingClusterCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCheating/CheatingClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode CheatingClusterCreationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    MCParticleToHitListMap mcParticleToHitListMap;

    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            if (!PandoraContentApi::IsAvailable(*this, *iter))
                continue;

            this->SimpleMCParticleCollection(*iter, mcParticleToHitListMap);
        }
        catch (StatusCodeException &)
        {
        }
    }

    this->CreateClusters(mcParticleToHitListMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingClusterCreationAlgorithm::SimpleMCParticleCollection(CaloHit *const pCaloHit, MCParticleToHitListMap &mcParticleToHitListMap) const
{
    const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

    if (!this->SelectMCParticlesForClustering(pMCParticle))
        return;

    mcParticleToHitListMap[pMCParticle].insert(pCaloHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingClusterCreationAlgorithm::SelectMCParticlesForClustering(const MCParticle *const pMCParticle) const
{
    if (m_particleIdList.empty())
        return true;

    for (IntVector::const_iterator iter = m_particleIdList.begin(), iterEnd = m_particleIdList.end(); iter != iterEnd; ++iter)
    {
        if (pMCParticle->GetParticleId() == *iter)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingClusterCreationAlgorithm::CreateClusters(const MCParticleToHitListMap &mcParticleToHitListMap) const
{
    for (MCParticleToHitListMap::const_iterator iter = mcParticleToHitListMap.begin(), iterEnd = mcParticleToHitListMap.end(); 
         iter != iterEnd; ++iter)
    {
        const MCParticle *pMCParticle = iter->first;
        const CaloHitList &caloHitList = iter->second;

        if (caloHitList.empty())
            continue;

        Cluster *pCluster = NULL;
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList = caloHitList;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

        switch (pMCParticle->GetParticleId())
        {
        case PHOTON:
            pCluster->SetIsFixedPhotonFlag(true);
            break;

        case E_PLUS:
        case E_MINUS:
            pCluster->SetIsFixedElectronFlag(true);
            break;

        case MU_PLUS:
        case MU_MINUS:
            pCluster->SetIsFixedMuonFlag(true);
            break;

        default:
            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_particleIdList.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ParticleIdList", m_particleIdList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
