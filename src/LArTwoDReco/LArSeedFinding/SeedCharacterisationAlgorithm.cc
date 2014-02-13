/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/SeedCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the seed characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArSeedFinding/SeedCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedCharacterisationAlgorithm::Run()
{
    // Get the seed association list
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    // Count branches for seed characterisation
    ClusterVector particleSeedVector(pSeedClusterList->begin(), pSeedClusterList->end());
    std::sort(particleSeedVector.begin(), particleSeedVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    SeedAssociationList seedAssociationList;
    this->GetSeedAssociationList(particleSeedVector, seedAssociationList);

    for (SeedAssociationList::const_iterator iter = seedAssociationList.begin(), iterEnd = seedAssociationList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pParent(iter->first);
        const ClusterVector &associatedClusters(iter->second);
        const unsigned int nAssociatedClusters(associatedClusters.size());

        // TODO Characterise seed using associated list of branches or prongs ...
    }

    // ATTN Use MC information for now...
    for (ClusterList::const_iterator iter = pSeedClusterList->begin(), iterEnd = pSeedClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
            const int pdgCode(pMCParticle->GetParticleId());

            // Maybe able to distinguish protons using reco properties? Unlikely to distinguish between mu+-/pi+- and e+-/photons?
            if ((std::abs(pdgCode) == PROTON) || (std::abs(pdgCode) == MU_MINUS) || (std::abs(pdgCode) == PI_PLUS))
            {
                pCluster->SetIsMipTrackFlag(true);
            }
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

void SeedCharacterisationAlgorithm::GetSeedAssociationList(const ClusterVector &particleSeedVector, SeedAssociationList &seedAssociationList) const
{
    const ClusterList *pNonSeedClusterList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this,
        m_nonSeedClusterListName, pNonSeedClusterList));

    ClusterUsageMap forwardUsageMap, backwardUsageMap;

    for (ClusterVector::const_iterator iter1 = particleSeedVector.begin(), iter1End = particleSeedVector.end(); iter1 != iter1End; ++iter1)
    {
        ClusterVector candidateClusters;

        if (NULL != pNonSeedClusterList)
            this->GetCandidateClusters(pNonSeedClusterList, candidateClusters);

        for (ClusterVector::const_iterator iter2 = particleSeedVector.begin(), iter2End = particleSeedVector.end(); iter2 != iter2End; ++iter2)
        {
            if (*iter1 != *iter2)
                candidateClusters.push_back(*iter2);
        }

        std::sort(candidateClusters.begin(), candidateClusters.end(), LArClusterHelper::SortByNOccupiedLayers);
        this->FindAssociatedClusters(*iter1, candidateClusters, forwardUsageMap, backwardUsageMap);
    }

    this->IdentifyClusterMerges(particleSeedVector, backwardUsageMap, seedAssociationList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SeedCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    

    return SeedBranchGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
