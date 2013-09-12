/**
 *  @file   LArContent/src/LArTwoDSeed/SeedCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the seed characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDSeed/SeedCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedCharacterisationAlgorithm::Run()
{
    // Get the seed association list
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    ClusterVector particleSeedVector(pSeedClusterList->begin(), pSeedClusterList->end());
    std::sort(particleSeedVector.begin(), particleSeedVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    SeedAssociationList seedAssociationList;
    this->GetSeedAssociationList(particleSeedVector, seedAssociationList);

    // Query the seed association list
    ClusterList trackSeeds, showerSeeds, relegatedSeeds;

    for (SeedAssociationList::const_iterator iter1 = seedAssociationList.begin(), iter1End = seedAssociationList.end(); iter1 != iter1End; ++iter1)
    {
        Cluster *pParticleSeed(iter1->first);
        const ClusterVector &associatedClusters(iter1->second);
        ClusterList associatedSeeds, associatedNonSeeds;

        for (ClusterVector::const_iterator iter2 = associatedClusters.begin(), iter2End = associatedClusters.end(); iter2 != iter2End; ++iter2)
        {
            Cluster *pAssociatedCluster(*iter2);

            if (pSeedClusterList->count(pAssociatedCluster))
            {
                associatedSeeds.insert(pAssociatedCluster);
            }
            else
            {
                associatedNonSeeds.insert(pAssociatedCluster);
            }
        }
std::cout << "SeedCharacterisation: NAssociatedSeeds " << associatedSeeds.size() << ", NAssociatedNonSeeds " << associatedNonSeeds.size() << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList tempList; tempList.insert(pParticleSeed);
PandoraMonitoringApi::VisualizeClusters(&tempList, "ParticleSeed", RED);
PandoraMonitoringApi::VisualizeClusters(&associatedSeeds, "AssociatedSeedClusters", GREEN);
PandoraMonitoringApi::VisualizeClusters(&associatedNonSeeds, "AssociatedNonSeedClusters", BLUE);
PandoraMonitoringApi::ViewEvent();
    }

    // Cluster list output
    if (!trackSeeds.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName, m_trackSeedClusterListName, trackSeeds));
    }

    if (!showerSeeds.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName, m_showerSeedClusterListName, showerSeeds));
    }

    if (!relegatedSeeds.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName, m_nonSeedClusterListName, relegatedSeeds));
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackSeedClusterListName", m_trackSeedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerSeedClusterListName", m_showerSeedClusterListName));

    return SeedBranchGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
