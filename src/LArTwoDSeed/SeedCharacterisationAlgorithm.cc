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
        Cluster *pParent(iter1->first);

        if (relegatedSeeds.count(pParent))
            continue;

        const ClusterVector &associatedClusters(iter1->second);
        ClusterList associatedSeeds, associatedNonSeeds;
        ClusterList seedProngs, seedBranches, nonSeedProngs, nonSeedBranches;

        for (ClusterVector::const_iterator iter2 = associatedClusters.begin(), iter2End = associatedClusters.end(); iter2 != iter2End; ++iter2)
        {
            Cluster *pDaughter(*iter2);

            bool isProng(false), isBranch(false);
            const CartesianVector pInnerCentroid(pParent->GetCentroid(pParent->GetInnerPseudoLayer()));
            const CartesianVector pOuterCentroid(pParent->GetCentroid(pParent->GetOuterPseudoLayer()));
            const CartesianVector dInnerCentroid(pDaughter->GetCentroid(pDaughter->GetInnerPseudoLayer()));
            const CartesianVector dOuterCentroid(pDaughter->GetCentroid(pDaughter->GetOuterPseudoLayer()));

            // Prongs
            const float rsqInnerInner((pInnerCentroid - dInnerCentroid).GetMagnitudeSquared());
            const float rsqInnerOuter((pInnerCentroid - dOuterCentroid).GetMagnitudeSquared());
            const float rsqOuterInner((pOuterCentroid - dInnerCentroid).GetMagnitudeSquared());
            const float rsqOuterOuter((pOuterCentroid - dOuterCentroid).GetMagnitudeSquared()); 

            if (rsqInnerInner < 2.5 * 2.5 || rsqInnerOuter < 2.5 * 2.5 ||
                rsqOuterInner < 2.5 * 2.5 || rsqOuterOuter < 2.5 * 2.5) 
            {
                isProng = true;
            }

            // Branches
            const float pOuter(LArClusterHelper::GetClosestDistance(pOuterCentroid, pDaughter));
            const float pInner(LArClusterHelper::GetClosestDistance(pInnerCentroid, pDaughter));
            const float dOuter(LArClusterHelper::GetClosestDistance(dOuterCentroid, pParent));
            const float dInner(LArClusterHelper::GetClosestDistance(dInnerCentroid, pParent));

            if ((pInner > 2.5) && (pOuter > 2.5) && ((dInner < 2.5) || (dOuter < 2.5)))
            {
                isBranch = true;
            }

            // Fill containers
            if (pSeedClusterList->count(pDaughter))
            {
                associatedSeeds.insert(pDaughter);
                if (isProng) seedProngs.insert(pDaughter);
                if (isBranch) seedBranches.insert(pDaughter);
            }
            else
            {
                associatedNonSeeds.insert(pDaughter);
                if (isProng) nonSeedProngs.insert(pDaughter);
                if (isBranch) nonSeedBranches.insert(pDaughter);
            }
        }

        // Classify the particle seed
std::cout << "SeedCharacterisation: NAssociatedSeeds " << associatedSeeds.size() << " (br: " << seedBranches.size() << ", pr: " << seedProngs.size()
          << "), NAssociatedNonSeeds " << associatedNonSeeds.size() << " (br: " << nonSeedBranches.size() << ", pr: " << nonSeedProngs.size() << ")" << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList tempList; tempList.insert(pParent);
PandoraMonitoringApi::VisualizeClusters(&tempList, "ParticleSeed", RED);
PandoraMonitoringApi::VisualizeClusters(&associatedSeeds, "AssociatedSeedClusters", GREEN);
PandoraMonitoringApi::VisualizeClusters(&associatedNonSeeds, "AssociatedNonSeedClusters", BLUE);
PandoraMonitoringApi::VisualizeClusters(&seedProngs, "seedProngs", CYAN);
PandoraMonitoringApi::VisualizeClusters(&seedBranches, "seedBranches", MAGENTA);
PandoraMonitoringApi::VisualizeClusters(&nonSeedProngs, "nonSeedProngs", CYAN);
PandoraMonitoringApi::VisualizeClusters(&nonSeedBranches, "nonSeedBranches", MAGENTA);
PandoraMonitoringApi::ViewEvent();
    }

//std::cout << "SeedCharacterisation: results " << std::endl;
//std::cout << " p.size() " << particleSeedVector.size() << " (t.size() + s.size() + r.size()) " << (trackSeeds.size() + showerSeeds.size() + relegatedSeeds.size()) << std::endl;
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//PandoraMonitoringApi::VisualizeClusters(&trackSeeds, "trackSeeds", RED);
//PandoraMonitoringApi::VisualizeClusters(&showerSeeds, "showerSeeds", BLUE);
//PandoraMonitoringApi::VisualizeClusters(&relegatedSeeds, "relegatedSeeds", GREEN);
//PandoraMonitoringApi::ViewEvent();

//    // Sanity check
//    if (particleSeedVector.size() != (trackSeeds.size() + showerSeeds.size() + relegatedSeeds.size()))
//        return STATUS_CODE_FAILURE;

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
