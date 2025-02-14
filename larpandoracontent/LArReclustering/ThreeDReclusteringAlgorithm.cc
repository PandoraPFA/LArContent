/**
 *  @file   larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.cc
 *
 *  @brief  Implementation file for the reclustering algorithm class.
 *
 *  $Log: $
 */

#include "Managers/ClusterManager.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

//Figure of merit type enum to string map
const std::unordered_map<std::string, ThreeDReclusteringAlgorithm::FigureOfMeritType> ThreeDReclusteringAlgorithm::m_stringToEnumMap = {
    {"cheated", ThreeDReclusteringAlgorithm::FigureOfMeritType::CHEATED}};

ThreeDReclusteringAlgorithm::ThreeDReclusteringAlgorithm() :
    m_pfoListName("ShowerParticles3D"),
    m_clusterListName("ShowerClusters3D"),
    m_fomThresholdForReclustering(0.3),
    m_minNumCaloHitsForReclustering(2),
    m_uClustersListName("ClustersU"),
    m_vClustersListName("ClustersV"),
    m_wClustersListName("ClustersW")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDReclusteringAlgorithm::~ThreeDReclusteringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDReclusteringAlgorithm::Run()
{
    //Only proceed if successfuly get shower pfo list, and it has at least one pfo
    const PfoList *pShowerPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_pfoListName, pShowerPfoList));
    if (!pShowerPfoList)
    {
        return STATUS_CODE_SUCCESS;
    }

    if (pShowerPfoList->empty())
    {
        return STATUS_CODE_SUCCESS;
    }

    m_pfosForReclusteringListName = "newShowerParticles3D";
    PfoList changedPfoList;

    //Save current pfo list name, so that this can be set as current again at the end of reclustering
    std::string initialPfoListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Pfo>(*this, initialPfoListName));

    //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
    const ClusterList *pShowerClusters(nullptr);
    PandoraContentApi::GetList(*this, m_clusterListName, pShowerClusters);

    if (!pShowerClusters)
    {
        return STATUS_CODE_NOT_FOUND;
    }

    for (const Pfo *const pShowerPfo : *pShowerPfoList)
    {
        ClusterList clusterList3D;
        LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);

        //Check if pfo passes cuts for reclustering. Also, some pfos are shower-like and yet are made of track-like 3D clusters. For the moment I don't want to deal with these.
        if ((!this->PassesCutsForReclustering(pShowerPfo)) ||
            (pShowerClusters->end() == std::find(pShowerClusters->begin(), pShowerClusters->end(), clusterList3D.front())))
            continue;

        CaloHitList initialCaloHitList;
        clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(initialCaloHitList);

        //Create a variable for the minimum figure of merit and initialize to initial FOM
        float minimumFigureOfMerit(this->GetFigureOfMerit(initialCaloHitList));

        //Free the hits in this cluster, so that they are not owned by the original pfo, and are available for reclustering
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, clusterList3D.front()));

        // Initialize reclustering with these local lists
        std::string currentClustersListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClustersListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListName));

        //Split the calo hit list into a new set of calo hit lists, taking the best outcome out of different algorithms
        std::vector<CaloHitList> newCaloHitListsVector, minimumFigureOfMeritCaloHitListsVector;
        minimumFigureOfMeritCaloHitListsVector.push_back(initialCaloHitList);

        for (auto *const pTool : m_algorithmToolVector)
        {
            try
            {
                pTool->Run(initialCaloHitList, newCaloHitListsVector);
            }
            catch (const StatusCodeException &)
            {
                std::cout << pTool->GetType() << ": Exception caught! Cannot run reclustering tool!: " << std::endl;
                continue;
            }

            //Calculate FOM for this vector of new CaloHitLists
            const float newFigureOfMerit(this->GetFigureOfMerit(newCaloHitListsVector));

            //Is this FOM smaller?
            if (newFigureOfMerit < minimumFigureOfMerit)
            {
                minimumFigureOfMerit = newFigureOfMerit;
                minimumFigureOfMeritCaloHitListsVector = newCaloHitListsVector;
            }
            newCaloHitListsVector.clear();
        }

        //If the new best calo hit lists outcome is equivalent to original, return the clusters to it. Else, create new vector of 3D clusters
        if ((minimumFigureOfMeritCaloHitListsVector.size() == 1) && (minimumFigureOfMeritCaloHitListsVector.at(0) == initialCaloHitList))
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pShowerPfo, clusterList3D.front()));
            continue;
        }

        const ClusterList reclusterClusterList(1, clusterList3D.front());
        std::string clusterListToSaveName, clusterListToDeleteName;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::InitializeFragmentation(*this, reclusterClusterList, clusterListToDeleteName, clusterListToSaveName));

        ClusterList newClustersList;
        for (CaloHitList &list : minimumFigureOfMeritCaloHitListsVector)
        {
            const Cluster *pCluster = nullptr;
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList = list;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
            newClustersList.push_back(pCluster);
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RebuildPfo(pShowerPfo, newClustersList));
        changedPfoList.push_back(pShowerPfo);
        newCaloHitListsVector.clear();
    }
    // Delete the original pfos from the reclustering
    for (const auto &pfo : changedPfoList)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Pfo>(*this, pfo, m_pfoListName));
    // Any remaining pfos (those that didnt undergo reclustering) get moved to the new pfo list
    if (!pShowerPfoList->empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_pfoListName, m_pfosForReclusteringListName));

    //Save list of Pfos after reclustering and save into list called m_pfoListName
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_pfosForReclusteringListName, m_pfoListName));

    //Set current list to be the same as before reclustering
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, initialPfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDReclusteringAlgorithm::BuildNewTwoDClusters(const Pfo *pPfoToRebuild, ClusterList &newClustersList)
{
    ClusterList clusterList2D;
    LArPfoHelper::GetTwoDClusterList(pPfoToRebuild, clusterList2D);

    for (const Cluster *const pTwoDCluster : clusterList2D)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfoToRebuild, pTwoDCluster));
        HitType hitType = LArClusterHelper::GetClusterHitType(pTwoDCluster);
        std::string clusterListName(hitType == TPC_VIEW_U ? m_uClustersListName
                : hitType == TPC_VIEW_V                   ? m_vClustersListName
                                                          : m_wClustersListName);

        std::string initialListName = "";
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, initialListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

        std::string originalListName, fragmentListName;
        ClusterList originalClusterList;
        originalClusterList.push_back(pTwoDCluster);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

        const OrderedCaloHitList &twoDClusterOrderedCaloHitList(pTwoDCluster->GetOrderedCaloHitList());
        OrderedCaloHitList leftoverCaloHitList = twoDClusterOrderedCaloHitList;

        int iCluster(0);
        for (const Cluster *pNewCluster3D : newClustersList)
        {
            if (!pNewCluster3D)
            {
                std::cout << "Error: found null pointer in list of new clusters!" << std::endl;
                continue;
            }
            PandoraContentApi::Cluster::Parameters parameters;
            CaloHitList newClusterCaloHitList3D;
            pNewCluster3D->GetOrderedCaloHitList().FillCaloHitList(newClusterCaloHitList3D);

            for (const CaloHit *const p3DCaloHit : newClusterCaloHitList3D)
            {
                for (const OrderedCaloHitList::value_type &mapEntry : twoDClusterOrderedCaloHitList)
                {
                    for (const CaloHit *const pCaloHit : *mapEntry.second)
                    {
                        if (pCaloHit == static_cast<const CaloHit *>(p3DCaloHit->GetParentAddress()))
                        {
                            parameters.m_caloHitList.push_back(static_cast<const CaloHit *>(pCaloHit));
                            leftoverCaloHitList.Remove(pCaloHit);
                            break;
                        }
                    }
                }
            }
            const Cluster *pNewTwoDCluster(nullptr);
            if (!parameters.m_caloHitList.empty())
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewTwoDCluster));
            }
            if (pNewTwoDCluster != nullptr && !parameters.m_caloHitList.empty() && hitType == TPC_VIEW_U)
            {
                m_newClustersUMap.insert(std::make_pair(iCluster, pNewTwoDCluster));
            }
            else if (pNewTwoDCluster != nullptr && !parameters.m_caloHitList.empty() && hitType == TPC_VIEW_V)
            {
                m_newClustersVMap.insert(std::make_pair(iCluster, pNewTwoDCluster));
            }
            else if (pNewTwoDCluster != nullptr && !parameters.m_caloHitList.empty() && hitType == TPC_VIEW_W)
            {
                m_newClustersWMap.insert(std::make_pair(iCluster, pNewTwoDCluster));
            }

            iCluster++;
        }

        //Check the leftover caloHits. Attach to the nearest cluster in the new cluster list.
        std::map<int, const Cluster *> clustersForLeftoverHitsMap(hitType == TPC_VIEW_U ? m_newClustersUMap
                : hitType == TPC_VIEW_V                                                 ? m_newClustersVMap
                                                                                        : m_newClustersWMap);

        if (!clustersForLeftoverHitsMap.empty())
        {
            for (const OrderedCaloHitList::value_type &mapEntry : leftoverCaloHitList)
            {
                for (const CaloHit *const pCaloHit : *mapEntry.second)
                {
                    const Cluster *pNearestCluster(nullptr);
                    double minimumDistance(std::numeric_limits<float>::max());
                    for (const auto &[clusterIndex, pNewTwoDCluster] : clustersForLeftoverHitsMap)
                    {
                        double dist = LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pNewTwoDCluster);
                        if (dist < minimumDistance)
                        {
                            minimumDistance = dist;
                            pNearestCluster = pNewTwoDCluster;
                        }
                    }
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pNearestCluster, pCaloHit));
                }
            }
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, initialListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDReclusteringAlgorithm::BuildNewPfos(const Pfo *pPfoToRebuild, ClusterList &newClustersList)
{
    const PfoList *pNewPfoList(nullptr);
    std::string newPfoListName = "changedShowerParticles3D";
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewPfoList, newPfoListName));

    const std::string originalClusterListName = "InitialCluster";
    int iCluster(0);

    for (const Cluster *pNewThreeDCluster : newClustersList)
    {
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        const bool isAvailableU((m_newClustersUMap.count(iCluster)) && m_newClustersUMap.at(iCluster)->IsAvailable());
        const bool isAvailableV((m_newClustersVMap.count(iCluster)) && m_newClustersVMap.at(iCluster)->IsAvailable());
        const bool isAvailableW((m_newClustersWMap.count(iCluster)) && m_newClustersWMap.at(iCluster)->IsAvailable());
        CaloHitList clusterUHits, clusterVHits, clusterWHits;
        if (isAvailableU)
        {
            m_newClustersUMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterUHits);
            pfoParameters.m_clusterList.push_back(m_newClustersUMap.at(iCluster));
        }
        if (isAvailableV)
        {
            m_newClustersVMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterVHits);
            pfoParameters.m_clusterList.push_back(m_newClustersVMap.at(iCluster));
        }
        if (isAvailableW)
        {
            m_newClustersWMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterWHits);
            pfoParameters.m_clusterList.push_back(m_newClustersWMap.at(iCluster));
        }
        pfoParameters.m_clusterList.push_back(pNewThreeDCluster);

        pfoParameters.m_particleId = pPfoToRebuild->GetParticleId();
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);

        const ParticleFlowObject *pNewPfo(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNewPfo));

        iCluster++;
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, newPfoListName, m_pfosForReclusteringListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDReclusteringAlgorithm::RebuildPfo(const Pfo *pPfoToRebuild, ClusterList &newClustersList)
{
    m_newClustersUMap.clear();
    m_newClustersVMap.clear();
    m_newClustersWMap.clear();

    //For each of the new 3D clusters, build a new 2D cluster in each view, and a new Pfo
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->BuildNewTwoDClusters(pPfoToRebuild, newClustersList));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->BuildNewPfos(pPfoToRebuild, newClustersList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDReclusteringAlgorithm::PassesCutsForReclustering(const pandora::ParticleFlowObject *const pShowerPfo)
{
    if (!LArPfoHelper::IsShower(pShowerPfo))
        return false;
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);

    if (clusterList3D.empty())
        return false;
    CaloHitList caloHitList3D;
    clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

    //Quality cuts
    if (caloHitList3D.size() < m_minNumCaloHitsForReclustering)
        return false;

    //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
    const ClusterList *pShowerClusters(nullptr);
    PandoraContentApi::GetList(*this, m_clusterListName, pShowerClusters);
    if (!pShowerClusters)
        return false;

    if (this->GetFigureOfMerit(caloHitList3D) < m_fomThresholdForReclustering)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDReclusteringAlgorithm::GetCheatedFigureOfMerit(const CaloHitList &mergedClusterCaloHitList3D)
{
    std::map<const pandora::MCParticle *, int> mainMcParticleMap;
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CaloHit *const pParentCaloHit2D = static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress());

        const MCParticle *pMainMCParticle{nullptr};
        try
        {
            pMainMCParticle = MCParticleHelper::GetMainMCParticle(pParentCaloHit2D);
        }
        catch (const StatusCodeException &err) // possible for some hits to have missing truth matching
        {
            if (err.GetStatusCode() == STATUS_CODE_NOT_INITIALIZED)
                continue;
            throw err;
        }

        std::map<const pandora::MCParticle *, int>::iterator it = mainMcParticleMap.find(pMainMCParticle);

        if (it != mainMcParticleMap.end())
            it->second++;
        else
            mainMcParticleMap.insert(std::make_pair(pMainMCParticle, 1));
    }
    if (mainMcParticleMap.size() == 0)
    {
        std::cout << "Cheated figure of merit is incalculable - all hits in a cluster have missing truth matching" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }
    const auto maxSharedHits =
        std::max_element(mainMcParticleMap.begin(), mainMcParticleMap.end(), [](const auto &x, const auto &y) { return x.second < y.second; });
    const int totalHits =
        std::accumulate(mainMcParticleMap.begin(), mainMcParticleMap.end(), 0, [](int tot, const auto &p) { return tot + p.second; });
    const float mainMcParticleFraction = (float)maxSharedHits->second / (float)totalHits;
    return (1.f - mainMcParticleFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(const std::string &figureOfMeritName, const CaloHitList &mergedClusterCaloHitList3D)
{
    float figureOfMerit(-999.f);

    FigureOfMeritType figureOfMeritType(FigureOfMeritType::CHEATED);
    auto stringToEnumIt = m_stringToEnumMap.find(figureOfMeritName);
    if (stringToEnumIt != m_stringToEnumMap.end())
    {
        figureOfMeritType = stringToEnumIt->second;
    }

    if (figureOfMeritType == FigureOfMeritType::CHEATED)
        figureOfMerit = this->GetCheatedFigureOfMerit(mergedClusterCaloHitList3D);
    else
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(const std::string &figureOfMeritName, const std::vector<CaloHitList> &newClustersCaloHitLists3D)
{
    std::vector<float> newClustersFigureOfMeritVector;
    for (auto clusterCaloHitLists3D : newClustersCaloHitLists3D)
    {
        newClustersFigureOfMeritVector.push_back(this->GetFigureOfMerit(figureOfMeritName, clusterCaloHitLists3D));
    }
    const float figureOfMerit(*std::min_element(newClustersFigureOfMeritVector.begin(), newClustersFigureOfMeritVector.end()));
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(const std::vector<CaloHitList> &newClustersCaloHitLists3D)
{
    std::vector<float> figureOfMeritVector;
    for (StringVector::const_iterator iter = m_figureOfMeritNames.begin(), iterEnd = m_figureOfMeritNames.end(); iter != iterEnd; ++iter)
    {
        figureOfMeritVector.push_back(this->GetFigureOfMerit(*iter, newClustersCaloHitLists3D));
    }

    const float figureOfMerit = *(std::min_element(figureOfMeritVector.begin(), figureOfMeritVector.end()));
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(const CaloHitList &mergedClusterCaloHitList3D)
{
    std::vector<float> figureOfMeritVector;
    for (StringVector::const_iterator iter = m_figureOfMeritNames.begin(), iterEnd = m_figureOfMeritNames.end(); iter != iterEnd; ++iter)
    {
        figureOfMeritVector.push_back(this->GetFigureOfMerit(*iter, mergedClusterCaloHitList3D));
    }
    const float figureOfMerit = *(std::min_element(figureOfMeritVector.begin(), figureOfMeritVector.end()));
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "FigureOfMeritNames", m_figureOfMeritNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FOMThresholdForReclustering", m_fomThresholdForReclustering));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinimumNumberCaloHitsForReclustering", m_minNumCaloHitsForReclustering));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UClustersListName", m_uClustersListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VClustersListName", m_vClustersListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WClustersListName", m_wClustersListName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "ClusteringTools", algorithmToolVector));

    for (auto algorithmTool : algorithmToolVector)
    {
        ClusteringTool *const pClusteringTool(dynamic_cast<ClusteringTool *>(algorithmTool));

        if (!pClusteringTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pClusteringTool);
    }
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
