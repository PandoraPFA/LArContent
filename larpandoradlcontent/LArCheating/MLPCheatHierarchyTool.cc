/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPCheatHierarchyTool.cc
 *
 *  @brief  Implementation of the cheat hierarchy tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoradlcontent/LArCheating/MLPCheatHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

MLPCheatHierarchyTool::MLPCheatHierarchyTool() :
    m_mcParticleListName("Input"),
    m_neutrinoPfoListName("NeutrinoParticles3D"),
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D"})
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPCheatHierarchyTool::Run(const PfoToMCParticleMap &pfoToMCParticleMap, const PfoToPfoMap &childToParentPfoMap,
    const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo, bool &isTrueLink, bool &trueParentOrientation, bool &trueChildOrientation)
{
    // If we dont know the Pfo->MCParticle match
    if ((pfoToMCParticleMap.find(parentPfo.GetPfo()) == pfoToMCParticleMap.end()) ||
        (pfoToMCParticleMap.find(childPfo.GetPfo()) == pfoToMCParticleMap.end()))
    {
        return STATUS_CODE_NOT_FOUND;
    }

    // If we dont know the true parentPfo
    if (childToParentPfoMap.find(childPfo.GetPfo()) == childToParentPfoMap.end())
        return STATUS_CODE_NOT_FOUND;

    isTrueLink = childToParentPfoMap.at(childPfo.GetPfo()) == parentPfo.GetPfo();

    // What is the true orientation of the parent?
    trueParentOrientation = this->IsUpstreamTrueVertex(pfoToMCParticleMap, parentPfo.GetPfo(),
        parentPfo.GetUpstreamVertex(), parentPfo.GetDownstreamVertex());
    
    // What is the true orientation of the child? 
    trueChildOrientation = this->IsUpstreamTrueVertex(pfoToMCParticleMap, childPfo.GetPfo(),
        childPfo.GetUpstreamVertex(), childPfo.GetDownstreamVertex());
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPCheatHierarchyTool::IsUpstreamTrueVertex(const PfoToMCParticleMap &pfoToMCParticleMap, const ParticleFlowObject *const pPfo,
    const CartesianVector &upstreamVertex, const CartesianVector &downstreamVertex)
{
    const MCParticle *const pMatchedMCParticle(pfoToMCParticleMap.at(pPfo));

    const float upstreamSepSq((pMatchedMCParticle->GetVertex() - upstreamVertex).GetMagnitudeSquared());
    const float downstreamSepSq((pMatchedMCParticle->GetVertex() - downstreamVertex).GetMagnitudeSquared());    

    return (upstreamSepSq < downstreamSepSq);
}
   
//------------------------------------------------------------------------------------------------------------------------------------------

void MLPCheatHierarchyTool::FillHierarchyMap(const Algorithm *const pAlgorithm) const
{
    // Get MCParticles
    const MCParticleList *pMCParticleList(nullptr);
    if (!this->GetMCParticleList(pAlgorithm, pMCParticleList))
        return;

    // Get reconstructed neutrino
    const ParticleFlowObject *pNeutrinoPfo(nullptr);
    if (!this->GetNeutrinoPfo(pAlgorithm, pNeutrinoPfo))
        return;
    
    // Match our PFPs to our MCParticles 
    PfoToMCParticleMap pfoToMCParticleMap;
    this->MatchPFParticles(pAlgorithm, pfoToMCParticleMap);

    // Do MCParticle visible hierarchy building
    MCToMCMap childToParentMCMap;
    this->GetVisibleMCHierarchy(pfoToMCParticleMap, childToParentMCMap);

    // Do Pfo visible hierarchy building
    PfoToPfoMap childToParentPfoMap;
    this->GetVisiblePfoHierarchy(pNeutrinoPfo, pfoToMCParticleMap, childToParentMCMap, childToParentPfoMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPCheatHierarchyTool::GetMCParticleList(const Algorithm *const pAlgorithm, const MCParticleList *& pMCParticleList) const
{
    if (PandoraContentApi::GetList(*pAlgorithm, m_mcParticleListName, pMCParticleList) != STATUS_CODE_SUCCESS)
        return false;

    if (!pMCParticleList || pMCParticleList->empty())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool MLPCheatHierarchyTool::GetNeutrinoPfo(const Algorithm *const pAlgorithm, const ParticleFlowObject *&pNeutrinoPfo) const
{
    const PfoList *pPfoList = nullptr;

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, 
        PandoraContentApi::GetList(*pAlgorithm, m_neutrinoPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
        return false;

    // ATTN Enforces that only one pfo, of neutrino-type, be in the specified input list
    pNeutrinoPfo = (1 == pPfoList->size()) ? *(pPfoList->begin()) : nullptr;

    if (!pNeutrinoPfo || !LArPfoHelper::IsNeutrino(pNeutrinoPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return true;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void MLPCheatHierarchyTool::MatchPFParticles(const Algorithm *const pAlgorithm, PfoToMCParticleMap &pfoToMCParticleMap) const
{
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*pAlgorithm, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        if (!pPfoList)
            continue;

        for (const ParticleFlowObject * pPfo : *pPfoList)
        {
            // Demand that we have 3D hits
            if (this->GetNSpacepoints(pPfo) == 0)
                continue;

            // Fill contribution map
            LArMCParticleHelper::MCContributionMap contributionMap;

            for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            {
                CaloHitList caloHitList;
                LArPfoHelper::GetCaloHits(pPfo, hitType, caloHitList);

                for (const CaloHit *const pCaloHit : caloHitList)
                {
                    try
                    {
                        const MCParticle *pHitMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

                        // Sometimes we have to fold back EM shower hierarchies
                        if (this->IsEMParticle(pHitMCParticle))
                            pHitMCParticle = this->GetLeadEMParticle(pHitMCParticle);

                        contributionMap[pHitMCParticle].emplace_back(pCaloHit);
                    }
                    catch (...) { continue; }
                }
            }

            // Work out who the match is
            unsigned int highestHit(0);
            float highestEnergy(-1.f);
            const MCParticle *pMatchedMCParticle(nullptr);

            for (const auto &entry : contributionMap)
            {
                if (entry.second.size() > highestHit)
                {
                    highestHit = entry.second.size();
                    highestEnergy = this->SumEnergy(entry.second);
                    pMatchedMCParticle = entry.first;
                }
                else if (entry.second.size() == highestHit)
                {
                    const float totalEnergy(this->SumEnergy(entry.second));

                    if (totalEnergy > highestEnergy)
                    {
                        highestHit = entry.second.size();
                        highestEnergy = totalEnergy;
                        pMatchedMCParticle = entry.first;
                    }
                }
            }

            if (!pMatchedMCParticle)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            pfoToMCParticleMap[pPfo] = pMatchedMCParticle;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPCheatHierarchyTool::GetNSpacepoints(const ParticleFlowObject *const pPfo) const
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    int total3DHits(0);

    for (const Cluster *const pCluster3D : clusterList3D)
        total3DHits += pCluster3D->GetNCaloHits();

    return total3DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPCheatHierarchyTool::IsEMParticle(const MCParticle *const pMCParticle) const
{
    return ((pMCParticle->GetParticleId() == 22) || (std::abs(pMCParticle->GetParticleId()) == 11));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* MLPCheatHierarchyTool::GetLeadEMParticle(const MCParticle *const pMCParticle) const
{
    const MCParticle *pThisMCParticle(pMCParticle);
    const MCParticle *pParentMCParticle(pMCParticle);

    while (this->IsEMParticle(pParentMCParticle)) 
    {
        const MCParticleList &parentList(pThisMCParticle->GetParentList());

        if (parentList.size() != 1)
        {
            std::cout << "EM shower has no parent?" << std::endl;
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }

        pParentMCParticle = *parentList.begin();

        if (this->IsEMParticle(pParentMCParticle))
            pThisMCParticle = pParentMCParticle;
    }

    return pThisMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPCheatHierarchyTool::SumEnergy(const CaloHitList &caloHitList) const
{
    float totalEnergy(0.f);

    for (const CaloHit *const pCaloHit : caloHitList)
        totalEnergy += pCaloHit->GetElectromagneticEnergy();

    return totalEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPCheatHierarchyTool::GetVisibleMCHierarchy(const PfoToMCParticleMap &pfoToMCParticleMap, MCToMCMap &childToParentMCMap) const
{
    MCParticleList reconstructedMCList;

    for (const auto &entry : pfoToMCParticleMap)
    {
        if (std::find(reconstructedMCList.begin(), reconstructedMCList.end(), entry.second) == reconstructedMCList.end())
            reconstructedMCList.emplace_back(entry.second);
    }

    for (const MCParticle *const pMCParticle : reconstructedMCList)
    {
        bool foundParent(false);
        const MCParticle *pThisMCParticle(pMCParticle);
        const MCParticle *pParentMCParticle(nullptr);

        while (!foundParent)
        {            
            const MCParticleList &parentList(pThisMCParticle->GetParentList());

            if (parentList.size() != 1)
            {
                std::cout << "No MCParent for MCParticle" << std::endl;
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);
            }

            pParentMCParticle = *parentList.begin();
            
            foundParent = (LArMCParticleHelper::IsNeutrino(pParentMCParticle) || 
                (std::find(reconstructedMCList.begin(), reconstructedMCList.end(), pParentMCParticle) != reconstructedMCList.end()));

            pThisMCParticle = pParentMCParticle;
        }

        childToParentMCMap[pMCParticle] = pParentMCParticle;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPCheatHierarchyTool::GetVisiblePfoHierarchy(const ParticleFlowObject *const pNeutrinoPfo, const PfoToMCParticleMap &pfoToMCParticleMap,
    const MCToMCMap &childToParentMCMap, PfoToPfoMap &childToParentPfoMap) const
{
    std::map<const MCParticle *, PfoList> mcParticleToPfoListMap;

    for (const auto &entry : pfoToMCParticleMap)
        mcParticleToPfoListMap[entry.second].emplace_back(entry.first);

    // First handle any split paricles
    for (const auto &entry : mcParticleToPfoListMap)
    {
        if (entry.second.size() == 1)
            continue;

        this->BuildSplitHierarchy(entry, childToParentPfoMap);
    }

    // Now assign remaining
    for (const auto &entry : mcParticleToPfoListMap)
    {
        const MCParticle *const pMatchedMCParticle(entry.first);

        if (childToParentMCMap.find(pMatchedMCParticle) == childToParentMCMap.end())
        {
            std::cout << "have found a mcparticle with no parent" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        const MCParticle *const pMatchedMCParent(childToParentMCMap.at(pMatchedMCParticle));

        for (const ParticleFlowObject *const pChildPfo : entry.second)
        {
            // If already assigned, move on
            if (childToParentPfoMap.find(pChildPfo) != childToParentPfoMap.end())
                continue;

            if (LArMCParticleHelper::IsNeutrino(pMatchedMCParent))
            {
                // If primary, great! Assign neutrino pfo
                childToParentPfoMap[pChildPfo] = pNeutrinoPfo;
            }
            else
            {
                // Else, let's find the parent pfo
                if (mcParticleToPfoListMap.find(pMatchedMCParent) == mcParticleToPfoListMap.end())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const PfoList &parentPfoList(mcParticleToPfoListMap.at(pMatchedMCParent));
           
                if (parentPfoList.size() == 1)
                    childToParentPfoMap[pChildPfo] = *parentPfoList.begin();
                else
                    childToParentPfoMap[pChildPfo] = this->BestParentInSplitHierarchy(pChildPfo, parentPfoList);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPCheatHierarchyTool::BuildSplitHierarchy(const std::pair<const MCParticle*, PfoList> &splitPfo, PfoToPfoMap &childToParentPfoMap) const
{
    // Order segments wrt distance to the true vertex
    std::vector<std::pair<const ParticleFlowObject*, float>> pfoComponents;

    for (const ParticleFlowObject *const pSplitPfo : splitPfo.second)
    {
        ClusterList clusterList;
        LArPfoHelper::GetThreeDClusterList(pSplitPfo, clusterList);

        if (clusterList.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);
              
        pfoComponents.emplace_back(std::make_pair(pSplitPfo, LArClusterHelper::GetClosestDistance(splitPfo.first->GetVertex(), clusterList)));
    }

    std::sort(pfoComponents.begin(), pfoComponents.end(), 
        [](const std::pair<const ParticleFlowObject*, float> &a, const std::pair<const ParticleFlowObject*, float> &b)
        { 
            return a.second < b.second;
        });

    // Now build pfo hierarchy, assuming closest is the leading pfo
    childToParentPfoMap[pfoComponents.at(1).first] = pfoComponents.at(0).first;

    for (unsigned int i = 2; i < pfoComponents.size(); ++i)
    {
        const ParticleFlowObject *const pChildPfo(pfoComponents.at(i).first);
        double closestSep = std::numeric_limits<double>::max();
        const ParticleFlowObject *pParentPfo(nullptr);

        for (unsigned int j = 0; j < i; ++j)
        {
            const ParticleFlowObject *const pPotentialParentPfo(pfoComponents.at(j).first);
            const float thisSep(this->GetClosestDistance(pChildPfo, pPotentialParentPfo));

            if (thisSep < closestSep)
            {
                closestSep = thisSep;
                pParentPfo = pPotentialParentPfo;
            }
        }

        if (!pParentPfo)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        childToParentPfoMap[pChildPfo] = pParentPfo;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPCheatHierarchyTool::GetClosestDistance(const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2) const
{
    ClusterList clusterList1, clusterList2;
    LArPfoHelper::GetThreeDClusterList(pPfo1, clusterList1);
    LArPfoHelper::GetThreeDClusterList(pPfo2, clusterList2);

    return LArClusterHelper::GetClosestDistance(clusterList1, clusterList2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject* MLPCheatHierarchyTool::BestParentInSplitHierarchy(const ParticleFlowObject *const pChildPfo, const PfoList &splitParticle) const
{
    float closestSep = std::numeric_limits<float>::max();
    const ParticleFlowObject *pParentPfo(nullptr);

    for (const ParticleFlowObject *const pPotentialParent : splitParticle)
    {
        float thisSep(this->GetClosestDistance(pPotentialParent, pPotentialParent));

        if (thisSep < closestSep)
        {
            closestSep = thisSep;
            pParentPfo = pPotentialParent;
        }
    }

    if (!pParentPfo)
    {
        std::cout << "No closest Pfo?" << std::endl; 
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    return pParentPfo;

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPCheatHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    
    return STATUS_CODE_SUCCESS;
}



//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
