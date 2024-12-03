/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPNeutrinoHierarchyAlgorithm.cc
 *
 *  @brief  Implementation of the MLP neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPNeutrinoHierarchyAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::HierarchyPfo() :
    m_pPfo(nullptr),
    m_pParentPfo(nullptr),
    m_pPredictedParentPfo(nullptr),
    m_primaryScore(-std::numeric_limits<float>::max()),
    m_laterTierScore(-std::numeric_limits<float>::max()),
    m_parentOrientation(-1),
    m_childOrientation(-1),
    m_isSet(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::HierarchyPfo(const ParticleFlowObject *pPfo) :
    HierarchyPfo()
{
    m_pPfo = pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::operator== (const HierarchyPfo &otherHierarchyPfo) const
{
    return this->GetPfo() == otherHierarchyPfo.GetPfo();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MLPNeutrinoHierarchyAlgorithm::MLPNeutrinoHierarchyAlgorithm() :
    m_primaryThresholdTrackPass1(0.8f),
    m_primaryThresholdShowerPass1(0.45f),
    m_laterTierThresholdTrackPass1(0.8f),
    m_laterTierThresholdShowerPass1(0.8f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPNeutrinoHierarchyAlgorithm::Run()
{
    srand(0.1);
    std::cout << "RUNNING THE NEW HIERARCHY ALG!" << std::endl;

    // Fill the track/shower vectors
    this->FillTrackShowerVectors();

    // Calculate primary scores
    this->SetPrimaryScores();

    // Set what primary things we can
    this->BuildPrimaryTierPass1();

    if (m_hierarchy.size() == 0)
    {
        std::cout << "NO PRIMARIES FOUND!" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    std::cout << "nPrimaries: " << m_hierarchy.at(0).size() << std::endl;

    // Set later tier scores
    this->SetLaterTierScores();

    // Build the later tier
    this->BuildLaterTierPass1();


    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::FillTrackShowerVectors()
{
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject * pPfo : *pPfoList)
        {
            // Apply hit cut
            ClusterList clusterList3D;
            LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

            int total3DHits(0);

            for (const Cluster *const pCluster3D : clusterList3D)
                total3DHits += pCluster3D->GetNCaloHits();

            if (total3DHits == 0)
                continue;

            // Create track/shower objects
            if (pPfo->GetParticleId() == 13)
                m_trackPfos[pPfo] = HierarchyPfo(pPfo);
            else if (pPfo->GetParticleId() == 11) 
                m_showerPfos[pPfo] = HierarchyPfo(pPfo);
            else
                std::cout << "IDK what this pfo is" << std::endl;
        }
    }

    std::cout << "We have " << m_trackPfos.size() << " track(s) and " << m_showerPfos.size() << " shower(s)" << std::endl;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::SetPrimaryScores()
{
    // For tracks
    for (auto& [pPfo, hierarchyPfo] : m_trackPfos)
        this->SetPrimaryScoreTrack(hierarchyPfo);

    // For showers
    for (auto& [pPfo, hierarchyPfo] : m_showerPfos)
        this->SetPrimaryScoreTrack(hierarchyPfo);

    //////////////////////////////////////////////////////////////////
    for (const auto& [pPfo, hierarchyPfo] : m_trackPfos)
        std::cout << "track primary score: " << hierarchyPfo.GetPrimaryScore() << std::endl;

    for (const auto& [pPfo, hierarchyPfo] : m_showerPfos)
        std::cout << "shower primary score: " << hierarchyPfo.GetPrimaryScore() << std::endl;
    //////////////////////////////////////////////////////////////////    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::SetPrimaryScoreTrack(HierarchyPfo &trackPfo)
{
    trackPfo.SetPrimaryScore(this->GetRandomNumber());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::SetPrimaryScoreShower(HierarchyPfo &showerPfo)
{
    showerPfo.SetPrimaryScore(this->GetRandomNumber());
}


//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::BuildPrimaryTierPass1()
{
    std::vector<const ParticleFlowObject*> primaryTier;

    for (bool isTrack : {true, false})
    {
        std::map<const ParticleFlowObject*, HierarchyPfo> &hierarchyPfoMap(isTrack ? m_trackPfos : m_showerPfos);

        for (auto& [pPfo, hierarchyPfo] : hierarchyPfoMap)
        {
            const float thisPrimaryScore(hierarchyPfo.GetPrimaryScore());

            std::cout << "thisPrimaryScore: " << thisPrimaryScore << std::endl;

            if ((thisPrimaryScore > 0.f) && ((isTrack && (thisPrimaryScore > m_primaryThresholdTrackPass1)) || 
                                           (!isTrack && (thisPrimaryScore > m_primaryThresholdShowerPass1))))
            {
                hierarchyPfo.SetIsSet(true);
                primaryTier.emplace_back(pPfo);
            }
        }
    }

    m_hierarchy.emplace_back(primaryTier);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::SetLaterTierScores()
{
    for (bool isTrack : {true, false})
    {
        std::map<const ParticleFlowObject*, HierarchyPfo> &hierarchyPfoMap(isTrack ? m_trackPfos : m_showerPfos);

        for (auto& [pChildPfo, childPfo] : hierarchyPfoMap)
        {
            // Don't bother if we have already found a parent
            if (childPfo.GetIsSet())
                continue;

            float highestScore(-std::numeric_limits<float>::max());

            for (const auto& [pParentPfo, parentPfo] : m_trackPfos)
            {
                if (pChildPfo == pParentPfo)
                    continue;

                int parentOrientation(-1), childOrientation(-1);
                const float thisScore(isTrack ? this->GetLaterTierScoreTrackToTrack(parentPfo, childPfo, parentOrientation, childOrientation) :
                    this->GetLaterTierScoreTrackToShower(parentPfo, childPfo, parentOrientation, childOrientation));

                if (thisScore > highestScore)
                {
                    childPfo.SetLaterTierScore(thisScore);
                    childPfo.SetParentOrientation(parentOrientation);
                    childPfo.SetChildOrientation(childOrientation);
                    childPfo.SetPredictedParentPfo(parentPfo.GetPfo());
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPNeutrinoHierarchyAlgorithm::BuildLaterTierPass1()
{
    while (m_hierarchy.at(m_hierarchy.size() - 1).size() != 0)
    {
        std::vector<const ParticleFlowObject*> &preceedingTier(m_hierarchy.at(m_hierarchy.size() - 1));
        std::vector<const ParticleFlowObject*> thisTier;

        for (bool isTrack : {true, false})
        {
            std::map<const ParticleFlowObject*, HierarchyPfo> &hierarchyPfoMap(isTrack ? m_trackPfos : m_showerPfos);

            for (auto& [pChildPfo, childPfo] : hierarchyPfoMap)
            {
                // Don't bother if we have already found a parent
                if (childPfo.GetIsSet())
                    continue;

                // Do we have a predicted parent
                if (!childPfo.GetPredictedParentPfo())
                    continue;

                // Is predicted parent in preceeding tier?
                const ParticleFlowObject *const pPredictedParent = childPfo.GetPredictedParentPfo();
                if (std::find(preceedingTier.begin(), preceedingTier.end(), pPredictedParent) == preceedingTier.end())
                    continue;

                // Does it pass tier cut?
                if ((isTrack && (childPfo.GetLaterTierScore() > m_laterTierThresholdTrackPass1)) ||
                    (!isTrack && (childPfo.GetLaterTierScore() > m_laterTierThresholdShowerPass1)))
                {
                    thisTier.emplace_back(pChildPfo);
                }
            }
        }

        m_hierarchy.emplace_back(thisTier);
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetLaterTierScoreTrackToTrack(const HierarchyPfo &/*parentPfo*/, const HierarchyPfo &/*childPfo*/, 
    int &parentOrientation, int &childOrientation) const
{
    return this->GetRandomNumber();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetLaterTierScoreTrackToShower(const HierarchyPfo &/*parentPfo*/, const HierarchyPfo &/*childPfo*/,
    int &parentOrientation, int &childOrientation) const
{
    return this->GetRandomNumber();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MLPNeutrinoHierarchyAlgorithm::GetRandomNumber() const
{
    int randomNumber1 = rand() % 10 + 1;
    int randomNumber2 = rand() % 10 + 1;

    return ((static_cast<float>(randomNumber1) / 10.f) + (static_cast<float>(randomNumber2) / 100.f));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPNeutrinoHierarchyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryThresholdTrackPass1", m_primaryThresholdTrackPass1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryThresholdShowerPass1", m_primaryThresholdShowerPass1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierThresholdTrackPass1", m_laterTierThresholdTrackPass1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LaterTierThresholdShowerPass1", m_laterTierThresholdShowerPass1));

    return STATUS_CODE_SUCCESS;
}







} // namespace lar_dl_content
