/**
 *  @file   larpandoracontent/LArThreeReco/LArEventBuilding/PfoCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the pfo characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/PfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoCharacterisationAlgorithm::PfoCharacterisationAlgorithm() :
    m_writeToTree(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoCharacterisationAlgorithm::~PfoCharacterisationAlgorithm()
{
    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationAlgorithm::Run()
{
    PfoList tracksToShowers, showersToTracks;

    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PfoCharacterisationAlgorithm: unable to find pfo list " << pfoListName << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            if (this->IsClearTrack(pPfo))
            {
                PandoraContentApi::ParticleFlowObject::Metadata metadata;
                metadata.m_particleId = MU_MINUS;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));

                if (m_showerPfoListName == pfoListName)
                    showersToTracks.push_back(pPfo);
            }
            else
            {
                PandoraContentApi::ParticleFlowObject::Metadata metadata;
                metadata.m_particleId = E_MINUS;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));

                if (m_trackPfoListName == pfoListName)
                    tracksToShowers.push_back(pPfo);
            }
        }
    }

    if (!tracksToShowers.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

    if (!showersToTracks.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCharacterisationAlgorithm::IsClearTrack(const ParticleFlowObject *const pPfo) const
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList);

    MCParticleWeightMap mcParticleWeightMap;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        for (const MCParticleWeightMap::value_type &mapEntry : pCaloHit->GetMCParticleWeightMap())
            mcParticleWeightMap[mapEntry.first] += mapEntry.second;
    }

    float bestWeight(0.f);
    const MCParticle *pBestMCParticle(nullptr);

    MCParticleList mcParticleList;
    for (const auto &mapEntry : mcParticleWeightMap) mcParticleList.push_back(mapEntry.first);
    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        const float weight(mcParticleWeightMap.at(pMCParticle));

        if (weight > bestWeight)
        {
            pBestMCParticle = pMCParticle;
            bestWeight = weight;
        }
    }

    bool isTrueTrack(false);

    if (pBestMCParticle)
    {
        const MCParticle *const pPrimaryMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(pBestMCParticle));
        const int absParticleId(std::abs(pPrimaryMCParticle->GetParticleId()));
        isTrueTrack = ((MU_MINUS == absParticleId) || (PROTON == absParticleId) || (PI_PLUS == absParticleId));
    }

    // Tree variables here
    //--------------------------------------------------------------------------------------------------------------------------------------
    const int trueTrack(isTrueTrack ? 1 : 0);
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueTrack", trueTrack));

    const int inputParticleId(pPfo->GetParticleId());
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "inputParticleId", inputParticleId));

    int nHitsU(0), nHitsV(0), nHitsW(0);
    int nClustersU(0), nClustersV(0), nClustersW(0);
    float vertexDistanceU(-1.f), vertexDistanceV(-1.f), vertexDistanceW(-1.f);

    ClusterList twoDClusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

    for (const Cluster *const pCluster : twoDClusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
            continue;

        const std::string hitTypeLabel((TPC_VIEW_U == hitType) ? "U" : (TPC_VIEW_V == hitType) ? "V" : "W");

        // nHits
        int &nHits((TPC_VIEW_U == hitType) ? nHitsU : (TPC_VIEW_V == hitType) ? nHitsV : nHitsW);
        nHits += pCluster->GetNCaloHits();
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits" + hitTypeLabel, nHits));

        // nClusters
        int &nClusters((TPC_VIEW_U == hitType) ? nClustersU : (TPC_VIEW_V == hitType) ? nClustersV : nClustersW);
        ++nClusters;
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nClusters" + hitTypeLabel, nClusters));

        // vertexDistance
        float &vertexDistance((TPC_VIEW_U == hitType) ? vertexDistanceU : (TPC_VIEW_V == hitType) ? vertexDistanceV : vertexDistanceW);

        const VertexList *pVertexList = nullptr;
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
        const Vertex *const pSelectedVertex((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

        if (pSelectedVertex)
        {
            const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), LArClusterHelper::GetClusterHitType(pCluster)));
            vertexDistance = LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster);
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertexDistance" + hitTypeLabel, vertexDistance));

        // Add more variables

    }

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    //--------------------------------------------------------------------------------------------------------------------------------------
    // End tree variables calculations

    return isTrueTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    m_inputPfoListNames.push_back(m_trackPfoListName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    m_inputPfoListNames.push_back(m_showerPfoListName);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
