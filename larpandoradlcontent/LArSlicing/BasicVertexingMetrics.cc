/**
 *  @file   src/BasicVertexingMetrics.cc
 *
 *  @brief  Implementation of the basic vertexing metrics.
 *
 *  $Log: $
 */

#include <cmath>
#include <map>
#include <string>

#include "Objects/CartesianVector.h"
#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArSlicing/BasicVertexingMetrics.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

BasicVertexingMetrics::BasicVertexingMetrics() :
    m_caloHitListName(""),
    m_inputVertexListName(""),
    m_treeName{"threeDVertexMetrics"},
    m_fileName{"threeDVertexMetrics.root"}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

BasicVertexingMetrics::~BasicVertexingMetrics()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "BasicVertexingMetrics: Unable to write to ROOT tree" << std::endl;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode BasicVertexingMetrics::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const VertexList *pInputVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pInputVertexList));

    if (nullptr == pCaloHitList || pCaloHitList->empty())
    {
        std::cout << "BasicVertexingMetrics::WriteMetrics - No CaloHits in this slice." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    if (nullptr == pInputVertexList || pInputVertexList->empty())
    {
        std::cout << "BasicVertexingMetrics::WriteMetrics - No vertices in this slice." << std::endl;
        return STATUS_CODE_SUCCESS;
    }


    std::map<const MCParticle *, CaloHitList> neutrinoToHitMap;
    std::map<const CaloHit*, const MCParticle*> hitToNeutrinoMap;
    unsigned int failureCount{0};

    // First, get a map of every hit to its neutrino, and the other way around.
    for (const auto pCaloHit : *pCaloHitList)
    {
        try {
            const MCParticle *const pMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);

            if (!pMCParticle)
                throw std::runtime_error("CaloHit with no associated MCParticle");

            const MCParticle *const pNeutrino(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

            if (nullptr == pNeutrino || !LArMCParticleHelper::IsNeutrino(pNeutrino))
                continue;

            neutrinoToHitMap[pNeutrino].push_back(pCaloHit);
            hitToNeutrinoMap[pCaloHit] = pNeutrino;
        } catch (StatusCodeException &e) {
            ++failureCount;
            continue;
        }
    }

    std::cout << "BasicVertexingMetrics::WriteMetrics - Found " << neutrinoToHitMap.size() << " neutrinos in this slice." << std::endl;
    std::cout << "BasicVertexingMetrics::WriteMetrics - Failed to find neutrino for " << failureCount << " hits, out of " << pCaloHitList->size() << std::endl;

    unsigned int mostHits{0};
    float mostEnergetic{0.f};

    for (const auto &neutrinoIt : neutrinoToHitMap)
    {
        const MCParticle *const pNeutrino(neutrinoIt.first);
        const CaloHitList allRecoHitsForNeutrino(neutrinoIt.second);

        if (allRecoHitsForNeutrino.size() > mostHits)
            mostHits = allRecoHitsForNeutrino.size();

        float totalEnergy{0.f};
        for (const auto &hit : allRecoHitsForNeutrino)
            totalEnergy += hit->GetInputEnergy();

        if (totalEnergy > mostEnergetic)
            mostEnergetic = totalEnergy;
    }

    // With that, we can start to write over the found vertices, and compare them to the MC truth.
    for (const auto &foundVertex : *pInputVertexList)
    {
        const CartesianVector &foundVertexPos = foundVertex->GetPosition();

        float closestDistance{std::numeric_limits<float>::max()};
        for (const auto &neutrinoIt : neutrinoToHitMap)
        {
            const MCParticle *const pNeutrino(neutrinoIt.first);
            const auto &trueVertex(pNeutrino->GetVertex());
            const float distSq = foundVertexPos.GetDistanceSquared(trueVertex);
            const float dist = std::sqrt(distSq);

            if (dist < closestDistance)
                closestDistance = dist;
        }

        for (const auto &neutrinoIt : neutrinoToHitMap)
        {
            const MCParticle *const pNeutrino(neutrinoIt.first);
            const CaloHitList allRecoHitsForNeutrino(neutrinoIt.second);

            const bool hasMostHits(allRecoHitsForNeutrino.size() >= mostHits);

            float totalEnergy{0.f};
            for (const auto &hit : allRecoHitsForNeutrino)
                totalEnergy += hit->GetInputEnergy();

            const bool hasMostEnergy(totalEnergy >= mostEnergetic);

            const float dist = std::sqrt(foundVertexPos.GetDistanceSquared(pNeutrino->GetVertex()));
            const bool isClosestVertex(dist <= closestDistance);

            const auto parentInstance = MultiPandoraApi::GetPrimaryPandoraInstance(&(*this).GetPandora());

            // Now, we can write out the reco and MC truth values for this vertex.
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "run", static_cast<int>(parentInstance->GetRun())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "subrun", static_cast<int>(parentInstance->GetSubrun())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "event", static_cast<int>(parentInstance->GetEvent())));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuEnergy", pNeutrino->GetEnergy()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuPDG", pNeutrino->GetParticleId()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNumNusInSlice", static_cast<int>(neutrinoToHitMap.size())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxX", pNeutrino->GetVertex().GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxY", pNeutrino->GetVertex().GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxZ", pNeutrino->GetVertex().GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNumHitsInSlice", static_cast<int>(allRecoHitsForNeutrino.size())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "biggestContributor", static_cast<int>(hasMostHits)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mostEnergetic", static_cast<int>(hasMostEnergy)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "closestVertex", static_cast<int>(isClosestVertex)));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVtxX", foundVertexPos.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVtxY", foundVertexPos.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVtxZ", foundVertexPos.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoNumSliceNumHits", static_cast<int>(pCaloHitList->size())));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxDx", foundVertexPos.GetX() - pNeutrino->GetVertex().GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxDy", foundVertexPos.GetY() - pNeutrino->GetVertex().GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxDz", foundVertexPos.GetZ() - pNeutrino->GetVertex().GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxDr", dist));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode BasicVertexingMetrics::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
