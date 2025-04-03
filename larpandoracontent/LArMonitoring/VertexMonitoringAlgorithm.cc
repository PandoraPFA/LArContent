/**
 *  @file   larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"
#include "larpandoracontent/LArObjects/LArEventTopology.h"

using namespace pandora;

namespace lar_content
{

VertexMonitoringAlgorithm::VertexMonitoringAlgorithm() :
    m_visualise{true},
    m_writeFile{false},
    m_useSecondaries{false},
    m_transparencyThresholdE{-1.f},
    m_energyScaleThresholdE{1.f},
    m_scalingFactor{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexMonitoringAlgorithm::~VertexMonitoringAlgorithm()
{
    if (m_writeFile)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename, m_filename, "UPDATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "VertexMonitoringAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::Run()
{
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(
            this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));
    }

    if (!m_useSecondaries)
        this->AssessVertices();
    else
        this->AssessSecondaryVertices();

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::AssessVertices() const
{
#ifdef MONITORING
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    const MCParticle *pTrueNeutrino{nullptr};
    const ParticleFlowObject *pRecoNeutrino{nullptr};
    if (!primaries.empty())
    {
        for (const MCParticle *primary : primaries)
        {
            const MCParticleList &parents{primary->GetParentList()};
            if (parents.size() == 1 && LArMCParticleHelper::IsNeutrino(parents.front()))
            {
                pTrueNeutrino = parents.front();
                break;
            }
        }
    }

    for (const ParticleFlowObject *pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            pRecoNeutrino = pPfo;
            break;
        }
    }

    MCParticleList primariesList(primaries.begin(), primaries.end());
    const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primariesList)};

    if (pRecoNeutrino && pTrueNeutrino)
    {
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};
        const CartesianVector &recoVertex{LArPfoHelper::GetVertex(pRecoNeutrino)->GetPosition()};
        if (m_visualise)
        {
            const CartesianVector tu(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ())));
            const CartesianVector tv(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ())));
            const CartesianVector tw(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ())));

            const CartesianVector ru(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(recoVertex.GetY(), recoVertex.GetZ())));
            const CartesianVector rv(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(recoVertex.GetY(), recoVertex.GetZ())));
            const CartesianVector rw(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(recoVertex.GetY(), recoVertex.GetZ())));

            const float du{(ru - tu).GetMagnitude()};
            const float dv{(rv - tv).GetMagnitude()};
            const float dw{(rw - tw).GetMagnitude()};

            std::cout << "delta(u, v, w): (" << du << ", " << dv << "," << dw << ")" << std::endl;

            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tu, "U true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tv, "V true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tw, "W true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &ru, "U reco vertex", RED, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &rv, "V reco vertex", RED, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &rw, "W reco vertex", RED, 2));
        }

        if (m_writeFile && LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, "dune_fd_hd"))
        {
            const CartesianVector delta{recoVertex - trueVertex};
            const float dx{delta.GetX()}, dy{delta.GetY()}, dz{delta.GetZ()}, dr{delta.GetMagnitude()};
            const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
            const int success{1};
            const int isCC{descriptor.IsCC()};
            const int isQE{descriptor.IsQE()};
            const int isRes{descriptor.IsResonant()};
            const int isDIS{descriptor.IsDIS()};
            const int isCoh{descriptor.IsCoherent()};
            const int isOther{!(isCC || isQE || isRes || isDIS)};
            const int isMu{descriptor.IsMuonNeutrino()};
            const int isElectron{descriptor.IsElectronNeutrino()};
            const int nPiZero{static_cast<int>(descriptor.GetNumPiZero())};
            const int nPiPlus{static_cast<int>(descriptor.GetNumPiPlus())};
            const int nPiMinus{static_cast<int>(descriptor.GetNumPiMinus())};
            const int nPhotons{static_cast<int>(descriptor.GetNumPhotons())};
            const int nProtons{static_cast<int>(descriptor.GetNumProtons())};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCC", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isQE", isQE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isRes", isRes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isDIS", isDIS));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCoh", isCoh));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isOther", isOther));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMu", isMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isElectron", isElectron));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiZero", nPiZero));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiPlus", nPiPlus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiMinus", nPiMinus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPhotons", nPhotons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nProtons", nProtons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", dx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", dy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", dz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", dr));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }
    else if (pTrueNeutrino)
    {
        const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};

        if (m_writeFile && LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, "dune_fd_hd"))
        {
            const int success{0};
            const float dx{-999.f}, dy{-999.f}, dz{-999.f}, dr{-999.f};
            const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
            const int isCC{descriptor.IsCC()};
            const int isQE{descriptor.IsQE()};
            const int isRes{descriptor.IsResonant()};
            const int isDIS{descriptor.IsDIS()};
            const int isCoh{descriptor.IsCoherent()};
            const int isOther{!(isCC || isQE || isRes || isDIS)};
            const int isMu{descriptor.IsMuonNeutrino()};
            const int isElectron{descriptor.IsElectronNeutrino()};
            const int nPiZero{static_cast<int>(descriptor.GetNumPiZero())};
            const int nPiPlus{static_cast<int>(descriptor.GetNumPiPlus())};
            const int nPiMinus{static_cast<int>(descriptor.GetNumPiMinus())};
            const int nPhotons{static_cast<int>(descriptor.GetNumPhotons())};
            const int nProtons{static_cast<int>(descriptor.GetNumProtons())};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCC", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isQE", isQE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isRes", isRes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isDIS", isDIS));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCoh", isCoh));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isOther", isOther));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMu", isMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isElectron", isElectron));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiZero", nPiZero));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiPlus", nPiPlus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiMinus", nPiMinus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPhotons", nPhotons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nProtons", nProtons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", dx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", dy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", dz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", dr));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }
#endif

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::AssessSecondaryVertices() const
{
#ifdef MONITORING
    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    LArEventTopology eventTopology(*pCaloHitList2D);
    eventTopology.ConstructVisibleHierarchy();
    //eventTopology.PruneHierarchy();
    CartesianPointVector trueVertices;
    eventTopology.GetVertices(trueVertices);
    if (trueVertices.empty())
        return STATUS_CODE_SUCCESS;

    const VertexList *pRecoVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_secVertexListName, pRecoVertexList));
    CartesianPointVector recoVertices;
    for (const Vertex *const pVertex : *pRecoVertexList)
        recoVertices.emplace_back(pVertex->GetPosition());

    std::vector<std::pair<CartesianVector, CartesianVector>> matches;
    for (const CartesianVector &recoVertex : recoVertices)
    {
        float best{std::numeric_limits<float>::max()};
        matches.emplace_back(std::make_pair(recoVertex, trueVertices.front()));
        for (const CartesianVector &trueVertex : trueVertices)
        {
            const float current{(recoVertex - trueVertex).GetMagnitudeSquared()};
            if (current < best)
            {
                best = current;
                matches.back().second = trueVertex;
            }
        }

        if (m_writeFile)
        {
            const CartesianVector delta{recoVertex - matches.back().second};
            const float dx{delta.GetX()}, dy{delta.GetY()}, dz{delta.GetZ()}, dr{delta.GetMagnitude()};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", dx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", dy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", dz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", dr));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }

        if (m_visualise)
        {
            std::cout << "(" << matches.back().first.GetX() << "," << matches.back().first.GetY() << "," << matches.back().first.GetZ() << ")"
                      << " == " << matches.back().second.GetX() << "," << matches.back().second.GetY() << ","
                      << matches.back().second.GetZ() << std::endl;
            const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
            const CartesianVector ur(matches.back().first.GetX(), 0.f,
                static_cast<float>(transform->YZtoU(matches.back().first.GetY(), matches.back().first.GetZ())));
            const CartesianVector vr(matches.back().first.GetX(), 0.f,
                static_cast<float>(transform->YZtoV(matches.back().first.GetY(), matches.back().first.GetZ())));
            const CartesianVector wr(matches.back().first.GetX(), 0.f,
                static_cast<float>(transform->YZtoW(matches.back().second.GetY(), matches.back().first.GetZ())));
            const CartesianVector ut(matches.back().second.GetX(), 0.f,
                static_cast<float>(transform->YZtoU(matches.back().second.GetY(), matches.back().second.GetZ())));
            const CartesianVector vt(matches.back().second.GetX(), 0.f,
                static_cast<float>(transform->YZtoV(matches.back().second.GetY(), matches.back().second.GetZ())));
            const CartesianVector wt(matches.back().second.GetX(), 0.f,
                static_cast<float>(transform->YZtoW(matches.back().second.GetY(), matches.back().second.GetZ())));

            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &ur, "Ur", RED, 1));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vr, "Vr", RED, 1));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &wr, "Wr", RED, 1));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &ut, "Ut", BLUE, 1));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vt, "Vt", BLUE, 1));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &wt, "Wt", BLUE, 1));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
#endif

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseSecondaries", m_useSecondaries));

    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Treename", m_treename));
    }
    if (m_useSecondaries)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SecondaryVertexListName", m_secVertexListName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
