/**
 *  @file   larpandoracontent/LArMonitoring/VertexAssessmentAlgorithm.cc
 *
 *  @brief  Implementation of the vertex assessment algorithm class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VertexAssessmentAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar_content
{

VertexAssessmentAlgorithm::VertexAssessmentAlgorithm() :
    m_visualize{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexAssessmentAlgorithm::~VertexAssessmentAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "vertex", "vertex.root", "RECREATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "VertexAssessmentAlgorithm: Unable to write to ROOT tree" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssessmentAlgorithm::Run()
{
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }
    this->Visualize();  // TODO: More generic name, as it does more than just visualize
    if (m_visualize)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexAssessmentAlgorithm::Visualize() const
{
    static int event{-1};
    ++event;
    const CaloHitList *pCaloHitListU{nullptr}, *pCaloHitListV{nullptr}, *pCaloHitListW{nullptr};
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, "CaloHitListU", pCaloHitListU))
        return;
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, "CaloHitListV", pCaloHitListV))
        return;
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitListW))
        return;
    if (m_visualize)
    {
        if (pCaloHitListU && !pCaloHitListU->empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitListU, "hits_u", BLACK));
        }
        if (pCaloHitListV && !pCaloHitListV->empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitListV, "hits_v", BLACK));
        }
        if (pCaloHitListW && !pCaloHitListW->empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitListW, "hits_w", BLACK));
        }
    }

    try
    {
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const MCParticleList *pMCParticleList{nullptr};
        CartesianVector trueVertex(0.f, 0.f, 0.f);
        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetCurrentList(*this, pMCParticleList))
        {
            if (pMCParticleList)
            {
                LArMCParticleHelper::MCContributionMap mcToHitsMap;
                this->GetMCToHitsMap(mcToHitsMap);

                MCParticleVector primaries;
                LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
                if (!primaries.empty())
                {
                    const MCParticle *primary{primaries.front()};
                    const MCParticleList &parents{primary->GetParentList()};
                    if (parents.size() == 1)
                    {
                        const MCParticle *trueNeutrino{parents.front()};
                        if (LArMCParticleHelper::IsNeutrino(trueNeutrino))
                        {
                            trueVertex = primaries.front()->GetVertex();
                            if (m_visualize)
                            {
                                const CartesianVector trueU(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ())));
                                const CartesianVector trueV(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ())));
                                const CartesianVector trueW(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ())));
                                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueU, "true_vtx_u", ORANGE, 2));
                                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueV, "true_vtx_v", ORANGE, 2));
                                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueW, "true_vtx_w", ORANGE, 2));
                            }
                        }
                    }
                }
            }
        }
    }
    catch (const StatusCodeException &)
    {   // No MC available
    }

    const VertexList *pVertexList{nullptr};
    if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_vertexListName, pVertexList))
    {
        const MCParticleList *pMCParticleList{nullptr};
        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetCurrentList(*this, pMCParticleList))
        {
            if (pMCParticleList)
            {
                LArMCParticleHelper::MCContributionMap mcToHitsMap;
                MCParticleVector primaries;
                LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
                if (!primaries.empty())
                {
                    const MCParticle *primary{primaries.front()};
                    const MCParticleList &parents{primary->GetParentList()};
                    if (parents.size() == 1)
                    {
                        const MCParticle *trueNeutrino{parents.front()};
                        if (LArMCParticleHelper::IsNeutrino(trueNeutrino))
                        {
                            const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
                            const CartesianVector &trueVertex{primaries.front()->GetVertex()};
                            if (LArVertexHelper::IsInFiducialVolume(trueVertex, "dune_fd_hd"))
                            {
                                if (pVertexList && !pVertexList->empty())
                                {
                                    for (const Vertex *pVertex : *pVertexList)
                                    {
                                        const CartesianVector &recoVertex(pVertex->GetPosition());
                                        if (pVertex->GetVertexLabel() == VERTEX_INTERACTION)
                                        {
                                            if (m_visualize)
                                            {
                                                const CartesianVector recoU(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(
                                                    recoVertex.GetY(), recoVertex.GetZ())));
                                                const CartesianVector recoV(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(
                                                    recoVertex.GetY(), recoVertex.GetZ())));
                                                const CartesianVector recoW(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(
                                                    recoVertex.GetY(), recoVertex.GetZ())));
                                                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &recoU, "reco_vtx_u", RED, 2));
                                                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &recoV, "reco_vtx_v", GREEN, 2));
                                                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &recoW, "reco_vtx_w", BLUE, 2));
                                            }

                                            const CartesianVector &dv{recoVertex - trueVertex};
                                            const float dr{dv.GetMagnitude()};
                                            const float dx{dv.GetX()}, dy{dv.GetY()}, dz{dv.GetZ()};
                                            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "vertex", "event", event));
                                            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "vertex", "pass", 1));
                                            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "vertex", "dr", dr));
                                            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "vertex", "dx", dx));
                                            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "vertex", "dy", dy));
                                            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "vertex", "dz", dz));
                                            PANDORA_MONITORING_API(FillTree(this->GetPandora(), "vertex"));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssessmentAlgorithm::GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList2D, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssessmentAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
