/**
 *  @file   LArContent/src/LArMonitoring/ParticleAnalysisAlgorithm.cc
 *
 *  @brief  Implementation of the particle analysis algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPfoHelper.h"

#include "LArMonitoring/ParticleAnalysisAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ParticleAnalysisAlgorithm::ParticleAnalysisAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParticleAnalysisAlgorithm::~ParticleAnalysisAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleAnalysisAlgorithm::Run()
{
#ifdef MONITORING
    // Tree elements
    IntVector indexVector;
    IntVector pdgCodeVector, primaryVector, neutrinoVector, finalStateVector, foundVertexVector;
    IntVector nClustersVector, nSpacePointsVector;
    FloatVector vtxXPosVector, vtxYPosVector, vtxZPosVector, vtxPxVector, vtxPyVector, vtxPzVector, vtxPTotVector;

    try
    {
        // Load List of Pfos
        const PfoList *pPfoList(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

        PfoVector connectedPfos;
        this->GetConnectedPfos(pPfoList, connectedPfos);

        if (connectedPfos.empty())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        // Loop over Pfos and fill vectors
        int index(0);

        for (PfoVector::const_iterator pIter = connectedPfos.begin(), pIterEnd = connectedPfos.end(); pIter != pIterEnd; ++pIter)
        {
            const ParticleFlowObject *const pPfo = *pIter;

            int pdgCode(0), primary(0), neutrino(0), finalState(0), foundVertex(0);
            int nClusters(0), nSpacePoints(0);
            float vtxXPos(0.f), vtxYPos(0.f), vtxZPos(0.f), vtxPx(0.f), vtxPy(0.f), vtxPz(0.f), vtxPTot(0.f);

            pdgCode = pPfo->GetParticleId();
            primary = (pPfo->GetParentPfoList().empty() ? 1 : 0);
            neutrino = LArPfoHelper::IsNeutrino(pPfo);
            finalState = LArPfoHelper::IsFinalState(pPfo);

            vtxPx = pPfo->GetMomentum().GetX();
            vtxPy = pPfo->GetMomentum().GetY();
            vtxPz = pPfo->GetMomentum().GetZ();
            vtxPTot = pPfo->GetMomentum().GetMagnitude();

            if (pPfo->GetVertexList().size() > 0)
            {
                if (pPfo->GetVertexList().size() != 1)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const Vertex *const pVertex = *(pPfo->GetVertexList().begin());

                foundVertex = 1;
                vtxXPos = pVertex->GetPosition().GetX();
                vtxYPos = pVertex->GetPosition().GetY();
                vtxZPos = pVertex->GetPosition().GetZ();
            }

            ClusterList clusterList2D;
            CaloHitList caloHitList3D;

            LArPfoHelper::GetClusters(pPfo, TPC_VIEW_U, clusterList2D);
            LArPfoHelper::GetClusters(pPfo, TPC_VIEW_V, clusterList2D);
            LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, clusterList2D);
            LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList3D);

            nClusters = clusterList2D.size();
            nSpacePoints = caloHitList3D.size();

            indexVector.push_back(index++);
            pdgCodeVector.push_back(pdgCode);
            primaryVector.push_back(primary);
            neutrinoVector.push_back(neutrino);
            finalStateVector.push_back(finalState);
            foundVertexVector.push_back(foundVertex);
            nClustersVector.push_back(nClusters);
            nSpacePointsVector.push_back(nSpacePoints);
            vtxXPosVector.push_back(vtxXPos);
            vtxYPosVector.push_back(vtxYPos);
            vtxZPosVector.push_back(vtxZPos);
            vtxPxVector.push_back(vtxPx);
            vtxPyVector.push_back(vtxPy);
            vtxPzVector.push_back(vtxPz);
            vtxPTotVector.push_back(vtxPTot);
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "index", &indexVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pdgcode", &pdgCodeVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "primary", &primaryVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "neutrino", &neutrinoVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "finalstate", &finalStateVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertex", &foundVertexVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "clusters", &nClustersVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "spacepoints", &nSpacePointsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxx", &vtxXPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxy", &vtxYPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxz", &vtxZPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "px", &vtxPxVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "py", &vtxPyVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pz", &vtxPzVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ptot", &vtxPTotVector));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
#endif
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleAnalysisAlgorithm::GetConnectedPfos(const PfoList *const pPfoList, PfoVector &pfoVector) const
{
    PfoList connectedPfos;

    for (PfoList::iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        ParticleFlowObject *pPfo = *pIter;
        LArPfoHelper::GetAllDownstreamPfos(pPfo, connectedPfos);
    }

    for (PfoList::iterator pIter = connectedPfos.begin(), pIterEnd = connectedPfos.end(); pIter != pIterEnd; ++pIter)
        pfoVector.push_back(*pIter);

    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleAnalysisAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
