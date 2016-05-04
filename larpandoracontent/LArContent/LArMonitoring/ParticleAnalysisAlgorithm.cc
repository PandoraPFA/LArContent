/**
 *  @file   LArContent/src/LArMonitoring/ParticleAnalysisAlgorithm.cc
 *
 *  @brief  Implementation of the particle analysis algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPfoHelper.h"

#include "LArObjects/LArTrackPfo.h"
#include "LArObjects/LArShowerPfo.h"

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
    IntVector pdgCodeVector, primaryVector, neutrinoVector, finalStateVector;
    IntVector foundVertexVector, foundTrackVector, foundShowerVector;
    IntVector nClustersVector, nSpacePointsVector, nCaloHitsVector, nIsolatedHitsVector;
    FloatVector pfoVtxXPosVector, pfoVtxYPosVector, pfoVtxZPosVector, pfoVtxXDirVector, pfoVtxYDirVector, pfoVtxZDirVector;
    FloatVector trkVtxXPosVector, trkVtxYPosVector, trkVtxZPosVector, trkVtxXDirVector, trkVtxYDirVector, trkVtxZDirVector;
    FloatVector trkEndXPosVector, trkEndYPosVector, trkEndZPosVector, trkEndXDirVector, trkEndYDirVector, trkEndZDirVector;

    try
    {
        // Load List of Pfos
        const PfoList *pPfoList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_pfoListName, 
            pPfoList));

        PfoVector connectedPfos;

        if (NULL != pPfoList)
            this->GetConnectedPfos(pPfoList, connectedPfos);

        if (connectedPfos.empty())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        // Loop over Pfos and fill vectors
        int index(0);

        for (PfoVector::const_iterator pIter = connectedPfos.begin(), pIterEnd = connectedPfos.end(); pIter != pIterEnd; ++pIter)
        {
            const ParticleFlowObject *const pPfo = *pIter;

            int pdgCode(0), primary(0), neutrino(0), finalState(0);
            int foundVertex(0), foundTrack(0), foundShower(0);
            int nClusters(0), nSpacePoints(0), nCaloHits(0), nIsolatedHits(0);
            float pfoVtxXPos(0.f), pfoVtxYPos(0.f), pfoVtxZPos(0.f), pfoVtxXDir(0.f), pfoVtxYDir(0.f), pfoVtxZDir(0.f);
            float trkVtxXPos(0.f), trkVtxYPos(0.f), trkVtxZPos(0.f), trkVtxXDir(0.f), trkVtxYDir(0.f), trkVtxZDir(0.f);
            float trkEndXPos(0.f), trkEndYPos(0.f), trkEndZPos(0.f), trkEndXDir(0.f), trkEndYDir(0.f), trkEndZDir(0.f);

            pdgCode = pPfo->GetParticleId();
            primary = (pPfo->GetParentPfoList().empty() ? 1 : 0);
            neutrino = LArPfoHelper::IsNeutrino(pPfo);
            finalState = LArPfoHelper::IsFinalState(pPfo);

            const float pTot(pPfo->GetMomentum().GetMagnitude());

            if (pTot > std::numeric_limits<float>::epsilon())
            {
                pfoVtxXDir = pPfo->GetMomentum().GetX() / pTot;
                pfoVtxYDir = pPfo->GetMomentum().GetY() / pTot;
                pfoVtxZDir = pPfo->GetMomentum().GetZ() / pTot;
            }

            if (!pPfo->GetVertexList().empty())
            {
                const Vertex *const pVertex = LArPfoHelper::GetVertex(pPfo);
                foundVertex = 1;
                pfoVtxXPos = pVertex->GetPosition().GetX();
                pfoVtxYPos = pVertex->GetPosition().GetY();
                pfoVtxZPos = pVertex->GetPosition().GetZ();
            }

            CaloHitList caloHitList2D, caloHitList3D, isolatedHitList2D;
            LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList3D);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitList2D);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_U, isolatedHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_V, isolatedHitList2D);
            LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_W, isolatedHitList2D);

            ClusterList clusterList2D;
            LArPfoHelper::GetTwoDClusterList(pPfo, clusterList2D);

            nClusters     = clusterList2D.size();
            nSpacePoints  = caloHitList3D.size();
            nCaloHits     = caloHitList2D.size();
            nIsolatedHits = isolatedHitList2D.size();

            const LArTrackPfo *const pLArTrackPfo = dynamic_cast<const LArTrackPfo*>(pPfo);
            if (NULL != pLArTrackPfo)
            {
                foundTrack = 1;
                trkVtxXPos = pLArTrackPfo->GetVertexPosition().GetX();
                trkVtxYPos = pLArTrackPfo->GetVertexPosition().GetY();
                trkVtxZPos = pLArTrackPfo->GetVertexPosition().GetZ();
                trkVtxXDir = pLArTrackPfo->GetVertexDirection().GetX();
                trkVtxYDir = pLArTrackPfo->GetVertexDirection().GetY();
                trkVtxZDir = pLArTrackPfo->GetVertexDirection().GetZ();
                trkEndXPos = pLArTrackPfo->GetEndPosition().GetX();
                trkEndYPos = pLArTrackPfo->GetEndPosition().GetY();
                trkEndZPos = pLArTrackPfo->GetEndPosition().GetZ();
                trkEndXDir = pLArTrackPfo->GetEndDirection().GetX();
                trkEndYDir = pLArTrackPfo->GetEndDirection().GetY();
                trkEndZDir = pLArTrackPfo->GetEndDirection().GetZ();
            }

            const LArShowerPfo *const pLArShowerPfo = dynamic_cast<const LArShowerPfo*>(pPfo);
            if (NULL != pLArShowerPfo)
            {
                foundShower = 1;
            }

            indexVector.push_back(index++);
            pdgCodeVector.push_back(pdgCode);
            primaryVector.push_back(primary);
            neutrinoVector.push_back(neutrino);
            finalStateVector.push_back(finalState);
            foundVertexVector.push_back(foundVertex);
            foundTrackVector.push_back(foundTrack);
            foundShowerVector.push_back(foundShower);
            nClustersVector.push_back(nClusters);
            nSpacePointsVector.push_back(nSpacePoints);
            nCaloHitsVector.push_back(nCaloHits);
            nIsolatedHitsVector.push_back(nIsolatedHits);
            pfoVtxXPosVector.push_back(pfoVtxXPos);
            pfoVtxYPosVector.push_back(pfoVtxYPos);
            pfoVtxZPosVector.push_back(pfoVtxZPos);
            pfoVtxXDirVector.push_back(pfoVtxXDir);
            pfoVtxYDirVector.push_back(pfoVtxYDir);
            pfoVtxZDirVector.push_back(pfoVtxZDir);
            trkVtxXPosVector.push_back(trkVtxXPos);
            trkVtxYPosVector.push_back(trkVtxYPos);
            trkVtxZPosVector.push_back(trkVtxZPos);
            trkVtxXDirVector.push_back(trkVtxXDir);
            trkVtxYDirVector.push_back(trkVtxYDir);
            trkVtxZDirVector.push_back(trkVtxZDir);
            trkEndXPosVector.push_back(trkEndXPos);
            trkEndYPosVector.push_back(trkEndYPos);
            trkEndZPosVector.push_back(trkEndZPos);
            trkEndXDirVector.push_back(trkEndXDir);
            trkEndYDirVector.push_back(trkEndYDir);
            trkEndZDirVector.push_back(trkEndZDir);
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
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "track", &foundTrackVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "shower", &foundShowerVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "clusters", &nClustersVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "spacepoints", &nSpacePointsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "calohits", &nCaloHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isolatedhits", &nIsolatedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfovtxx", &pfoVtxXPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfovtxy", &pfoVtxYPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfovtxz", &pfoVtxZPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfodirx", &pfoVtxXDirVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfodiry", &pfoVtxYDirVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfodirz", &pfoVtxZDirVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkvtxx", &trkVtxXPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkvtxy", &trkVtxYPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkvtxz", &trkVtxZPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkvtxdirx", &trkVtxXDirVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkvtxdiry", &trkVtxYDirVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkvtxdirz", &trkVtxZDirVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkendx", &trkEndXPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkendy", &trkEndYPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkendz", &trkEndZPosVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkenddirx", &trkEndXDirVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkenddiry", &trkEndYDirVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trkenddirz", &trkEndZDirVector));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
#endif
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleAnalysisAlgorithm::GetConnectedPfos(const PfoList *const pPfoList, PfoVector &pfoVector) const
{
    PfoList connectedPfos;

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pPfo = *pIter;
        LArPfoHelper::GetAllDownstreamPfos(pPfo, connectedPfos);
    }

    for (PfoList::const_iterator pIter = connectedPfos.begin(), pIterEnd = connectedPfos.end(); pIter != pIterEnd; ++pIter)
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
