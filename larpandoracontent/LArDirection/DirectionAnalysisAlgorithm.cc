/**
 *  @file   ExampleContent/src/ExampleAlgorithms/DirectionAnalysisAlgorithm.cc
 * 
 *  @brief  Implementation of the access lists algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

//#include "larpandoracontent/LArHelpers/LArSpaceChargeHelper.h"

#include "DirectionAnalysisAlgorithm.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

std::string mcTreeName("MC"), clusterTreeName("Cluster"), pfoTreeName("PFO"), hitTreeName("Hit"), vertexTreeName("Vertex");

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionAnalysisAlgorithm::DirectionAnalysisAlgorithm():
    m_targetParticlePDG(13),
    m_particleContained(true),
    m_cosmic(true),
    m_data(false),
    m_drawFit(false),
    m_fileIdentifier(0),
    m_eventNumber(-1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DirectionAnalysisAlgorithm::~DirectionAnalysisAlgorithm()
{
    if (m_writeToTree)
    {    
        if (!m_data)
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), mcTreeName.c_str(), m_fileName.c_str(), "UPDATE"));

        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), clusterTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), pfoTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), hitTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), vertexTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionAnalysisAlgorithm::Run()
{
    m_eventNumber++;

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));


    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));


    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    pandora::MCParticleVector mcParticleVector(pMCParticleList->begin(), pMCParticleList->end());
    pandora::ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());

    pandora::PfoList allPfos;
    LArPfoHelper::GetAllConnectedPfos(*pPfoList, allPfos);

    pandora::PfoVector pfoVector(allPfos.begin(), allPfos.end());

    if (pVertexList->size() != 1)
        std::cout << "TOO MANY VERTICES" << std::endl;

    const pandora::Vertex* const pNeutrinoVertex(pVertexList->front());

    pandora::MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    if (trueNeutrinos.size() != 1)
        return STATUS_CODE_SUCCESS;

    if (LArMCParticleHelper::GetNuanceCode((*(trueNeutrinos.begin()))) != 1001) 
        return STATUS_CODE_SUCCESS;

    //pandora::CaloHitVector hitVector(pCaloHitList->begin(), pCaloHitList->end());
    //for (auto pCaloHit : hitVector)
    //    std::cout << "Hit position: (" << pCaloHit->GetPositionVector().GetX() << "," << pCaloHit->GetPositionVector().GetY() << "," << pCaloHit->GetPositionVector().GetZ() << ")" << std::endl;
    
    // LArSpaceChargeHelper::Configure("/usera/jjd49/pandora_direction/PandoraPFA/LArContent-origin/vertex_direction/larpandoracontent/LArDirection/SCEoffsets_MicroBooNE_E273.root");

    if (!m_data)
        this->WriteMCInformation(pMCParticleList, pCaloHitList);

    this->WritePfoInformation(pfoVector);
    this->WriteClusterAndHitInformation(clusterVector);
    this->WriteVertexInformation(pMCParticleList, pCaloHitList, pNeutrinoVertex, pfoVector);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteMCInformation(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList)
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

    LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

    LArMCParticleHelper::MCContributionMap cosmicPrimaryMCParticles, neutrinoPrimaryMCParticles;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, cosmicPrimaryMCParticles);
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, neutrinoPrimaryMCParticles);

    int nCosmicPrimaryMCParticles(cosmicPrimaryMCParticles.size()), nNeutrinoPrimaryMCParticles(neutrinoPrimaryMCParticles.size());

    for (const auto entry : mcToTrueHitListMap)
        std::cout << entry.first << std::endl;

    for (const auto pMCParticle : *pMCParticleList)
    {
        float energy(pMCParticle->GetEnergy()), length((pMCParticle->GetVertex() - pMCParticle->GetEndpoint()).GetMagnitude());
        int numberHits(mcToTrueHitListMap.count(pMCParticle) != 0 ? mcToTrueHitListMap.at(pMCParticle).size() : 0);

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "Energy", energy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "Length", length));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "NumberHits", numberHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "NCosmicPrimaryMCParticles", nCosmicPrimaryMCParticles));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "NNeutrinoPrimaryMCParticles", nNeutrinoPrimaryMCParticles));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "FileIdentifier", m_fileIdentifier));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "EventNumber", m_eventNumber));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), mcTreeName.c_str()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::CheckEventType(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, pandora::PfoVector &pfoVector, int &nMuons, int &nProtons, int &nOthers, int &nPFOs)
{
    LArMCParticleHelper::MCContributionMap selectedMCParticles;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, selectedMCParticles);

    for (auto mcParticle : selectedMCParticles)
    {
        if (mcParticle.first->GetParticleId() == 13)
            ++nMuons;
        else if (mcParticle.first->GetParticleId() == 2212)
            ++nProtons;
        else
            ++nOthers;
    }

    for (auto pPfo : pfoVector)
    {
        if (!LArPfoHelper::IsFinalState(pPfo))
            continue;

        ++nPFOs;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteVertexInformation(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, const pandora::Vertex* const pVertex, pandora::PfoVector &pfoVector)
{
    int nMuons(0), nProtons(0), nOthers(0), nPFOs(0);
    this->CheckEventType(pMCParticleList, pCaloHitList, pfoVector, nMuons, nProtons, nOthers, nPFOs);

    float vertexDR(this->GetVertexDR(pMCParticleList, pVertex, false));
    float sccVertexDR(this->GetVertexDR(pMCParticleList, pVertex, true));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "nMuons", nMuons));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "nProtons", nProtons));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "nOthers", nOthers));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "nPFOs", nPFOs));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "VertexDR", vertexDR));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "SCCVertexDR", sccVertexDR));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "FileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), vertexTreeName.c_str(), "EventNumber", m_eventNumber));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), vertexTreeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DirectionAnalysisAlgorithm::GetVertexDR(const pandora::MCParticleList *pMCParticleList, const pandora::Vertex* const pVertex, bool enableSpaceChargeCorrection)
{
    CartesianVector trueVertexPosition(0.f, 0.f, 0.f);

    for (auto pMCParticle : *pMCParticleList)
    {
        if (pMCParticle->GetParticleId() == 13 && LArMCParticleHelper::IsPrimary(pMCParticle) && pMCParticle->GetVertex().GetY() != 0)
            trueVertexPosition = pMCParticle->GetVertex();
    }

    CartesianVector vertexCandidatePosition(pVertex->GetPosition());
    
    if (enableSpaceChargeCorrection)
      //vertexCandidatePosition = LArSpaceChargeHelper::GetSpaceChargeCorrectedPosition(vertexCandidatePosition);
      std::cout << "space charge correction here" << std::endl;

    return (vertexCandidatePosition - trueVertexPosition).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WritePfoInformation(pandora::PfoVector &pfoVector)
{
    for (const pandora::ParticleFlowObject *const pPfo : pfoVector)
    {
        try
        {
            if (!this->IsGoodPfo(pPfo))
                continue;

            const Cluster *const pCluster = this->GetTargetClusterFromPFO(pPfo); 
            TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetPfoDirection(pPfo);
            this->WriteToTree(pCluster, fitResult, pfoTreeName);

            if (m_drawFit)
            {
                std::cout << "Split applied: " << fitResult.GetSplitObject().GetSplitApplied() << std::endl;
                fitResult.DrawFit();
                PANDORA_MONITORING_API(Pause(this->GetPandora()));
            }
        }

        catch (...)
        {
            std::cout << "Skipping PFO..." << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IsGoodPfo(const pandora::ParticleFlowObject* pPfo)
{
    //const MCParticle* const pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
    const Cluster *const pCluster = this->GetTargetClusterFromPFO(pPfo); 
    const MCParticle* const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

    if (!LArPfoHelper::IsFinalState(pPfo))
        return false;

    if ((m_particleContained && !m_data && IsParticleContained(pMCParticle)) || (m_particleContained && m_data && !IsRecoParticleContained(pPfo)))
        return false;

    if (pMCParticle->GetParticleId() != m_targetParticlePDG)
        return false;

    if (m_cosmic && !LArMCParticleHelper::IsCosmicRay(pMCParticle))
        return false;

    if (!m_cosmic && LArMCParticleHelper::IsCosmicRay(pMCParticle)) 
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteClusterAndHitInformation(pandora::ClusterVector &clusterVector)
{
    for (const pandora::Cluster* const pCluster : clusterVector)
    {
        try
        {
            if (LArClusterHelper::GetClusterHitType(pCluster) != TPC_VIEW_W)
                continue;

            const MCParticle* const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

            if ((m_particleContained && !m_data && IsParticleContained(pMCParticle)) || (m_particleContained && m_data && !IsRecoParticleContained(pCluster)))
                continue;

            if (pMCParticle->GetParticleId() != m_targetParticlePDG)
                continue;

            if (m_cosmic && !LArMCParticleHelper::IsCosmicRay(pMCParticle))
                continue;

            if (!m_cosmic && LArMCParticleHelper::IsCosmicRay(pMCParticle)) 
                continue;

            TrackDirectionTool::DirectionFitObject fitResult = m_pTrackDirectionTool->GetClusterDirection(pCluster);
            this->WriteToTree(pCluster, fitResult, clusterTreeName);

            for (TrackDirectionTool::HitCharge &hitCharge : fitResult.GetHitChargeVector())
                this->WriteHitToTree(hitCharge, hitTreeName, fitResult.GetMCDirection());

            if (m_drawFit)
            {
                //fitResult.DrawFit();
                //PANDORA_MONITORING_API(Pause(this->GetPandora()));
            }
        }

        catch (...)
        {
            std::cout << "Skipping cluster..." << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster* DirectionAnalysisAlgorithm::GetTargetClusterFromPFO(const ParticleFlowObject* pPfo)
{
    HitType hitType(TPC_VIEW_W);
    ClusterList clusterListW;
    //LArPfoHelper::GetTwoDClusterList(pPfo, clusterListW);
    LArPfoHelper::GetClusters(pPfo, hitType, clusterListW);

    if (clusterListW.size() == 0)
    {    
        //std::cout << "No W cluster in PFO." << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }    

    float currentLength(std::numeric_limits<float>::min());
    ClusterList longestClusterList;

    for (const Cluster *const pCluster : clusterListW)
    {    
        if (LArClusterHelper::GetLength(pCluster) > currentLength)
        {    
            currentLength = LArClusterHelper::GetLength(pCluster);
            longestClusterList.clear();
            longestClusterList.push_back(pCluster);
        }    
    }    

    const Cluster *const pCluster(*(longestClusterList.begin()));
    return pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::IsClusterTwoParticles(const Cluster *const pCluster, TrackDirectionTool::HitChargeVector forwardsFitCharges, TrackDirectionTool::HitChargeVector backwardsFitCharges, bool &isTwoParticles)
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    int nContributingParticles(0), contributionThreshold(10);
    int nSecondaryHits(0);
    float chargeContributionThreshold(0.5);
    std::map<int, int> primaryPDGToContributionMap;

    const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

    for (const CaloHit* pCaloHit : caloHitList)
    {    
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {    
            int primaryPDG(mapIter->first->GetParticleId());
            if (primaryPDG == 11 && ((pCaloHit->GetPositionVector().GetZ() > pMCParticle->GetVertex().GetZ() && pCaloHit->GetPositionVector().GetZ() < pMCParticle->GetEndpoint().GetZ()) || (pCaloHit->GetPositionVector().GetZ() < pMCParticle->GetVertex().GetZ() && pCaloHit->GetPositionVector().GetZ() > pMCParticle->GetEndpoint().GetZ())))
                primaryPDG = 13;

            float contribution(mapIter->second);

            if (primaryPDGToContributionMap.find(primaryPDG) == primaryPDGToContributionMap.end())
            {    
                if (contribution > chargeContributionThreshold)
                    primaryPDGToContributionMap[primaryPDG] = 1; 
            }    
            else 
            {    
                if (contribution > chargeContributionThreshold)
                    primaryPDGToContributionMap.at(primaryPDG)++;
            }    
        }    
    }    

    for (auto &entry : primaryPDGToContributionMap)
    {    
        if (entry.first != m_targetParticlePDG)
            nSecondaryHits += entry.second;

        if (entry.second >= contributionThreshold)
            nContributingParticles++;
    } 

    float chargeContributionThreshold2(0.25);
    int backwardsNonPrimaryHits(0), backwardsPrimaryHits(0);
    for (TrackDirectionTool::HitCharge &hitCharge : backwardsFitCharges)
    {    
        const pandora::CaloHit* pCaloHit(hitCharge.GetCaloHit());
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {    
            int primaryPDG(mapIter->first->GetParticleId());

            float contribution(mapIter->second);
            if (contribution >= chargeContributionThreshold2 && primaryPDG != m_targetParticlePDG)
                backwardsNonPrimaryHits++;
            else if (contribution >= chargeContributionThreshold2 && primaryPDG == m_targetParticlePDG)
                backwardsPrimaryHits++;
        }    
    }    

    int forwardsNonPrimaryHits(0), forwardsPrimaryHits(0);
    for (TrackDirectionTool::HitCharge &hitCharge : forwardsFitCharges)
    {    
        const pandora::CaloHit* pCaloHit(hitCharge.GetCaloHit());
        MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

        for (MCParticleWeightMap::const_iterator mapIter = mcParticleWeightMap.begin(), mapIterEnd = mcParticleWeightMap.end(); mapIter != mapIterEnd; ++mapIter)
        {    
            int primaryPDG(mapIter->first->GetParticleId());

            float contribution(mapIter->second);
            if (contribution >= chargeContributionThreshold2 && primaryPDG != m_targetParticlePDG)
                forwardsNonPrimaryHits++;
            else if (contribution >= chargeContributionThreshold2 && primaryPDG == m_targetParticlePDG)
                forwardsPrimaryHits++;
        }    
    }    

    float forwardsImpurityFraction(((float)forwardsNonPrimaryHits)/forwardsFitCharges.size()), backwardsImpurityFraction(((float)backwardsNonPrimaryHits)/backwardsFitCharges.size());

    if (nSecondaryHits >= 10 || forwardsImpurityFraction > 0.8 || backwardsImpurityFraction > 0.8)
        isTwoParticles= true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IsParticleContained(const MCParticle* pMCParticle)
{
    const CartesianVector mcVertex(pMCParticle->GetVertex());
    const CartesianVector mcEndpoint(pMCParticle->GetEndpoint());

    if (pMCParticle->GetEnergy() < 0.05)
        return false;

    /*
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.FillCaloHitList(caloHitList);

    float nHitsInGap(0.f);

    for (const CaloHit* pCaloHit : caloHitList)
    {
        if (LArGeometryHelper::IsInGap(this->GetPandora(), pCaloHit->GetPositionVector(), TPC_VIEW_W, 2.5))
            nHitsInGap += 1.0;
    }

    if (nHitsInGap/caloHitList.size() >= 0.3)
        return false;
    */

    const float eVx(256.35), eVy(233.), eVz(1036.8);
    const float xBorder(10.), yBorder(20.), zBorder(10.);

    if ((mcEndpoint.GetX() < (eVx - xBorder)) && (mcEndpoint.GetX() > xBorder) && (mcEndpoint.GetY() < (eVy / 2. - yBorder)) && (mcEndpoint.GetY() > (-eVy / 2. + yBorder)) && (mcEndpoint.GetZ() < (eVz - zBorder)) && (mcEndpoint.GetZ() > zBorder))
    {
        if (!LArGeometryHelper::IsInGap(this->GetPandora(), mcVertex, TPC_VIEW_W, 0.05*((mcEndpoint - mcVertex).GetMagnitude())) && !LArGeometryHelper::IsInGap(this->GetPandora(), mcEndpoint, TPC_VIEW_W, 0.05*((mcEndpoint - mcVertex).GetMagnitude())))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IsRecoParticleContained(const ParticleFlowObject* pPfo)
{
    const pandora::Vertex *const pVertex = LArPfoHelper::GetVertex(pPfo);
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    LArTrackStateVector trackStateVector;
    LArPfoHelper::GetSlidingFitTrajectory(pPfo, pVertex, 10, slidingFitPitch, trackStateVector);

    if (trackStateVector.size() < 2)
        return false;

    TrackState firstTrackState(*(trackStateVector.begin())), lastTrackState(trackStateVector.back());
    const pandora::CartesianVector initialPosition(firstTrackState.GetPosition());
    const pandora::CartesianVector endPosition(lastTrackState.GetPosition());

    const float trackLength3D((endPosition - initialPosition).GetMagnitude());

    const pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector lowZVector(initialPosition.GetZ() < endPosition.GetZ() ? initialPosition : endPosition);
    const pandora::CartesianVector highZVector(initialPosition.GetZ() > endPosition.GetZ() ? initialPosition : endPosition);

    const float eVx(256.35), eVy(233.), eVz(1036.8);
    const float xBorder(10.), yBorder(20.), zBorder(10.);

    if (!m_cosmic)
    {
        if ((lowZVector.GetX() < (eVx / 2. - xBorder)) && (lowZVector.GetX() > (-eVx / 2. + xBorder)) && (lowZVector.GetY() < (eVy / 2. - yBorder)) && (lowZVector.GetY() > (-eVy / 2. + yBorder)) && (lowYVector.GetZ() < (eVz - zBorder)) && (lowZVector.GetZ() > zBorder) && (highZVector.GetX() < (eVx / 2. - xBorder)) && (highZVector.GetX() > (-eVx / 2. + xBorder)) && (highZVector.GetY() < (eVy / 2. - yBorder)) && (highZVector.GetY() > (-eVy / 2. + yBorder)) && (lowYVector.GetZ() < (eVz - zBorder)) && (highZVector.GetZ() > zBorder))
        {
            if (!LArGeometryHelper::IsInGap3D(this->GetPandora(), lowZVector, TPC_VIEW_W, 0.05*trackLength3D) && !LArGeometryHelper::IsInGap3D(this->GetPandora(), highZVector, TPC_VIEW_W, 0.05*trackLength3D))
                return true;
        }
    }
    else
    { 
        if ((lowYVector.GetX() < (eVx / 2. - xBorder)) && (lowYVector.GetX() > (-eVx / 2. + xBorder)) && (lowYVector.GetY() < (eVy / 2. - yBorder)) && (lowYVector.GetY() > (-eVy / 2. + yBorder)) && (lowYVector.GetZ() < (eVz - zBorder)) && (lowYVector.GetZ() > zBorder))
        {
            if (!LArGeometryHelper::IsInGap3D(this->GetPandora(), lowYVector, TPC_VIEW_W, 0.05*trackLength3D))
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IsRecoParticleContained(const Cluster* pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFit(pCluster, 20, slidingFitPitch);

    const CartesianVector lowZVector(slidingFit.GetGlobalMinLayerPosition());
    const CartesianVector highZVector(slidingFit.GetGlobalMaxLayerPosition());

    const float eVx(256.35), eVz(1036.8);
    const float xBorder(10.), zBorder(10.);

    if ((lowZVector.GetX() < (eVx / 2. - xBorder)) && (lowZVector.GetX() > (-eVx / 2. + xBorder)) && (lowZVector.GetZ() < (eVz - zBorder)) && (lowZVector.GetZ() > zBorder) && (highZVector.GetX() < (eVx / 2. - xBorder)) && (highZVector.GetX() > (-eVx / 2. + xBorder)) && (lowZVector.GetZ() < (eVz - zBorder)) && (highZVector.GetZ() > zBorder))
    {
        if (!LArGeometryHelper::IsInGap3D(this->GetPandora(), lowZVector, TPC_VIEW_W, 0.05*(highZVector.GetZ() - lowZVector.GetZ())) && !LArGeometryHelper::IsInGap3D(this->GetPandora(), highZVector, TPC_VIEW_W, 0.05*(highZVector.GetZ() - lowZVector.GetZ())))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::IntersectsYFace(const MCParticle* pMCParticle)
{
    CartesianVector mcBeginPoint(pMCParticle->GetVertex()), mcEndpoint(pMCParticle->GetEndpoint());
    CartesianVector lowestPoint(mcBeginPoint.GetY() < mcEndpoint.GetY() ? mcBeginPoint : mcEndpoint), highestPoint(mcBeginPoint.GetY() > mcEndpoint.GetY() ? mcBeginPoint : mcEndpoint);
    float xExtent(highestPoint.GetX() - lowestPoint.GetX()), yExtent(highestPoint.GetY() - lowestPoint.GetY()), zExtent(highestPoint.GetZ() - lowestPoint.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - lowestPoint.GetY());
    CartesianVector yFaceIntersection(lowestPoint.GetX() + xSlope*yDistanceToTravel, 116.5, lowestPoint.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 128.0 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1037.0)
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DirectionAnalysisAlgorithm::RecoIntersectsYFace(TrackDirectionTool::DirectionFitObject &fitResult)
{
    const pandora::CartesianVector initialPosition(fitResult.GetBeginpoint());
    const pandora::CartesianVector endPosition(fitResult.GetEndpoint());

    const pandora::CartesianVector lowYVector(initialPosition.GetY() < endPosition.GetY() ? initialPosition : endPosition);
    const pandora::CartesianVector highYVector(initialPosition.GetY() > endPosition.GetY() ? initialPosition : endPosition);

    float xExtent(highYVector.GetX() - lowYVector.GetX()), yExtent(highYVector.GetY() - lowYVector.GetY()), zExtent(highYVector.GetZ() - lowYVector.GetZ());
    float xSlope(xExtent/yExtent), zSlope(zExtent/yExtent);
    float yDistanceToTravel(116.5 - lowYVector.GetY());
    CartesianVector yFaceIntersection(lowYVector.GetX() + xSlope*yDistanceToTravel, 116.5, lowYVector.GetZ() + zSlope*yDistanceToTravel);

    if (yFaceIntersection.GetX() > 0.0 && yFaceIntersection.GetX() < 128.0 && yFaceIntersection.GetZ() > 0.0 && yFaceIntersection.GetZ() < 1037.0)
        return true;
    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteToTree(const Cluster* pCluster, TrackDirectionTool::DirectionFitObject &fitResult, std::string &treeName)
{
    if (!m_writeToTree)
        return;

    CartesianVector xAxis(1.f, 0.f, 0.f), yAxis(0.f, 1.f, 0.f);

    if (!m_data)
    {
        const MCParticle* pMainMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

        CartesianVector mcEndpoint(pMainMCParticle->GetEndpoint());
        CartesianVector mcBeginpoint(pMainMCParticle->GetVertex());
        int neutrinoInduced(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMainMCParticle));

        int mcForwards(mcBeginpoint.GetZ() < mcEndpoint.GetZ() ? 1 : 0);
        int mcDownwards(mcBeginpoint.GetY() > mcEndpoint.GetY() ? 1 : 0);

        CartesianVector mcDirection((mcEndpoint - mcBeginpoint).GetUnitVector());
        float mcPhi(mcDirection.GetOpeningAngle(xAxis)), mcTheta(mcDirection.GetOpeningAngle(yAxis));

        int mcIntersectsYFace(this->IntersectsYFace(MCParticleHelper::GetMainMCParticle(pCluster)) ? 1 : 0);

        bool isClusterTwoParticles(false);
        this->IsClusterTwoParticles(pCluster, fitResult.GetForwardsFitCharges(), fitResult.GetBackwardsFitCharges(), isClusterTwoParticles);
        int mcIsClusterTwoParticles(isClusterTwoParticles? 1 : 0);

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCDownwards", mcDownwards));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCForwards", mcForwards)); 
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCIntersectsYFace", mcIntersectsYFace));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCPhi", mcPhi));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCTheta", mcTheta));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCIsClusterTwoParticles", mcIsClusterTwoParticles));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "NeutrinoInduced", neutrinoInduced));
    }

    CartesianVector recoDirectionVector((fitResult.GetEndpoint() - fitResult.GetBeginpoint()).GetUnitVector());
    float recoPhi(recoDirectionVector.GetOpeningAngle(xAxis)), recoTheta(recoDirectionVector.GetOpeningAngle(yAxis));

    int recoForwards(fitResult.GetEndpoint().GetZ() > fitResult.GetBeginpoint().GetZ() ? 1 : 0); //1 by definition, fill just as sanity check
    int recoDownwards(fitResult.GetEndpoint().GetY() < fitResult.GetBeginpoint().GetY() ? 1 : 0);

    float UpwardsChiSquared(recoDownwards == 0 ? fitResult.GetForwardsChiSquared() : fitResult.GetBackwardsChiSquared()), DownwardsChiSquared(recoDownwards == 1 ? fitResult.GetForwardsChiSquared() : fitResult.GetBackwardsChiSquared());
    float UpwardsChiSquaredPerHit(recoDownwards == 0 ? fitResult.GetForwardsChiSquaredPerHit() : fitResult.GetBackwardsChiSquaredPerHit()), DownwardsChiSquaredPerHit(recoDownwards == 1 ? fitResult.GetForwardsChiSquaredPerHit() : fitResult.GetBackwardsChiSquaredPerHit());
    float UpDownDeltaChiSquaredPerHit(DownwardsChiSquaredPerHit - UpwardsChiSquaredPerHit); 

    float recoLength((fitResult.GetEndpoint() - fitResult.GetBeginpoint()).GetMagnitude());
    int recoIntersectsYFace(this->RecoIntersectsYFace(fitResult) ? 1 : 0);

    //------------------------------------------------------------

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoDownwards", recoDownwards));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoForwards", recoForwards)); 
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoIntersectsYFace", recoIntersectsYFace));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoPhi", recoPhi));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "RecoTheta", recoTheta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "recoLength", recoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BeginX", fitResult.GetBeginpoint().GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BeginY", fitResult.GetBeginpoint().GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BeginZ", fitResult.GetBeginpoint().GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EndX", fitResult.GetEndpoint().GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EndY", fitResult.GetEndpoint().GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EndZ", fitResult.GetEndpoint().GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "Probability", fitResult.GetProbability()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "Hypothesis", fitResult.GetHypothesis()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MinChiSquaredPerHit", fitResult.GetMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DeltaChiSquaredPerHit", fitResult.GetDeltaChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsChiSquared", fitResult.GetForwardsChiSquared()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsChiSquared", fitResult.GetBackwardsChiSquared()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsChiSquaredPerHit", fitResult.GetForwardsChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsChiSquaredPerHit", fitResult.GetBackwardsChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpDownDeltaChiSquaredPerHit", UpDownDeltaChiSquaredPerHit));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpwardsChiSquared", UpwardsChiSquared));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DownwardsChiSquared", DownwardsChiSquared));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "UpwardsChiSquaredPerHit", UpwardsChiSquaredPerHit));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "DownwardsChiSquaredPerHit", DownwardsChiSquaredPerHit));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "NumberHits", fitResult.GetNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "SplittingBeforeChiSquaredPerHit", fitResult.GetSplitObject().GetBeforeMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "SplittingAfterChiSquaredPerHit", fitResult.GetSplitObject().GetAfterMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "SplittingBeforeNumberHits", fitResult.GetSplitObject().GetBeforeNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "SplittingAfterNumberHits", fitResult.GetSplitObject().GetAfterNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TEFBeforeChiSquaredPerHit", fitResult.GetTEFObject().GetBeforeMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TEFAfterChiSquaredPerHit", fitResult.GetTEFObject().GetAfterMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TEFBeforeNumberHits", fitResult.GetSplitObject().GetBeforeNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "TEFAfterNumberHits", fitResult.GetSplitObject().GetAfterNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FRBeforeChiSquaredPerHit", fitResult.GetFRObject().GetBeforeMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FRAfterChiSquaredPerHit", fitResult.GetFRObject().GetAfterMinChiSquaredPerHit()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FRBeforeNumberHits", fitResult.GetSplitObject().GetBeforeNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FRAfterNumberHits", fitResult.GetSplitObject().GetAfterNHits()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "FileIdentifier", m_fileIdentifier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "EventNumber", m_eventNumber));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DirectionAnalysisAlgorithm::WriteHitToTree(TrackDirectionTool::HitCharge &hitCharge, std::string &treeName, int mcDirection)
{
    if (!m_writeToTree)
        return;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "MCDirection", mcDirection));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ChargeOverWidth", hitCharge.GetChargeOverWidth()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "HitCharge", hitCharge.GetCharge()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "HitWidth", hitCharge.GetHitWidth()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsFitCharge", hitCharge.GetForwardsFitCharge()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsSigma", hitCharge.GetForwardsSigma()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsDelta", hitCharge.GetForwardsDelta()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "ForwardsChiSquared", hitCharge.GetForwardsChiSquared()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsFitCharge", hitCharge.GetBackwardsFitCharge()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsSigma", hitCharge.GetBackwardsSigma()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsDelta", hitCharge.GetBackwardsDelta()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName.c_str(), "BackwardsChiSquared", hitCharge.GetBackwardsChiSquared()));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DirectionAnalysisAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputHitListName", m_inputHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TargetParticlePDG", m_targetParticlePDG));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ParticleContained", m_particleContained));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Cosmic", m_cosmic));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Data", m_data));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DrawFit", m_drawFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FileIdentifier", m_fileIdentifier));

    if (m_writeToTree)
    {    
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    } 

    AlgorithmTool *pAlgorithmTool(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackDirection", pAlgorithmTool));

    if (!(this->m_pTrackDirectionTool = dynamic_cast<TrackDirectionTool *>(pAlgorithmTool)))
        throw STATUS_CODE_FAILURE;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

