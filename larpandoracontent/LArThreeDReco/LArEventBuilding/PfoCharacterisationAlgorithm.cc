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

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/PfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoCharacterisationAlgorithm::PfoCharacterisationAlgorithm() :
    m_updateClusterIds(true),
    m_postBranchAddition(false),
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
    m_clusterDirectionMap.clear();
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
            PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;

            if (this->IsClearTrack(pPfo))
            {
                pfoMetadata.m_particleId = MU_MINUS;

                if (m_showerPfoListName == pfoListName)
                    showersToTracks.push_back(pPfo);
            }
            else
            {
                pfoMetadata.m_particleId = E_MINUS;

                if (m_trackPfoListName == pfoListName)
                    tracksToShowers.push_back(pPfo);
            }

            if (pPfo->GetParticleId() != pfoMetadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, pfoMetadata));

            if (!m_updateClusterIds)
                continue;

            ClusterList twoDClusterList;
            LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

            for (const Cluster *const pCluster : twoDClusterList)
            {
                if (pCluster->GetParticleId() == pfoMetadata.m_particleId.Get())
                    continue;

                PandoraContentApi::Cluster::Metadata clusterMetadata;
                clusterMetadata.m_particleId = pfoMetadata.m_particleId.Get();
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, clusterMetadata));
            }
        }
    }

    if (!tracksToShowers.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

    if (!showersToTracks.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));

    m_clusterDirectionMap.clear();
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
    float bestMCParticleLength(-1.f);
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
        const int absParticleId(std::abs(pBestMCParticle->GetParticleId()));
        isTrueTrack = ((MU_MINUS == absParticleId) || (PROTON == absParticleId) || (PI_PLUS == absParticleId));
        bestMCParticleLength = (pBestMCParticle->GetEndpoint() - pBestMCParticle->GetVertex()).GetMagnitude();
    }

    // Tree variables here
    //--------------------------------------------------------------------------------------------------------------------------------------
    const int trueTrack(isTrueTrack ? 1 : 0);
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueTrack", trueTrack));

    const int inputParticleId(pPfo->GetParticleId());
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "inputParticleId", inputParticleId));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticleLength", bestMCParticleLength));

    int nHitsU(0), nHitsV(0), nHitsW(0);
    int nClustersU(0), nClustersV(0), nClustersW(0);
    float vertexDistanceU(-1.f), vertexDistanceV(-1.f), vertexDistanceW(-1.f);

    FloatVector xPositionsU, xPositionsV, xPositionsW, zPositionsU, zPositionsV, zPositionsW;
    float mipEnergyU(0.f), mipEnergyV(0.f), mipEnergyW(0.f);
    int nHitsOutsideLimitsU(0), nHitsOutsideLimitsV(0), nHitsOutsideLimitsW(0);

    float straightLineLengthU(-1.f), straightLineLengthV(-1.f), straightLineLengthW(-1.f);
    float integratedPathLengthU(-1.f), integratedPathLengthV(-1.f), integratedPathLengthW(-1.f);

    float widthDirectionXU(-1.f), widthDirectionXV(-1.f), widthDirectionXW(-1.f);
    float widthDirectionZU(-1.f), widthDirectionZV(-1.f), widthDirectionZW(-1.f);

    float rTWidthU(-1.f), rTWidthV(-1.f), rTWidthW(-1.f);
    float rLWidthU(-1.f), rLWidthV(-1.f), rLWidthW(-1.f);
    float dTdLWidthU(-1.f), dTdLWidthV(-1.f), dTdLWidthW(-1.f);
    float rmsWidthU(-1.f), rmsWidthV(-1.f), rmsWidthW(-1.f);

    float rTMeanU(0.f), rTMeanV(0.f), rTMeanW(0.f);
    float rLMeanU(0.f), rLMeanV(0.f), rLMeanW(0.f);
    float dTdLMeanU(0.f), dTdLMeanV(0.f), dTdLMeanW(0.f);
    float rmsMeanU(0.f), rmsMeanV(0.f), rmsMeanW(0.f);

    float rTSigmaU(0.f), rTSigmaV(0.f), rTSigmaW(0.f);
    float rLSigmaU(0.f), rLSigmaV(0.f), rLSigmaW(0.f);
    float dTdLSigmaU(0.f), dTdLSigmaV(0.f), dTdLSigmaW(0.f);
    float rmsSigmaU(0.f), rmsSigmaV(0.f), rmsSigmaW(0.f);

    float showerFitWidthU(-1.f), showerFitWidthV(-1.f), showerFitWidthW(-1.f);
    float showerFitGapLengthU(-1.f), showerFitGapLengthV(-1.f), showerFitGapLengthW(-1.f);

    int nPointsOfContactU(0), nPointsOfContactV(0), nPointsOfContactW(0);
    int nHitsInBranchesU(0), nHitsInBranchesV(0), nHitsInBranchesW(0);

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

        // straight line length and integrated pathlength
        float &straightLineLength((TPC_VIEW_U == hitType) ? straightLineLengthU : (TPC_VIEW_V == hitType) ? straightLineLengthV : straightLineLengthW);
        float &integratedPathLength((TPC_VIEW_U == hitType) ? integratedPathLengthU : (TPC_VIEW_V == hitType) ? integratedPathLengthV : integratedPathLengthW);

        bool slidingFitSuccess(false);
        float xMin(0.f), xMax(0.f), zMin(0.f), zMax(0.f);

        bool widthSet(false);
        float minXDir(+std::numeric_limits<float>::max()), minZDir(+std::numeric_limits<float>::max());
        float maxXDir(-std::numeric_limits<float>::max()), maxZDir(-std::numeric_limits<float>::max());
        float rLMin(+std::numeric_limits<float>::max()), rTMin(+std::numeric_limits<float>::max()), dTdLMin(+std::numeric_limits<float>::max()), rmsMin(+std::numeric_limits<float>::max());
        float rLMax(-std::numeric_limits<float>::max()), rTMax(-std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max()), rmsMax(-std::numeric_limits<float>::max());
        FloatVector rLVector, rTVector, dTdLVector, rmsVector;

        try
        {
            const TwoDSlidingFitResult slidingFitResult(pCluster, 5, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
            const CartesianVector globalMaxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());
            straightLineLength = (globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude();

            slidingFitSuccess = true;
            xMin = std::min(globalMinLayerPosition.GetX(), globalMaxLayerPosition.GetX());
            xMax = std::max(globalMinLayerPosition.GetX(), globalMaxLayerPosition.GetX());
            zMin = std::min(globalMinLayerPosition.GetZ(), globalMaxLayerPosition.GetZ());
            zMax = std::max(globalMinLayerPosition.GetZ(), globalMaxLayerPosition.GetZ());

            integratedPathLength = 0.f;
            CartesianVector previousFitPosition(globalMinLayerPosition);
            const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

            for (const auto &mapEntry : layerFitResultMap)
            {
                rLMin = std::min(rLMin, static_cast<float>(mapEntry.second.GetL()));
                rLMax = std::max(rLMax, static_cast<float>(mapEntry.second.GetL()));
                rTMin = std::min(rTMin, static_cast<float>(mapEntry.second.GetFitT()));
                rTMax = std::max(rTMax, static_cast<float>(mapEntry.second.GetFitT()));

                dTdLMin = std::min(dTdLMin, static_cast<float>(mapEntry.second.GetGradient()));
                dTdLMax = std::max(dTdLMax, static_cast<float>(mapEntry.second.GetGradient()));
                rmsMin = std::min(rmsMin, static_cast<float>(mapEntry.second.GetRms()));
                rmsMax = std::max(rmsMax, static_cast<float>(mapEntry.second.GetRms()));

                rTVector.push_back(mapEntry.second.GetFitT());
                rLVector.push_back(mapEntry.second.GetL());
                dTdLVector.push_back(mapEntry.second.GetGradient());
                rmsVector.push_back(mapEntry.second.GetRms());

                CartesianVector thisFitPosition(0.f, 0.f, 0.f);
                slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
                integratedPathLength += (thisFitPosition - previousFitPosition).GetMagnitude();
                previousFitPosition = thisFitPosition;

                CartesianVector thisFitDirection(0.f, 0.f, 0.f);
                if (STATUS_CODE_SUCCESS == slidingFitResult.GetGlobalFitDirection(mapEntry.second.GetL(), thisFitDirection))
                {
                    minXDir = std::min(minXDir, thisFitDirection.GetX());
                    maxXDir = std::max(maxXDir, thisFitDirection.GetX());
                    minZDir = std::min(minZDir, thisFitDirection.GetZ());
                    maxZDir = std::max(maxZDir, thisFitDirection.GetZ());
                    widthSet = true;
                }
            }
        }
        catch (const StatusCodeException &)
        {
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength" + hitTypeLabel, straightLineLength));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "integratedPathLength" + hitTypeLabel, integratedPathLength));

        float &widthDirectionX((TPC_VIEW_U == hitType) ? widthDirectionXU : (TPC_VIEW_V == hitType) ? widthDirectionXV : widthDirectionXW);
        float &widthDirectionZ((TPC_VIEW_U == hitType) ? widthDirectionZU : (TPC_VIEW_V == hitType) ? widthDirectionZV : widthDirectionZW);

        widthDirectionX = widthSet ? maxXDir - minXDir : -1.f;
        widthDirectionZ = widthSet ? maxZDir - minZDir : -1.f;

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "widthDirectionX" + hitTypeLabel, widthDirectionX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "widthDirectionZ" + hitTypeLabel, widthDirectionZ));

        float &rLWidth((TPC_VIEW_U == hitType) ? rLWidthU : (TPC_VIEW_V == hitType) ? rLWidthV : rLWidthW);
        float &rTWidth((TPC_VIEW_U == hitType) ? rTWidthU : (TPC_VIEW_V == hitType) ? rTWidthV : rTWidthW);
        float &dTdLWidth((TPC_VIEW_U == hitType) ? dTdLWidthU : (TPC_VIEW_V == hitType) ? dTdLWidthV : dTdLWidthW);
        float &rmsWidth((TPC_VIEW_U == hitType) ? rmsWidthU : (TPC_VIEW_V == hitType) ? rmsWidthV : rmsWidthW);

        rLWidth = slidingFitSuccess ? rLMax - rLMin : -1.f;
        rTWidth = slidingFitSuccess ? rTMax - rTMin : -1.f;
        dTdLWidth = slidingFitSuccess ? dTdLMax - dTdLMin : -1.f;
        rmsWidth = slidingFitSuccess ? rmsMax - rmsMin : -1.f;

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rLWidth" + hitTypeLabel, rLWidth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rTWidth" + hitTypeLabel, rTWidth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLWidth" + hitTypeLabel, dTdLWidth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rmsWidth" + hitTypeLabel, rmsWidth));

        float &rLMean((TPC_VIEW_U == hitType) ? rLMeanU : (TPC_VIEW_V == hitType) ? rLMeanV : rLMeanW);
        float &rTMean((TPC_VIEW_U == hitType) ? rTMeanU : (TPC_VIEW_V == hitType) ? rTMeanV : rTMeanW);
        float &dTdLMean((TPC_VIEW_U == hitType) ? dTdLMeanU : (TPC_VIEW_V == hitType) ? dTdLMeanV : dTdLMeanW);
        float &rmsMean((TPC_VIEW_U == hitType) ? rmsMeanU : (TPC_VIEW_V == hitType) ? rmsMeanV : rmsMeanW);

        for (const float value : rLVector) rLMean += value;
        for (const float value : rTVector) rTMean += value;
        for (const float value : dTdLVector) dTdLMean += value;
        for (const float value : rmsVector) rmsMean += value;

        rLMean = !rLVector.empty() ? std::fabs(rLMean) / static_cast<float>(rLVector.size()) : -1.f;
        rTMean = !rTVector.empty() ? std::fabs(rTMean) / static_cast<float>(rTVector.size()) : -1.f;
        dTdLMean = !dTdLVector.empty() ? std::fabs(dTdLMean) / static_cast<float>(dTdLVector.size()) : -1.f;
        rmsMean = !rmsVector.empty() ? std::fabs(rmsMean) / static_cast<float>(rmsVector.size()) : -1.f;

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rLMean" + hitTypeLabel, rLMean));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rTMean" + hitTypeLabel, rTMean));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLMean" + hitTypeLabel, dTdLMean));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rmsMean" + hitTypeLabel, rmsMean));

        float &rLSigma((TPC_VIEW_U == hitType) ? rLSigmaU : (TPC_VIEW_V == hitType) ? rLSigmaV : rLSigmaW);
        float &rTSigma((TPC_VIEW_U == hitType) ? rTSigmaU : (TPC_VIEW_V == hitType) ? rTSigmaV : rTSigmaW);
        float &dTdLSigma((TPC_VIEW_U == hitType) ? dTdLSigmaU : (TPC_VIEW_V == hitType) ? dTdLSigmaV : dTdLSigmaW);
        float &rmsSigma((TPC_VIEW_U == hitType) ? rmsSigmaU : (TPC_VIEW_V == hitType) ? rmsSigmaV : rmsSigmaW);

        for (const float value : rLVector) rLSigma += (value - rLMean) * (value - rLMean);
        for (const float value : rTVector) rTSigma += (value - rTMean) * (value - rTMean);
        for (const float value : dTdLVector) dTdLSigma += (value - dTdLMean) * (value - dTdLMean);
        for (const float value : rmsVector) rmsSigma += (value - rmsMean) * (value - rmsMean);

        rLSigma = !rLVector.empty() ? std::sqrt(rLSigma / static_cast<float>(rLVector.size())) : -1.f;
        rTSigma = !rTVector.empty() ? std::sqrt(rTSigma / static_cast<float>(rTVector.size())) : -1.f;
        dTdLSigma = !dTdLVector.empty() ? std::sqrt(dTdLSigma / static_cast<float>(dTdLVector.size())) : -1.f;
        rmsSigma = !rmsVector.empty() ? std::sqrt(rmsSigma / static_cast<float>(rmsVector.size())) : -1.f;

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rLSigma" + hitTypeLabel, rLSigma));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rTSigma" + hitTypeLabel, rTSigma));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLSigma" + hitTypeLabel, dTdLSigma));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rmsSigma" + hitTypeLabel, rmsSigma));

        // hit positions and energy
        FloatVector &xPositions((TPC_VIEW_U == hitType) ? xPositionsU : (TPC_VIEW_V == hitType) ? xPositionsV : xPositionsW);
        FloatVector &zPositions((TPC_VIEW_U == hitType) ? zPositionsU : (TPC_VIEW_V == hitType) ? zPositionsV : zPositionsW);
        float &mipEnergy((TPC_VIEW_U == hitType) ? mipEnergyU : (TPC_VIEW_V == hitType) ? mipEnergyV : mipEnergyW);
        int &nHitsOutsideLimits((TPC_VIEW_U == hitType) ? nHitsOutsideLimitsU : (TPC_VIEW_V == hitType) ?  nHitsOutsideLimitsV : nHitsOutsideLimitsW);

        CaloHitList clusterHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHitList);

        for (const CaloHit *const pCaloHit : clusterHitList)
        {
            xPositions.push_back(pCaloHit->GetPositionVector().GetX());
            zPositions.push_back(pCaloHit->GetPositionVector().GetZ());
            mipEnergy += pCaloHit->GetMipEquivalentEnergy();

            if (slidingFitSuccess &&
               ((pCaloHit->GetPositionVector().GetX() < xMin) || (pCaloHit->GetPositionVector().GetX() > xMax) ||
                (pCaloHit->GetPositionVector().GetZ() < zMin) || (pCaloHit->GetPositionVector().GetX() > zMax)) )
            {
                ++nHitsOutsideLimits;
            }
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xPositions" + hitTypeLabel, &xPositions));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zPositions" + hitTypeLabel, &zPositions));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mipEnergy" + hitTypeLabel, mipEnergy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsOutsideLimits" + hitTypeLabel, nHitsOutsideLimits));

        // shower fit width and gap length
        float &showerFitWidth((TPC_VIEW_U == hitType) ? showerFitWidthU : (TPC_VIEW_V == hitType) ? showerFitWidthV : showerFitWidthW);
        float &showerFitGapLength((TPC_VIEW_U == hitType) ? showerFitGapLengthU : (TPC_VIEW_V == hitType) ? showerFitGapLengthV : showerFitGapLengthW);

        try
        {
            const TwoDSlidingShowerFitResult showerFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
            const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
            const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

            if (layerFitResultMapS.size() > 1)
            {
                CartesianVector globalMinLayerPositionOnAxis(0.f, 0.f, 0.f), globalMaxLayerPositionOnAxis(0.f, 0.f, 0.f);
                showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, globalMinLayerPositionOnAxis);
                showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, globalMaxLayerPositionOnAxis);
                const float straightLinePathLength((globalMaxLayerPositionOnAxis - globalMinLayerPositionOnAxis).GetMagnitude());

                if (straightLinePathLength > std::numeric_limits<float>::epsilon())
                {
                    showerFitWidth = 0.f;
                    showerFitGapLength = 0.f;
                    CartesianVector previousLayerPosition(globalMinLayerPositionOnAxis);

                    for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
                    {
                        CartesianVector thisLayerPosition(0.f, 0.f, 0.f);
                        showerFitResult.GetShowerFitResult().GetGlobalPosition(iterS->second.GetL(), 0.f, thisLayerPosition);
                        const float thisGapLength((thisLayerPosition - previousLayerPosition).GetMagnitude());

                        if (thisGapLength > showerFitGapLength)
                        {
                            const float minZ(std::min(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));
                            const float maxZ(std::max(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));

                            if ((maxZ - minZ) < std::numeric_limits<float>::epsilon())
                                throw StatusCodeException(STATUS_CODE_FAILURE);

                            const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, LArClusterHelper::GetClusterHitType(pCluster)));
                            const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));

                            if (correctedGapLength > showerFitGapLength)
                                showerFitGapLength = correctedGapLength;
                        }

                        previousLayerPosition = thisLayerPosition;

                        LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
                        LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

                        if ((layerFitResultMapP.end() == iterP) || (layerFitResultMapN.end() == iterN))
                            continue;

                        showerFitWidth += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
                    }
                }
            }
        }
        catch (StatusCodeException &)
        {
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitWidth" + hitTypeLabel, showerFitWidth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitGapLength" + hitTypeLabel, showerFitGapLength));

        // nPointsOfContact, nHitsInPrimaryBranches
        int &nPointsOfContact((TPC_VIEW_U == hitType) ? nPointsOfContactU : (TPC_VIEW_V == hitType) ? nPointsOfContactV : nPointsOfContactW);
        int &nHitsInBranches((TPC_VIEW_U == hitType) ? nHitsInBranchesU : (TPC_VIEW_V == hitType) ? nHitsInBranchesV : nHitsInBranchesW);

        const std::string &clusterListName((TPC_VIEW_U == hitType) ? m_clusterListNameU : (TPC_VIEW_V == hitType) ? m_clusterListNameV : m_clusterListNameW);
        const CartesianVector vertexPosition2D(pSelectedVertex ? LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), hitType) : CartesianVector(0.f, 0.f, 0.f));

        const ClusterList *pClusterList(nullptr);
        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, clusterListName, pClusterList))
        {
            ClusterVector candidateClusters;
            for (const Cluster *const pCandidateCluster : *pClusterList)
            {
                if ((pCandidateCluster == pCluster) || (pCandidateCluster->GetNCaloHits() < 5))
                    continue;

                try
                {
                    if (pSelectedVertex && this->IsVertexAssociated(LArPointingCluster(pCluster), vertexPosition2D))
                        continue;
                }
                catch (const StatusCodeException &)
                {
                }

                candidateClusters.push_back(pCandidateCluster);
            }
            std::sort(candidateClusters.begin(), candidateClusters.end(), ShowerGrowingAlgorithm::SortClusters);

            ClusterUsageMap forwardUsageMap, backwardUsageMap;
            this->FindAssociatedClusters(pCluster, candidateClusters, forwardUsageMap, backwardUsageMap);

            SeedAssociationList seedAssociationList;
            this->IdentifyClusterMerges(ClusterVector(1, pCluster), backwardUsageMap, seedAssociationList);

            for (const Cluster *const pBranchCluster : seedAssociationList.at(pCluster))
            {
                ++nPointsOfContact;
                nHitsInBranches += pBranchCluster->GetNCaloHits();
            }
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPointsOfContact" + hitTypeLabel, nPointsOfContact));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBranches" + hitTypeLabel, nHitsInBranches));
    }

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    //--------------------------------------------------------------------------------------------------------------------------------------
    // End tree variables calculations

    if (straightLineLengthW < std::numeric_limits<float>::epsilon())
        return false;

    if (m_postBranchAddition)
    {
        if (straightLineLengthW > 80.f)
            return true;

        if (showerFitWidthW < 0.f || showerFitWidthW / straightLineLengthW > 0.3f)
            return false;

        if (vertexDistanceW / straightLineLengthW > 1.0f)
            return false;

        if (dTdLWidthW < 0.f || dTdLWidthW / straightLineLengthW > 0.03f)
            return false;
    }
    else
    {
        if (straightLineLengthW > 40.f && nPointsOfContactW < 4)
            return true;

        if (showerFitWidthW < 0.f || showerFitWidthW / straightLineLengthW > 0.2f)
            return false;

        if (vertexDistanceW / straightLineLengthW > 0.6f)
            return false;

        if (dTdLWidthW < 0.f || dTdLWidthW / straightLineLengthW > 0.04f)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    m_inputPfoListNames.push_back(m_trackPfoListName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    m_inputPfoListNames.push_back(m_showerPfoListName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UpdateClusterIds", m_updateClusterIds));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PostBranchAddition", m_postBranchAddition));

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
