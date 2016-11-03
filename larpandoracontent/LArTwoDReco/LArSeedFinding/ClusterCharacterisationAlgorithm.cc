/**
 *  @file   larpandoracontent/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm() :
    m_slidingFitWindow(10),
    m_minHitsInCluster(20),
    m_maxLayerGapFraction(0.2f),
    m_maxWidthPerUnitLength(0.15f),
    m_maxShowerLength(1000.f),
    m_useDetectorGaps(true),
    m_overwriteExistingId(false),
    m_useUnavailableClusters(false),
    m_writeToTree(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterCharacterisationAlgorithm::~ClusterCharacterisationAlgorithm()
{
    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    m_clusterDirectionMap.clear();

    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ClusterCharacterisationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!m_overwriteExistingId && (UNKNOWN_PARTICLE_TYPE != pCluster->GetParticleId()))
                continue;

            if (!m_useUnavailableClusters && !PandoraContentApi::IsAvailable(*this, pCluster))
                continue;

            PandoraContentApi::Cluster::Metadata metadata;

            if (this->IsClearTrack(pCluster, pClusterList))
            {
                metadata.m_particleId = MU_MINUS;
            }
            else
            {
                metadata.m_particleId = E_MINUS;
            }

            if (pCluster->GetParticleId() != metadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
        }
    }

    m_clusterDirectionMap.clear();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster, const ClusterList *const pClusterList) const
{
    bool isTrueTrack(false);

    try
    {
        // ATTN Slightly curious definition of a clear track, but this is most-likely what is needed for shower-growing
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
        isTrueTrack = (PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId()));
    }
    catch (StatusCodeException &)
    {
    }

    // Tree variables here
    //--------------------------------------------------------------------------------------------------------------------------------------
    const int trueTrack(isTrueTrack ? 1 : 0);
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueTrack", trueTrack));

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const int hitTypeInt(static_cast<int>(hitType));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "hitType", hitTypeInt));

    const int nHits(pCluster->GetNCaloHits());
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits", nHits));

    // vertexDistance
    float vertexDistance(-1.f);
    const VertexList *pVertexList = nullptr;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pSelectedVertex((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

    if (pSelectedVertex)
    {
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), LArClusterHelper::GetClusterHitType(pCluster)));
        vertexDistance = LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertexDistance", vertexDistance));

    // straight line length and integrated pathlength
    float straightLineLength(-1.f), integratedPathLength(-1.f);

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

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength", straightLineLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "integratedPathLength", integratedPathLength));

    const float widthDirectionX = widthSet ? maxXDir - minXDir : -1.f;
    const float widthDirectionZ = widthSet ? maxZDir - minZDir : -1.f;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "widthDirectionX", widthDirectionX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "widthDirectionZ", widthDirectionZ));

    const float rLWidth = slidingFitSuccess ? rLMax - rLMin : -1.f;
    const float rTWidth = slidingFitSuccess ? rTMax - rTMin : -1.f;
    const float dTdLWidth = slidingFitSuccess ? dTdLMax - dTdLMin : -1.f;
    const float rmsWidth = slidingFitSuccess ? rmsMax - rmsMin : -1.f;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rLWidth", rLWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rTWidth", rTWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLWidth", dTdLWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rmsWidth", rmsWidth));

    float rLMean(0.f), rTMean(0.f), dTdLMean(0.f), rmsMean(0.f);
    for (const float value : rLVector) rLMean += value;
    for (const float value : rTVector) rTMean += value;
    for (const float value : dTdLVector) dTdLMean += value;
    for (const float value : rmsVector) rmsMean += value;

    rLMean = !rLVector.empty() ? std::fabs(rLMean) / static_cast<float>(rLVector.size()) : -1.f;
    rTMean = !rTVector.empty() ? std::fabs(rTMean) / static_cast<float>(rTVector.size()) : -1.f;
    dTdLMean = !dTdLVector.empty() ? std::fabs(dTdLMean) / static_cast<float>(dTdLVector.size()) : -1.f;
    rmsMean = !rmsVector.empty() ? std::fabs(rmsMean) / static_cast<float>(rmsVector.size()) : -1.f;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rLMean", rLMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rTMean", rTMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLMean", dTdLMean));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rmsMean", rmsMean));

    float rLSigma(0.f), rTSigma(0.f), dTdLSigma(0.f), rmsSigma(0.f);
    for (const float value : rLVector) rLSigma += (value - rLMean) * (value - rLMean);
    for (const float value : rTVector) rTSigma += (value - rTMean) * (value - rTMean);
    for (const float value : dTdLVector) dTdLSigma += (value - dTdLMean) * (value - dTdLMean);
    for (const float value : rmsVector) rmsSigma += (value - rmsMean) * (value - rmsMean);

    rLSigma = !rLVector.empty() ? std::sqrt(rLSigma / static_cast<float>(rLVector.size())) : -1.f;
    rTSigma = !rTVector.empty() ? std::sqrt(rTSigma / static_cast<float>(rTVector.size())) : -1.f;
    dTdLSigma = !dTdLVector.empty() ? std::sqrt(dTdLSigma / static_cast<float>(dTdLVector.size())) : -1.f;
    rmsSigma = !rmsVector.empty() ? std::sqrt(rmsSigma / static_cast<float>(rmsVector.size())) : -1.f;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rLSigma", rLSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rTSigma", rTSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLSigma", dTdLSigma));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "rmsSigma", rmsSigma));

    // hit positions and energy
    FloatVector xPositions, zPositions;
    float mipEnergy(0.f);
    int nHitsOutsideLimits(0);

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
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

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xPositions", &xPositions));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zPositions", &zPositions));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mipEnergy", mipEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsOutsideLimits", nHitsOutsideLimits));

    // shower fit width and gap length
    float showerFitWidth(-1.f), showerFitGapLength(-1.f);

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

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitWidth", showerFitWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitGapLength", showerFitGapLength));

    // nPointsOfContact, nHitsInPrimaryBranches
    int nPointsOfContact(0), nHitsInBranches(0);
    const CartesianVector vertexPosition2D(pSelectedVertex ? LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), hitType) : CartesianVector(0.f, 0.f, 0.f));

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

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPointsOfContact", nPointsOfContact));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBranches", nHitsInBranches));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    //--------------------------------------------------------------------------------------------------------------------------------------
    // End tree variables calculations

    if (nHits < 6)
        return false;

    if (straightLineLength < std::numeric_limits<float>::epsilon())
        return false;

    if (straightLineLength > 80.f)
        return true;

    if (vertexDistance / straightLineLength > 0.5f)
        return false;

    if (showerFitWidth < 0.f || showerFitWidth / straightLineLength > 0.35f)
        return false;

    if (rTWidth < 0.f || rTWidth / straightLineLength > 0.05f)
        return false;

    if (integratedPathLength < 0.f || integratedPathLength / straightLineLength > 1.005f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsInCluster", m_minHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLayerGapFraction", m_maxLayerGapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxWidthPerUnitLength", m_maxWidthPerUnitLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerLength", m_maxShowerLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDetectorGaps", m_useDetectorGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverwriteExistingId", m_overwriteExistingId));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnavailableClusters", m_useUnavailableClusters));

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
